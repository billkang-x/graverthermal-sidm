#include "gadgetconfig.h"
#ifdef SIDM

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>
#include <mpi.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "sidm.h"
#include "../data/allvars.h"
#include "../data/simparticles.h"
#include "../system/system.h"

sidm Sidm;

static gsl_rng *sidm_rng = NULL;
static int sidm_initialized = 0;

// Linked-cell grid (每步重建一次，所有粒子共用)
static const int NGRID = 256;
static const double CELL_SIZE = 2.0;  // kpc
static const double BOXSIZE = 512.0;
static const double BOXHALF = 256.0;
static std::vector<int> grid_cells[NGRID*NGRID*NGRID];
static double grid_pos[3][81000];  // 缓存位置，避免重复intpos_to_pos

sidm::sidm()
{
  SigmaOverMass = SIDM_SIGMA_OVER_MASS;
  KNeighbors = SIDM_K_NEIGHBORS;
  PMax = SIDM_PMAX;
}

sidm::~sidm()
{
  if(sidm_rng) gsl_rng_free(sidm_rng);
}

static void sidm_init_rng(void)
{
  if(sidm_initialized) return;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  sidm_rng = gsl_rng_alloc(gsl_rng_ranlxd2);
  gsl_rng_set(sidm_rng, 42 + rank * 137);
  sidm_initialized = 1;
}

int sidm::poisson_sample(double lambda)
{
  if(lambda <= 0) return 0;
  if(lambda < 12) return (int)gsl_ran_poisson(sidm_rng, lambda);
  return (int)(lambda + gsl_ran_gaussian(sidm_rng, sqrt(lambda)) + 0.5);
}

void sidm::random_unit_vector(double *v)
{
  double z = 2.0 * gsl_rng_uniform(sidm_rng) - 1.0;
  double phi = 2.0 * M_PI * gsl_rng_uniform(sidm_rng);
  double s = sqrt(std::max(0.0, 1.0 - z*z));
  v[0] = s * cos(phi);
  v[1] = s * sin(phi);
  v[2] = z;
}

/*
 * 每步构建一次linked-cell网格（关键修复：只建一次！）
 */
static void build_grid(simparticles *Sp, int numpart)
{
  // 清空网格
  for(int c = 0; c < NGRID*NGRID*NGRID; c++) grid_cells[c].clear();

  double pos[3];
  for(int p = 0; p < numpart; p++) {
    if(Sp->P[p].getType() != 1) continue;  // 只放暗物质粒子
    Sp->intpos_to_pos(Sp->P[p].IntPos, pos);

    // 缓存位置
    grid_pos[0][p] = pos[0];
    grid_pos[1][p] = pos[1];
    grid_pos[2][p] = pos[2];

    int ix = (int)((pos[0] + BOXHALF) / CELL_SIZE);
    int iy = (int)((pos[1] + BOXHALF) / CELL_SIZE);
    int iz = (int)((pos[2] + BOXHALF) / CELL_SIZE);
    ix = std::max(0, std::min(NGRID-1, ix));
    iy = std::max(0, std::min(NGRID-1, iy));
    iz = std::max(0, std::min(NGRID-1, iz));
    grid_cells[ix*NGRID*NGRID + iy*NGRID + iz].push_back(p);
  }
}

/*
 * 从已构建的网格中找邻居
 */
static int find_neighbors_from_grid(int i, int *ngb_indices, double *ngb_dist2, int k)
{
  double xi = grid_pos[0][i];
  double yi = grid_pos[1][i];
  double zi = grid_pos[2][i];

  int ix = (int)((xi + BOXHALF) / CELL_SIZE);
  int iy = (int)((yi + BOXHALF) / CELL_SIZE);
  int iz = (int)((zi + BOXHALF) / CELL_SIZE);
  ix = std::max(0, std::min(NGRID-1, ix));
  iy = std::max(0, std::min(NGRID-1, iy));
  iz = std::max(0, std::min(NGRID-1, iz));

  std::vector<std::pair<double, int>> dms;
  dms.reserve(k * 10);

  for(int dx = -1; dx <= 1; dx++)
    for(int dy = -1; dy <= 1; dy++)
      for(int dz = -1; dz <= 1; dz++) {
        int cx = ix+dx, cy = iy+dy, cz = iz+dz;
        if(cx<0||cx>=NGRID||cy<0||cy>=NGRID||cz<0||cz>=NGRID) continue;
        int cell = cx*NGRID*NGRID + cy*NGRID + cz;
        for(int j : grid_cells[cell]) {
          if(j == i) continue;
          double dxj = xi - grid_pos[0][j];
          double dyj = yi - grid_pos[1][j];
          double dzj = zi - grid_pos[2][j];
          double d2 = dxj*dxj + dyj*dyj + dzj*dzj;
          dms.push_back(std::make_pair(d2, j));
        }
      }

  if((int)dms.size() < k) k = dms.size();
  if(k <= 0) return 0;

  std::nth_element(dms.begin(), dms.begin() + k, dms.end());

  int found = 0;
  for(int idx = 0; idx < k; idx++) {
    ngb_indices[idx] = dms[idx].second;
    ngb_dist2[idx] = dms[idx].first;
    found++;
  }
  return found;
}

void sidm::scatter_pair(simparticles *Sp, int i, int j)
{
  double mi = Sp->P[i].getMass();
  double mj = Sp->P[j].getMass();
  double M = mi + mj;
  if(M <= 0) return;

  double dv[3] = {
    Sp->P[i].Vel[0] - Sp->P[j].Vel[0],
    Sp->P[i].Vel[1] - Sp->P[j].Vel[1],
    Sp->P[i].Vel[2] - Sp->P[j].Vel[2]
  };
  double vrel2 = dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2];
  double vrel = sqrt(vrel2);
  if(vrel < 1e-6) return;

  double vcm[3] = {
    (mi*Sp->P[i].Vel[0] + mj*Sp->P[j].Vel[0])/M,
    (mi*Sp->P[i].Vel[1] + mj*Sp->P[j].Vel[1])/M,
    (mi*Sp->P[i].Vel[2] + mj*Sp->P[j].Vel[2])/M
  };

  double n[3];
  random_unit_vector(n);

  Sp->P[i].Vel[0] = vcm[0] + (mj/M) * vrel * n[0];
  Sp->P[i].Vel[1] = vcm[1] + (mj/M) * vrel * n[1];
  Sp->P[i].Vel[2] = vcm[2] + (mj/M) * vrel * n[2];

  Sp->P[j].Vel[0] = vcm[0] - (mi/M) * vrel * n[0];
  Sp->P[j].Vel[1] = vcm[1] - (mi/M) * vrel * n[1];
  Sp->P[j].Vel[2] = vcm[2] - (mi/M) * vrel * n[2];
}

void sidm::do_sidm_scattering(simparticles *Sp, double dt)
{
  if(dt <= 0) return;
  sidm_init_rng();

  int numpart = Sp->NumPart;
  if(numpart <= 0) return;

  // 每步构建一次网格（O(N)）
  build_grid(Sp, numpart);

  std::vector<char> scattered(numpart, 0);

  double L_unit = All.UnitLength_in_cm;
  double M_unit = All.UnitMass_in_g;
  double sigma_m = SigmaOverMass * M_unit / (L_unit * L_unit);

  for(int i = 0; i < numpart; i++) {
    if(Sp->P[i].getType() != 1) continue;
    if(scattered[i]) continue;

    int ngb_idx[64];
    double ngb_dist2[64];
    int nngb = find_neighbors_from_grid(i, ngb_idx, ngb_dist2, KNeighbors);
    if(nngb < KNeighbors) continue;

    double r_k = sqrt(ngb_dist2[nngb-1]);
    double volume = (4.0/3.0) * M_PI * r_k * r_k * r_k;
    if(volume <= 0) continue;

    double rates[64];
    double total_rate = 0.0;

    for(int k = 0; k < nngb; k++) {
      int j = ngb_idx[k];
      if(scattered[j]) {
        rates[k] = 0.0;
        continue;
      }
      double mj = Sp->P[j].getMass();

      double dv0 = Sp->P[i].Vel[0] - Sp->P[j].Vel[0];
      double dv1 = Sp->P[i].Vel[1] - Sp->P[j].Vel[1];
      double dv2 = Sp->P[i].Vel[2] - Sp->P[j].Vel[2];
      double vrel = sqrt(dv0*dv0 + dv1*dv1 + dv2*dv2);

      rates[k] = 0.5 * sigma_m * mj * vrel / volume;
      total_rate += rates[k];
    }

    double lambda = total_rate * dt;
    if(lambda > PMax) lambda = PMax;
    if(lambda <= 0) continue;

    int n_scatter = poisson_sample(lambda);
    if(n_scatter <= 0) continue;

    for(int s = 0; s < n_scatter; s++) {
      double r = gsl_rng_uniform(sidm_rng) * total_rate;
      double cum = 0.0;
      int selected = -1;
      for(int k = 0; k < nngb; k++) {
        cum += rates[k];
        if(cum >= r) {
          selected = ngb_idx[k];
          break;
        }
      }
      if(selected < 0 || selected == i || scattered[selected]) continue;

      scatter_pair(Sp, i, selected);
      scattered[i] = 1;
      scattered[selected] = 1;
      break;
    }
  }
}

#endif
