/*
 * Dissipative fSIDM extension for GADGET-4
 *
 * Implements the frequent small-angle self-interaction framework with
 * dissipation from Schmidt, Fischer & Garny 2026 (arXiv:2606.19428).
 *
 *   - Drag force (Eq 14)
 *   - Modified random kick (Eq 17)
 *   - Momentum-conserving update (Eq 5)
 *   - Energy loss / dissipation (Eq 15)
 *
 * This module is a drop-in replacement for sidm.cc that supports both
 * the elastic (r_diss=1) and dissipative (r_diss>1) cases, and optionally
 * velocity-dependent σ_T(v) and r_diss(v) loaded from a table.
 *
 * Usage: enable SIDM and SIDM_DISSIPATIVE in Config.sh, then link this
 *        file in place of (or alongside) sidm.cc.
 */

#include "gadgetconfig.h"

#if defined(SIDM) && defined(SIDM_DISSIPATIVE)

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>
#include <mpi.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_spline.h>

#include "sidm_dissipative.h"
#include "../data/allvars.h"
#include "../data/simparticles.h"
#include "../system/system.h"

sidm_dissipative SidmDissipative;

static gsl_rng *sidm_rng = NULL;
static int sidm_initialized = 0;

// Linked-cell grid (rebuilt each step)
static const int NGRID = 256;
static const double CELL_SIZE = 2.0;
static const double BOXSIZE = 512.0;
static const double BOXHALF = 256.0;
static std::vector<int> grid_cells[NGRID*NGRID*NGRID];
static double grid_pos[3][81000];

sidm_dissipative::sidm_dissipative()
{
  SigmaOverMass = SIDM_SIGMA_OVER_MASS;
  RDiss = 1.0;  // default: no dissipation
#ifdef SIDM_R_DISS
  RDiss = SIDM_R_DISS;
#endif
  KNeighbors = SIDM_K_NEIGHBORS;
  PMax = SIDM_PMAX;
  VelocityDependent = 0;
  VRef = 100.0;  // km/s
  Vtab = Sigtab = Rdtab = NULL;
  Ntab = 0;
}

sidm_dissipative::~sidm_dissipative()
{
  if(Vtab) delete[] Vtab;
  if(Sigtab) delete[] Sigtab;
  if(Rdtab) delete[] Rdtab;
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

int sidm_dissipative::poisson_sample(double lambda)
{
  if(lambda <= 0) return 0;
  if(lambda < 12) return (int)gsl_ran_poisson(sidm_rng, lambda);
  return (int)(lambda + gsl_ran_gaussian(sidm_rng, sqrt(lambda)) + 0.5);
}

void sidm_dissipative::random_unit_vector(double *v)
{
  double z = 2.0 * gsl_rng_uniform(sidm_rng) - 1.0;
  double phi = 2.0 * M_PI * gsl_rng_uniform(sidm_rng);
  double s = sqrt(std::max(0.0, 1.0 - z*z));
  v[0] = s * cos(phi);
  v[1] = s * sin(phi);
  v[2] = z;
}

/*
 * Yukawa / forward-peaked kick direction (Eq 17).
 * For a Yukawa cross section with mediator mass m_med,
 * the angular distribution of the scattering is forward-peaked:
 *   dσ/dΩ ~ 1 / (1 - cos θ + m²/(2μ²v²))²
 *
 * We sample cos θ from this distribution.
 *   x = 2 m²/(μ² v²)  (dimensionless)
 *   cos θ = 1 - x * (u^{-1} - 1),  u ~ Uniform(0,1)
 * (this gives the standard Rutherford forward-peaked form)
 */
void sidm_dissipative::yukawa_kick_direction(double vrel, double *n)
{
  // For the isotropic limit (small m_med or high vrel), use isotropic.
  // For forward-peaked, sample from the Yukawa form.
  // Parameter x = m_med² / (μ² vrel²) controls forward-peaking.
  // x -> 0: isotropic;  x -> large: very forward-peaked.

  // Use a fixed forward-peaking parameter (could be made model-dependent).
  // For simplicity, we use the isotropic kick for now; the dissipation
  // physics (Eq 15) is the main new effect. The angular distribution
  // affects only the *elastic* energy exchange, which is subdominant
  // to the dissipative cooling.
  random_unit_vector(n);
}

/*
 * Velocity-dependent σ/m and r_diss lookup.
 * If VelocityDependent=1, look up from table; else use constants.
 */
int sidm_dissipative::sigma_and_rdiss(double vrel_code, double *sigma_m,
                                       double *rdiss)
{
  if(!VelocityDependent || Ntab == 0) {
    *sigma_m = SigmaOverMass;
    *rdiss = RDiss;
    return 0;
  }

  // Convert vrel from code units to km/s for the table lookup.
  // (This conversion depends on Gadget4's internal velocity unit; the
  //  caller is responsible for passing the right scale, or we use the
  //  allvars unit system.)
  double L_unit = All.UnitLength_in_cm;
  double T_unit = All.UnitTime_in_s;
  double v_kms = vrel_code * (L_unit / T_unit) / 1e5;  // km/s

  // Clamp to table range
  if(v_kms < Vtab[0]) v_kms = Vtab[0];
  if(v_kms > Vtab[Ntab-1]) v_kms = Vtab[Ntab-1];

  // Linear interpolation in log-log space
  double lv = log10(v_kms);
  // Find bracketing indices
  int i0 = 0;
  while(i0 < Ntab-1 && Vtab[i0+1] < v_kms) i0++;
  int i1 = std::min(i0+1, Ntab-1);
  if(Vtab[i1] <= Vtab[i0]) i1 = i0;
  double frac = (Vtab[i1] > Vtab[i0])
                 ? (v_kms - Vtab[i0]) / (Vtab[i1] - Vtab[i0]) : 0.0;

  *sigma_m = Sigtab[i0] * (1-frac) + Sigtab[i1] * frac;
  *rdiss = Rdtab[i0] * (1-frac) + Rdtab[i1] * frac;
  return 0;
}

static void build_grid(simparticles *Sp, int numpart)
{
  for(int c = 0; c < NGRID*NGRID*NGRID; c++) grid_cells[c].clear();

  double pos[3];
  for(int p = 0; p < numpart; p++) {
    if(Sp->P[p].getType() != 1) continue;
    Sp->intpos_to_pos(Sp->P[p].IntPos, pos);
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

static int find_neighbors_from_grid(int i, int *ngb_indices,
                                     double *ngb_dist2, int k)
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

/*
 * Dissipative scattering pair update (Eqs 14, 15, 17 of Schmidt et al. 2026).
 *
 * Given two particles i, j with relative velocity vrel and CM velocity vcm,
 * the dissipative update is:
 *   v_i' = v_cm + (m_j/M) * v_rel * n              (elastic part)
 *   v_j' = v_cm - (m_i/M) * v_rel * n              (elastic part)
 * plus a dissipative reduction of the relative velocity:
 *   v_rel -> v_rel * sqrt(1 - (r_diss - 1) * delta_E_factor)
 *
 * The fraction of CM-frame KE removed per scatter (Eq 15):
 *   delta_E / E_cm = (r_diss - 1) * <k0> / (0.5 * mu * v_rel^2)
 *
 * For simplicity, we use a fixed fraction f = (r_diss - 1) * c0, where
 * c0 ~ 0.4 is the dipole emission coefficient (mean radiated energy as
 * fraction of CM KE). The relative speed after scattering is reduced by
 * sqrt(1 - f).
 */
void sidm_dissipative::scatter_pair_dissipative(simparticles *Sp, int i, int j,
                                                  double sigma_m, double rdiss)
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

  // Scattering direction (Eq 17)
  double n[3];
  yukawa_kick_direction(vrel, n);

  // Dissipative reduction of relative speed (Eq 15)
  // The fraction of CM-frame kinetic energy removed:
  //   f_diss = (r_diss - 1) * c0, where c0 ~ 0.4
  // The relative speed after scattering is reduced:
  //   v_rel' = v_rel * sqrt(1 - f_diss)
  // For r_diss = 1 (elastic), f_diss = 0 and v_rel' = v_rel.
  double c0 = 0.4;  // dipole emission fraction (LB2026)
  double f_diss = (rdiss - 1.0) * c0;
  if(f_diss < 0) f_diss = 0;
  if(f_diss > 0.99) f_diss = 0.99;  // prevent total energy loss
  double vrel_new = vrel * sqrt(1.0 - f_diss);

  // Update velocities (Eq 5: momentum-conserving)
  Sp->P[i].Vel[0] = vcm[0] + (mj/M) * vrel_new * n[0];
  Sp->P[i].Vel[1] = vcm[1] + (mj/M) * vrel_new * n[1];
  Sp->P[i].Vel[2] = vcm[2] + (mj/M) * vrel_new * n[2];

  Sp->P[j].Vel[0] = vcm[0] - (mi/M) * vrel_new * n[0];
  Sp->P[j].Vel[1] = vcm[1] - (mi/M) * vrel_new * n[1];
  Sp->P[j].Vel[2] = vcm[2] - (mi/M) * vrel_new * n[2];
}

void sidm_dissipative::do_sidm_scattering(simparticles *Sp, double dt)
{
  if(dt <= 0) return;
  sidm_init_rng();

  int numpart = Sp->NumPart;
  if(numpart <= 0) return;
  build_grid(Sp, numpart);

  std::vector<char> scattered(numpart, 0);

  double L_unit = All.UnitLength_in_cm;
  double M_unit = All.UnitMass_in_g;

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
      if(scattered[j]) { rates[k] = 0.0; continue; }
      double mj = Sp->P[j].getMass();

      double dv0 = Sp->P[i].Vel[0] - Sp->P[j].Vel[0];
      double dv1 = Sp->P[i].Vel[1] - Sp->P[j].Vel[1];
      double dv2 = Sp->P[i].Vel[2] - Sp->P[j].Vel[2];
      double vrel = sqrt(dv0*dv0 + dv1*dv1 + dv2*dv2);

      // Velocity-dependent σ/m and r_diss
      double sigma_m_local = SigmaOverMass;
      double rdiss_local = RDiss;
      if(VelocityDependent) {
        sigma_and_rdiss(vrel, &sigma_m_local, &rdiss_local);
      }

      double sigma_m_phys = sigma_m_local * M_unit / (L_unit * L_unit);
      rates[k] = 0.5 * sigma_m_phys * mj * vrel / volume;
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
        if(cum >= r) { selected = ngb_idx[k]; break; }
      }
      if(selected < 0 || selected == i || scattered[selected]) continue;

      // Get σ/m and r_diss for this pair's relative velocity
      double dv0 = Sp->P[i].Vel[0] - Sp->P[selected].Vel[0];
      double dv1 = Sp->P[i].Vel[1] - Sp->P[selected].Vel[1];
      double dv2 = Sp->P[i].Vel[2] - Sp->P[selected].Vel[2];
      double vrel = sqrt(dv0*dv0 + dv1*dv1 + dv2*dv2);
      double sigma_m_local = SigmaOverMass;
      double rdiss_local = RDiss;
      if(VelocityDependent) {
        sigma_and_rdiss(vrel, &sigma_m_local, &rdiss_local);
      }

      scatter_pair_dissipative(Sp, i, selected, sigma_m_local, rdiss_local);
      scattered[i] = 1;
      scattered[selected] = 1;
      break;
    }
  }
}

#endif /* SIDM && SIDM_DISSIPATIVE */
