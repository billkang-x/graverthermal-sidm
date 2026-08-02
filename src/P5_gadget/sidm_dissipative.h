#ifndef SIDM_DISSIPATIVE_H
#define SIDM_DISSIPATIVE_H

/*
 * Dissipative fSIDM extension for GADGET-4
 *
 * Implements the frequent small-angle self-interaction framework with
 * dissipation from Schmidt, Fischer & Garny 2026 (arXiv:2606.19428).
 *
 * Key additions over the elastic isotropic SIDM (sidm.cc):
 *   - Drag force (Eq 14): F_drag = -mu_ij * <sigma v> * (v_i - v_j)
 *   - Modified random kick (Eq 17): kicks are forward-peaked (Yukawa-like)
 *       instead of isotropic, in the frequent small-angle regime.
 *   - Momentum-conserving update (Eq 5): ensures momentum conservation
 *       between paired particles.
 *   - Energy loss / dissipation (Eq 15): a fraction (r_diss - 1) of the
 *       kinetic energy in the center-of-mass frame is removed per scatter.
 *
 * Velocity-dependent cross sections are supported via a tabulated
 * sigma_T(v) and r_diss(v) loaded at startup.
 */

#include "gadgetconfig.h"

#if defined(SIDM) && defined(SIDM_DISSIPATIVE)

#include "../data/simparticles.h"

class sidm_dissipative
{
 public:
  sidm_dissipative();
  ~sidm_dissipative();

  void do_sidm_scattering(simparticles *Sp, double dt);
  double get_sidm_timestep(int i) { return 0; }

  // Velocity-dependent σ_T/m and r_diss at relative speed vrel (code units).
  // Returns 0 on success.
  int sigma_and_rdiss(double vrel_code, double *sigma_m, double *rdiss);

 private:
  int find_neighbors(simparticles *Sp, int i, double *pos_i,
                    int *ngb_indices, double *ngb_dist2, int k);
  void scatter_pair_dissipative(simparticles *Sp, int i, int j,
                                double sigma_m, double rdiss);
  int poisson_sample(double lambda);
  void random_unit_vector(double *v);
  // Forward-peaked (Yukawa) random scattering angle
  void yukawa_kick_direction(double vrel, double *n);

  double SigmaOverMass;     // fiducial σ_T/m at v_ref (code units)
  double RDiss;             // dissipation parameter r_diss (>=1)
  int KNeighbors;
  double PMax;
  int VelocityDependent;   // 0 = constant, 1 = use sigma(v) table
  double VRef;              // reference velocity for sigma normalization (km/s)

  // Velocity-dependent cross section table (loaded from file if VelocityDependent=1)
  double *Vtab;             // v axis (km/s)
  double *Sigtab;           // σ/m at v (cm^2/g)
  double *Rdtab;            // r_diss at v (dimensionless)
  int Ntab;
};

extern sidm_dissipative SidmDissipative;

#endif /* SIDM && SIDM_DISSIPATIVE */
#endif /* SIDM_DISSIPATIVE_H */
