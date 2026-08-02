#ifndef SIDM_H
#define SIDM_H

#include "gadgetconfig.h"

#ifdef SIDM

#include "../data/simparticles.h"

class sidm
{
 public:
  sidm();
  ~sidm();

  void do_sidm_scattering(simparticles *Sp, double dt);
  double get_sidm_timestep(int i) { return 0; }

 private:
  int find_neighbors(simparticles *Sp, int i, double *pos_i, int *ngb_indices, double *ngb_dist2, int k);
  void scatter_pair(simparticles *Sp, int i, int j);
  int poisson_sample(double lambda);
  void random_unit_vector(double *v);

  double SigmaOverMass;
  int KNeighbors;
  double PMax;
};

extern sidm Sidm;

#endif /* SIDM */

#endif /* SIDM_H */
