# P1 Elastic Relaxation Rerun

## Purpose

The original P1 (elastic control) run showed a -18.8% drop in the projected mass
ratio M_2D(<R_inner)/M_2D(<R_outer), far exceeding the fluid-model prediction of
-2.6%. The drop is dominated by initial-condition disequilibrium: the IC
generator samples positions from NFW and assigns isotropic Gaussian velocities
with sigma = v_circ/sqrt(3), which is NOT the Jeans-equilibrium solution for
NFW. This produces a virial-ratio offset and a transient mass redistribution in
the first ~0.5 code-time units.

This rerun isolates that effect by:

  Phase A — Relaxation: run the SAME initial conditions with PURE ELASTIC
            SIDM (sigma/m = 0.1, r_diss = 1.0, exactly as the original P1
            Config) for a short time t_relax = 0.5 code units (≈ 0.49 Gyr,
            ~5 dynamical times at r_s). No snapshots are written except the
            final state, which becomes the new equilibrium IC ("restart
            snapshot").
  Phase B — Production: restart from the Phase-A snapshot as t=0 and run the
            production P1 (elastic) for t_prod = 1.0 code units, with
            snapshots at t=0, 0.5, 1.0 — mirroring the original run's output
            schedule.

## Expected outcome

If the -18.8% drop was IC disequilibrium, Phase B should show a drop of at
most a few percent (consistent with the fluid model's -2.6% within the ~4%
projection error). If Phase B still drops by ~15% or more, the cause is NOT
IC disequilibrium but a code/physics modelling issue, and the current "P3
pass" conclusion must be revisited.

## Files in this directory

  relax_params.txt            Phase A parameter file (relaxation run)
  relax_Config.sh              Phase A Config.sh (elastic, no dissipation)
  relax_submit.sh             Phase A SLURM submission script
  prod_params.txt              Phase B parameter file (production, restart)
  prod_Config.sh              Phase B Config.sh (identical to original P1)
  prod_submit.sh              Phase B SLURM submission script
  convert_restart_to_ic.py    Convert Phase-A snapshot into a restart IC
  analyze_relaxation.py       Compare new P1 ratios with original + fluid
  relaxation_workflow.md      Step-by-step execution guide

## Key parameters

  lambda = 0.085 / 3.6 = 0.023611          (rescaling length)
  mu     = (10/7.09e-3) * lambda^3          (rescaling mass)
  T_SCALE = sqrt(lambda^3/mu) = 0.026627     (t_phys = T_SCALE * t_sim)
  t_relax_code = 0.5   ->  t_relax_phys ≈ 0.0133 Gyr (≈ 5 t_dyn at r_s)
  t_prod_code  = 1.0   ->  t_prod_phys  ≈ 0.0266 Gyr  (matches original)

  IC: same as original run — D:/graverthermal-sidm/data/P5_nbody_verify/ics_sim_space/ic.dat
      (r_s=3.6 kpc, rho_0=7.09e-3 Msun/pc^3, N=500000, PartType1 only)

  Projection radii (simulation space, invariant under rescaling):
    R_inner_sim = 0.847 kpc  (from 20 pc)
    R_outer_sim = 3.812 kpc  (from 90 pc)

## Original P1 result (for reference)

  t_code  N-body ratio   Fluid ratio   N-body rel   Fluid rel
  0.0     0.084074       0.146506      1.000        1.000
  0.5     0.068563       0.140624      0.8155       0.9599
  1.0     0.068248       0.142738      0.8118       0.9743

  Drop at t=1.0: -18.8% (N-body) vs -2.6% (fluid).
