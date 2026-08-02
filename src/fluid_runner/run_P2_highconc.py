"""
P2 high-concentration NFW run to probe the B1938+666 extended-solution branch.

Motivation:
    Our default P2 runs use r_s=3.6 kpc, rho_0=7.09e-3 Msun/pc^3, which gives
    c_200 ~ 8.5 for r_vir ~ 30 kpc. This probes the COMPACT-solution branch
    of B1938+666 (r_2D/r_s ~ 0.02-0.7, sigma_phys ~ 0.01-0.1 cm^2/g). The
    EXTENDED-solution branch (Row #2 of Table 2 in Schmidt 2026, sigma_phys
    ~ 5 cm^2/g, t_evo ~ 6 Gyr) requires probing larger r_2D/r_s, which
    corresponds to a higher-concentration initial NFW profile.

    To address review #9 (B1938 sampling bias), we run a M1 dark photon
    simulation with c_200 ~ 12 (r_s ~ 2.5 kpc, rho_0 ~ 1.5e-2 Msun/pc^3,
    scaling the density up by ~2x). The simulation should produce matches
    at larger r_2D/r_s, where v_phys can be higher, allowing the
    symmetry-breaking correction to be measured.

Outputs:
    data/P2_runs/M1_dark_photon_massive_highconc/ : snapshots
    (Then re-run P3 rescaling on this directory.)
"""
from __future__ import annotations

import os, sys, time, shutil
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJ_ROOT = os.path.normpath(os.path.join(_HERE, '..', '..'))
sys.path.insert(0, os.path.join(_PROJ_ROOT, 'src', 'cross_sections'))
sys.path.insert(0, os.path.join(_PROJ_ROOT, 'src', 'fluid_runner'))
sys.path.insert(0, os.path.join(_PROJ_ROOT, 'external', 'gravothermalsidm'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from dsidm_models import benchmark_models, sigma_T_born, r_diss, DSIDMParameters
from thermal_avg import effective_sigma_m_and_rdiss
from dissipative_halo import DissipativeHalo
from SourcePy.evolve import Halo
from SourcePy.record import HaloRecord


# ----------------------------------------------------------------------
# High-concentration parameters (c_200 ~ 12 vs default c_200 ~ 8.5)
# ----------------------------------------------------------------------
# Strategy: keep r_vir fixed at the default ~30 kpc, but reduce r_s to
# increase concentration. To keep the same velocity dispersion scale
# (v ~ sqrt(G rho_s r_s^2) ~ constant) we scale rho_s inversely with
# r_s^2.
#
# Default: r_s = 3.6 kpc, rho_0 = 7.09e-3 Msun/pc^3 (c ~ 8.5)
# Highconc: r_s = 2.5 kpc (c ~ 12), rho_0 ~ 7.09e-3 * (3.6/2.5)^2 ~ 1.47e-2
# This keeps the velocity scale and time scale fixed while shifting the
# matching radii.
RS_KPC = 2.5
RHO0 = 1.47e-2  # Msun/pc^3
SIGMA_M = 50.0  # cm^2/g (same as default)
W_UNITS = 100.0  # km/s
N_SHELLS = 100
T_EPSILON = 1e-2
R_EPSILON = 1e-12
T_END = 150.0
RHO_FACTOR_END = 1000.0
MAX_STEPS = 500000
SAVE_EVERY = 500


def evolve_safe(haloevo, t_end=T_END, rho_factor_end=RHO_FACTOR_END,
                max_steps=MAX_STEPS, save_every=SAVE_EVERY,
                t_epsilon=T_EPSILON, r_epsilon=R_EPSILON, label=""):
    haloevo.t_epsilon = t_epsilon
    haloevo.r_epsilon = r_epsilon

    if haloevo.t == 0:
        if haloevo.flag_hydrostatic_initial:
            haloevo.hydrostatic_adjustment()
        haloevo.save_halo()

    t0 = time.time()
    last_report = time.time()

    for step in range(max_steps):
        try:
            haloevo.conduct_heat()
            haloevo.hydrostatic_adjustment()
        except Exception as e:
            print(f"  [{label}] Error at step {step}: {e}")
            break

        if np.any(np.isnan(haloevo.rho)) or np.any(np.isnan(haloevo.r)):
            print(f"  [{label}] NaN at step {step}, t={haloevo.t:.4f}")
            break

        current_rho = haloevo.get_central_quantity(haloevo.rho)
        if haloevo.t >= t_end:
            haloevo.save_halo()
            t_gyr = haloevo.t * haloevo.scale_t.to('Gyr').value
            print(f"  [{label}] Done: t={haloevo.t:.2f} ({t_gyr:.2f} Gyr)")
            break
        elif current_rho > rho_factor_end * haloevo.rho_center:
            haloevo.save_halo()
            print(f"  [{label}] Collapse: rho>{rho_factor_end}*rho0 at t={haloevo.t:.2f}")
            break

        if (step + 1) % save_every == 0:
            haloevo.save_halo()

        if time.time() - last_report > 10.0:
            t_gyr = haloevo.t * haloevo.scale_t.to('Gyr').value
            rho_now = current_rho * haloevo.scale_rho.to('Msun/pc**3').value
            print(f"  [{label}] step {step+1}: t={haloevo.t:.1f} ({t_gyr:.2f} Gyr), rho_c={rho_now:.2e}")
            last_report = time.time()

    elapsed = time.time() - t0
    print(f"  [{label}] {haloevo.n_conduction} steps in {elapsed:.1f}s")
    return haloevo


def calibrate_model(p, target=50.0, v_ref=100.0):
    sig = sigma_T_born(np.array([v_ref]), p)[0]
    ratio = target / sig
    import dataclasses
    return dataclasses.replace(p, alpha_D=p.alpha_D * np.sqrt(ratio))


def main():
    dir_output = os.path.join(_PROJ_ROOT, 'data', 'P2_runs')
    os.makedirs(dir_output, exist_ok=True)

    print("=" * 80)
    print("P2 high-concentration NFW run (c_200 ~ 12)")
    print(f"NFW: r_s={RS_KPC}, rho_0={RHO0}, sigma_m={SIGMA_M}, w={W_UNITS}")
    print("=" * 80)

    # Run only M1 (the most interesting model for symmetry breaking)
    bm = benchmark_models()
    p = bm['M1_dark_photon_massive']
    p = calibrate_model(p, target=SIGMA_M, v_ref=W_UNITS)
    print(f"\nCalibrated M1: alpha_D = {p.alpha_D:.4e}")

    name = 'M1_dark_photon_massive_highconc'
    dir_data = os.path.join(dir_output, name)
    pickles = [f for f in os.listdir(dir_data) if f.endswith('.pickle')] \
        if os.path.exists(dir_data) else []
    if len(pickles) > 10:
        print(f"  [{name}] Already done ({len(pickles)} snapshots), skipping.")
        return

    if os.path.exists(dir_data):
        shutil.rmtree(dir_data)
    rec = HaloRecord(dir_data)

    sigma_m_func, rdiss_func, _, _ = effective_sigma_m_and_rdiss(
        lambda v: sigma_T_born(np.atleast_1d(v), p),
        lambda v: r_diss(np.atleast_1d(v), p),
        np.logspace(1.5, 6.5, 40)
    )
    evo = DissipativeHalo(rec,
                          sigma_m_eff_callable=sigma_m_func,
                          rdiss_eff_callable=rdiss_func,
                          flag_dissipation=True,
                          profile='NFW', r_s=RS_KPC, rho_s=RHO0,
                          sigma_m_with_units=SIGMA_M, w_units=W_UNITS,
                          n_shells=N_SHELLS, r_max=50.0, r_min=0.02,
                          flag_hydrostatic_initial=True,
                          flag_timestep_use_relaxation=True,
                          flag_timestep_use_energy=True)
    evolve_safe(evo, label=name)
    print(f"\n[Done] {name}")


if __name__ == '__main__':
    main()
