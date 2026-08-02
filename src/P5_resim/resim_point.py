"""
Self-consistent resimulation runner (bypasses rescaling symmetry).

For each candidate physical point (r_s_phys, rho_0_phys, sigma_phys, t_evo_phys)
selected from P3, we run a fresh DissipativeHalo simulation with these
PHYSICAL parameters — no rescaling — and check:

  1. The projected enclosed mass ratio M(r=20pc)/M(r=90pc) at the end of
     the simulation (t = t_evo_phys). If the rescaling symmetry holds
     exactly, this ratio should match the observed 0.364±0.022.
  2. The offset Δ = (M_resim - M_obs)/M_obs is the actual symmetry-breaking
     correction to the mass ratio, NOT a closed-form estimate.

Usage:
    python resim_point.py <point_csv> <out_dir> [--model_key M1]
                                              [--max_seconds 600]
                                              [--save_snapshots]

Outputs:
    <out_dir>/<model_label>_<pick_label>.pickle : per-point results
    <out_dir>/summary.csv : all points + Δ mass ratio
"""
from __future__ import annotations

import os, sys, time, argparse
import numpy as np
import pandas as pd
import pickle

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJ_ROOT = os.path.normpath(os.path.join(_HERE, '..', '..'))
sys.path.insert(0, os.path.join(_PROJ_ROOT, 'src', 'cross_sections'))
sys.path.insert(0, os.path.join(_PROJ_ROOT, 'src', 'fluid_runner'))
sys.path.insert(0, os.path.join(_PROJ_ROOT, 'external', 'gravothermalsidm'))
# P3 rescale module for the projected_enclosed_mass function
sys.path.insert(0, os.path.join(_PROJ_ROOT, 'src', 'P3_rescaling'))

import matplotlib
matplotlib.use('Agg')

from astropy import units as ut
from astropy import constants as ct

from dsidm_models import benchmark_models, sigma_T_born, r_diss, DSIDMParameters
from thermal_avg import effective_sigma_m_and_rdiss
from dissipative_halo import DissipativeHalo
from SourcePy.evolve import Halo
from SourcePy.record import HaloRecord
from rescale import (
    projected_enclosed_mass,
    M_20PC, M_90PC, MASS_RATIO_OBS, MASS_RATIO_ERR,
)


# ----------------------------------------------------------------------
# Simulation parameters (fixed across all resim runs)
# ----------------------------------------------------------------------
N_SHELLS = 100
T_EPSILON = 1e-2
R_EPSILON = 1e-12
RHO_FACTOR_END = 1000.0   # collapse trigger
MAX_STEPS = 200000
SAVE_EVERY = 500
W_UNITS = 100.0  # km/s reference; only used for sigma_m_with_units (calibrated separately)


def evolve_safe(haloevo, t_end_dimless, label, max_steps=MAX_STEPS,
                save_every=SAVE_EVERY, rho_factor_end=RHO_FACTOR_END):
    """Custom evolution loop with NaN check, terminating at t_end_dimless."""
    haloevo.t_epsilon = T_EPSILON
    haloevo.r_epsilon = R_EPSILON

    if haloevo.t == 0:
        if haloevo.flag_hydrostatic_initial:
            haloevo.hydrostatic_adjustment()
        haloevo.save_halo()

    t0 = time.time()
    last_report = time.time()
    n_steps = 0

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

        n_steps = step + 1
        current_rho = haloevo.get_central_quantity(haloevo.rho)
        if haloevo.t >= t_end_dimless:
            haloevo.save_halo()
            print(f"  [{label}] Done: t={haloevo.t:.4f} "
                  f"({(haloevo.t*haloevo.scale_t).to('Gyr').value:.4e} Gyr) "
                  f"after {n_steps} steps")
            break
        elif current_rho > rho_factor_end * haloevo.rho_center:
            haloevo.save_halo()
            t_gyr = (haloevo.t * haloevo.scale_t).to('Gyr').value
            print(f"  [{label}] Collapse trigger at t={haloevo.t:.4f} "
                  f"({t_gyr:.4e} Gyr); stopped before t_end")
            break

        if (step + 1) % save_every == 0:
            haloevo.save_halo()

        if time.time() - last_report > 10.0:
            t_gyr = (haloevo.t * haloevo.scale_t).to('Gyr').value
            rho_now = current_rho * haloevo.scale_rho.to('Msun/pc**3').value
            print(f"  [{label}] step {n_steps}: t={haloevo.t:.2f} "
                  f"({t_gyr:.4e} Gyr), rho_c={rho_now:.2e}")
            last_report = time.time()

    elapsed = time.time() - t0
    print(f"  [{label}] {n_steps} steps in {elapsed:.1f}s")
    return haloevo, n_steps, elapsed


def calibrate_alpha_D_for_sigma(p_model, target_sigma, v_ref=W_UNITS):
    """Adjust alpha_D so sigma_T_born(v_ref) = target_sigma (cm²/g)."""
    sig_now = sigma_T_born(np.array([v_ref]), p_model)[0]
    if sig_now <= 0:
        return p_model
    ratio = target_sigma / sig_now
    # sigma ∝ alpha_D^2 → alpha_D → alpha_D * sqrt(ratio)
    import dataclasses
    return dataclasses.replace(p_model, alpha_D=p_model.alpha_D * np.sqrt(ratio))


def run_one_point(row, out_dir, max_seconds=600, save_snapshots=False):
    """Run a single resimulation for one match point.

    Strategy:
      - Set up an NFW halo with the PHYSICAL r_s_phys and rho_0_phys.
      - Calibrate alpha_D so sigma_T_born(v_ref=100 km/s) = sigma_phys.
        For elastic / const_rdiss_1p05: use Halo with sigma_m_with_units.
      - Evolve to t_evo_phys (converted to dimensionless time).
      - Extract final projected enclosed masses at r=20 pc and r=90 pc.
      - Compute mass ratio and compare to observed.
    """
    model_label = row['model_label']
    pick_label = row['pick_label']
    model_key = row['model']
    label = f"{model_label}_{pick_label}"

    r_s_phys = float(row['r_s_kpc'])         # kpc
    rho_0_phys = float(row['rho0_msun_pc3']) # Msun/pc^3
    sigma_phys = float(row['sigma_m_cm2_g']) # cm^2/g
    t_evo_phys = float(row['t_evo_gyr'])     # Gyr

    print(f"\n[{label}] model={model_key}")
    print(f"  r_s={r_s_phys:.4e} kpc, ρ₀={rho_0_phys:.4e} Msun/pc³")
    print(f"  σ/m={sigma_phys:.4e} cm²/g, t_evo={t_evo_phys:.4e} Gyr")

    # Set up the halo record directory
    out_subdir = os.path.join(out_dir, label)
    if os.path.exists(out_subdir):
        import shutil
        shutil.rmtree(out_subdir)
    os.makedirs(out_subdir, exist_ok=True)
    rec = HaloRecord(out_subdir)

    # Choose halo evolution class based on model
    if model_key in ('M1_dark_photon_massive', 'M2_scalar_phi_massive'):
        # Velocity-dependent: use DissipativeHalo with calibrated model
        bm = benchmark_models()
        p_model = bm[model_key]
        p_cal = calibrate_alpha_D_for_sigma(p_model, sigma_phys)

        sigma_m_eff, rdiss_eff, _, _ = effective_sigma_m_and_rdiss(
            lambda v: sigma_T_born(np.atleast_1d(v), p_cal),
            lambda v: r_diss(np.atleast_1d(v), p_cal),
            np.logspace(1.5, 6.5, 40),
        )
        evo = DissipativeHalo(rec,
                              sigma_m_eff_callable=sigma_m_eff,
                              rdiss_eff_callable=rdiss_eff,
                              flag_dissipation=True,
                              profile='NFW', r_s=r_s_phys, rho_s=rho_0_phys,
                              sigma_m_with_units=sigma_phys, w_units=W_UNITS,
                              n_shells=N_SHELLS, r_max=50.0, r_min=0.02,
                              flag_hydrostatic_initial=True,
                              flag_timestep_use_relaxation=True,
                              flag_timestep_use_energy=True)
    elif model_key in ('elastic', 'const_rdiss_1p05'):
        # Elastic / constant r_diss: use plain Halo with constant sigma_m
        # For const_rdiss_1p05, we use DissipativeHalo with constant r_diss=1.05
        if model_key == 'const_rdiss_1p05':
            sig_fn = lambda T: np.full(np.atleast_1d(T).shape,
                                       sigma_phys, dtype=float)
            rd_fn = lambda T: np.full(np.atleast_1d(T).shape,
                                      1.05, dtype=float)
            evo = DissipativeHalo(rec,
                                  sigma_m_eff_callable=sig_fn,
                                  rdiss_eff_callable=rd_fn,
                                  flag_dissipation=True,
                                  profile='NFW', r_s=r_s_phys, rho_s=rho_0_phys,
                                  sigma_m_with_units=sigma_phys,
                                  w_units=W_UNITS,
                                  n_shells=N_SHELLS, r_max=50.0, r_min=0.02,
                                  flag_hydrostatic_initial=True,
                                  flag_timestep_use_relaxation=True,
                                  flag_timestep_use_energy=True)
        else:
            evo = Halo(rec, profile='NFW', r_s=r_s_phys, rho_s=rho_0_phys,
                       sigma_m_with_units=sigma_phys, w_units=W_UNITS,
                       n_shells=N_SHELLS, r_max=50.0, r_min=0.02,
                       flag_hydrostatic_initial=True,
                       flag_timestep_use_relaxation=True,
                       flag_timestep_use_energy=True)
    else:
        raise ValueError(f"Unknown model_key: {model_key}")

    # Compute t_end in dimensionless units
    # t_dimless = t_phys / scale_t   where scale_t = 1/sqrt(4*pi*G*rho_s)
    rho_s_astropy = rho_0_phys * ut.M_sun / ut.pc**3
    scale_t_gyr = (1.0 / np.sqrt(4.0 * np.pi * rho_s_astropy * ct.G)).to('Gyr').value
    t_end_dimless = t_evo_phys / scale_t_gyr
    print(f"  scale_t={scale_t_gyr:.4e} Gyr → t_end_dimless={t_end_dimless:.4e}")

    # Set up time budget
    t_start = time.time()

    # Run the evolution
    evo, n_steps, elapsed = evolve_safe(evo, t_end_dimless, label=label)

    # Extract final snapshot
    list_files, list_times = rec.glob_pickle_files()
    if len(list_files) == 0:
        print(f"  [{label}] No snapshots saved")
        return None

    last_file = list_files[-1]
    data = rec.get_halo_state_pickled(file_halo=last_file)
    scale_r_kpc = evo.scale_r.to('kpc').value
    scale_rho = evo.scale_rho.to('Msun/pc**3').value

    r_kpc = data['r'] * scale_r_kpc
    rho_phys_arr = data['rho'] * scale_rho  # Msun/pc^3

    # Compute projected enclosed masses at r=20 pc and r=90 pc
    # Note: r=20 pc = 0.02 kpc, r=90 pc = 0.09 kpc
    M_inner = projected_enclosed_mass(r_kpc, rho_phys_arr, 0.02,
                                       r_unit='kpc', rho_unit='Msun_pc3')
    M_outer = projected_enclosed_mass(r_kpc, rho_phys_arr, 0.09,
                                       r_unit='kpc', rho_unit='Msun_pc3')

    if M_outer > 0:
        mass_ratio_resim = M_inner / M_outer
    else:
        mass_ratio_resim = np.nan

    delta_ratio = mass_ratio_resim - MASS_RATIO_OBS
    delta_sigma_pct = 100.0 * delta_ratio / MASS_RATIO_OBS

    # Symmetry-breaking correction = the residual offset of the resim mass
    # ratio from the observed value. Under perfect elastic symmetry, this
    # should be 0. Non-zero values indicate the symmetry-breaking correction
    # that must be applied to the naive σ/m to match the observation.
    # Specifically, if Δ > 0 (resim ratio > observed), the halo has over-cooled
    # → σ/m should be reduced; the σ/m_correction factor is approximately
    # sqrt(M_obs / M_resim) under the LMFP limit (where t_cool ∝ 1/σ).

    result = {
        'model_label': model_label,
        'pick_label': pick_label,
        'model_key': model_key,
        'r_s_kpc': r_s_phys,
        'rho0_msun_pc3': rho_0_phys,
        'sigma_m_cm2_g': sigma_phys,
        't_evo_gyr': t_evo_phys,
        'M_inner_resim': float(M_inner),
        'M_outer_resim': float(M_outer),
        'mass_ratio_resim': float(mass_ratio_resim),
        'mass_ratio_obs': MASS_RATIO_OBS,
        'mass_ratio_obs_err': MASS_RATIO_ERR,
        'delta_ratio': float(delta_ratio),
        'delta_ratio_pct': float(delta_sigma_pct),
        'n_steps': n_steps,
        'elapsed_s': elapsed,
        't_final_gyr': float(data['t'] * scale_t_gyr),
        't_final_dimless': float(data['t']),
    }

    # Save result
    out_pickle = os.path.join(out_dir, f"{label}.pickle")
    with open(out_pickle, 'wb') as f:
        pickle.dump(result, f)
    print(f"  [{label}] M(20pc)/M(90pc)_resim = {mass_ratio_resim:.4f} "
          f"(obs {MASS_RATIO_OBS:.4f}±{MASS_RATIO_ERR:.4f}) "
          f"Δ={delta_ratio:+.4f} ({delta_sigma_pct:+.2f}%)")
    print(f"  [{label}] Saved → {out_pickle}")

    # Clean up snapshots unless explicitly saved
    if not save_snapshots:
        import shutil
        shutil.rmtree(out_subdir, ignore_errors=True)

    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('point_csv', help='CSV of resim points (from select_points.py)')
    parser.add_argument('out_dir', help='Output directory for results')
    parser.add_argument('--model_label', default=None,
                        help='Filter by model_label (M1, M2, M3). Default: all.')
    parser.add_argument('--pick_label', default=None,
                        help='Filter by pick_label (P1_lower, P2_p25, etc.). Default: all.')
    parser.add_argument('--max_seconds', type=int, default=600,
                        help='Per-point wall-clock budget (default 600s)')
    parser.add_argument('--save_snapshots', action='store_true',
                        help='Keep snapshot directories (default: delete)')
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    df = pd.read_csv(args.point_csv)
    print(f"Loaded {len(df)} resim points from {args.point_csv}")

    if args.model_label:
        df = df[df['model_label'] == args.model_label]
    if args.pick_label:
        df = df[df['pick_label'] == args.pick_label]
    print(f"Will run {len(df)} simulations")

    results = []
    summary_csv = os.path.join(args.out_dir, 'summary.csv')
    for i, row in df.iterrows():
        print(f"\n=== Run {i+1}/{len(df)} ===")
        try:
            r = run_one_point(row, args.out_dir,
                              max_seconds=args.max_seconds,
                              save_snapshots=args.save_snapshots)
            if r is not None:
                results.append(r)
        except Exception as e:
            import traceback
            print(f"  ERROR: {e}")
            traceback.print_exc()
            results.append({
                'model_label': row['model_label'],
                'pick_label': row['pick_label'],
                'model_key': row['model'],
                'error': str(e),
            })

        # Periodically write summary CSV so progress is tracked
        if results:
            pd.DataFrame(results).to_csv(summary_csv, index=False)
            print(f"  → updated {summary_csv}")

    print(f"\n=== Done: {len(results)} runs ===")
    if results:
        df_out = pd.DataFrame(results)
        df_out.to_csv(summary_csv, index=False)
        print(f"Final summary → {summary_csv}")
        print("\nSummary table:")
        cols = ['model_label', 'pick_label', 'sigma_m_cm2_g', 't_evo_gyr',
                'mass_ratio_resim', 'delta_ratio_pct']
        cols = [c for c in cols if c in df_out.columns]
        print(df_out[cols].to_string(index=False))


if __name__ == '__main__':
    main()
