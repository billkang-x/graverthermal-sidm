"""
Self-consistent B1938+666 viable region scan via 2D grid simulation.

Scans (sigma/m, r_s) directly, with fixed rho_0 chosen to give
dynamical timescales matching B1938. The initial NFW mass ratio
M(20pc)/M(90pc) depends only on r_s (for NFW, ratio is independent
of rho_0), so r_s < ~0.1 kpc gives init_ratio > 0.364, and collapse
drives it down through the observed value.

Grid:
  - sigma/m: 20 points, log-spaced 1e-3 to 1.0 cm²/g
  - r_s: 15 points, log-spaced 0.01 to 0.15 kpc
  - rho_0 fixed at 10 Msun/pc^3 (P3 median)
  - 3 models × 300 = 900 simulations

For each point: evolve, record mass ratio at checkpoints, find
t_cross where ratio = 0.364. Viable if t_cross <= 6.37 Gyr.
"""
from __future__ import annotations

import os, sys, time, argparse, traceback
import numpy as np
import pandas as pd
import pickle

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJ_ROOT = os.path.normpath(os.path.join(_HERE, '..', '..'))
sys.path.insert(0, os.path.join(_PROJ_ROOT, 'src', 'cross_sections'))
sys.path.insert(0, os.path.join(_PROJ_ROOT, 'src', 'fluid_runner'))
sys.path.insert(0, os.path.join(_PROJ_ROOT, 'external', 'gravothermalsidm'))
sys.path.insert(0, os.path.join(_PROJ_ROOT, 'src', 'P3_rescaling'))

import matplotlib
matplotlib.use('Agg')

from astropy import units as ut
from astropy import constants as ct

from dsidm_models import benchmark_models, sigma_T_born, r_diss
from thermal_avg import effective_sigma_m_and_rdiss
from dissipative_halo import DissipativeHalo
from SourcePy.evolve import Halo
from SourcePy.record import HaloRecord
from rescale import (
    projected_enclosed_mass,
    M_20PC, M_90PC, MASS_RATIO_OBS, MASS_RATIO_ERR,
)

# ----------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------
T_OBS_GYR = 6.37
RHO_0_FIXED = 10.0  # Msun/pc^3 (P3 median; gives t_dyn ~ 0.1 Gyr)

N_SHELLS = 100
T_EPSILON = 1e-2
R_EPSILON = 1e-12
RHO_FACTOR_END = 1000.0
MAX_STEPS = 30000  # Reduced: 30k steps is enough to reach t~40
SAVE_EVERY = 1000
W_UNITS = 100.0

# Time checkpoints (dimensionless), denser near t=0
# Extended for M3: weak dissipation means slower evolution.
# t=1280 dimless ~ 1.7 Gyr, t=5120 ~ 6.8 Gyr (covers t_obs=6.37)
T_CHECKPOINTS = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0,
                 80.0, 160.0, 320.0, 640.0, 1280.0, 2560.0, 5120.0]


def compute_initial_ratio(r_s_kpc, rho_0):
    """Initial NFW projected mass ratio M(20pc)/M(90pc).
    For NFW, ratio is independent of rho_0 (cancels)."""
    r_arr = np.logspace(-2, 1, 200)
    x = r_arr / r_s_kpc
    rho_arr = rho_0 / (x * (1.0 + x)**2)
    M_in = projected_enclosed_mass(r_arr, rho_arr, 0.02,
                                    r_unit='kpc', rho_unit='Msun_pc3')
    M_out = projected_enclosed_mass(r_arr, rho_arr, 0.09,
                                     r_unit='kpc', rho_unit='Msun_pc3')
    return M_in / M_out if M_out > 0 else np.nan


SIGMA_MIN, SIGMA_MAX = 0.005, 1.0
RS_MIN, RS_MAX = 0.01, 0.12

MODELS = {
    'M1': {'key': 'M1_dark_photon_massive', 'label': 'M1: dark photon'},
    'M2': {'key': 'M2_scalar_phi_massive',  'label': 'M2: scalar $\\phi$'},
    'M3': {'key': 'const_rdiss_1p05',       'label': 'M3: const $r_{\\rm diss}=1.05$'},
}


def calibrate_alpha_D(p_model, target_sigma, v_ref=W_UNITS):
    sig_now = sigma_T_born(np.array([v_ref]), p_model)[0]
    if sig_now <= 0:
        return p_model
    ratio = target_sigma / sig_now
    import dataclasses
    return dataclasses.replace(p_model, alpha_D=p_model.alpha_D * np.sqrt(ratio))


def evolve_with_checkpoints(haloevo, checkpoints, label, max_steps=MAX_STEPS,
                             stop_after_cross=True, wall_budget_s=60):
    """Evolve halo, recording mass ratio at each checkpoint.

    stop_after_cross: stop once crossing found.
    wall_budget_s: hard wall-clock limit (default 60s).
    """
    haloevo.t_epsilon = T_EPSILON
    haloevo.r_epsilon = R_EPSILON
    rec = haloevo.record

    if haloevo.t == 0:
        if haloevo.flag_hydrostatic_initial:
            haloevo.hydrostatic_adjustment()
        haloevo.save_halo()

    checkpoints_sorted = sorted(checkpoints)
    results = []
    next_ckpt_idx = 0
    found_cross = False
    t_start_wall = time.time()

    for step in range(max_steps):
        # Wall-clock budget check
        if time.time() - t_start_wall > wall_budget_s:
            break

        try:
            haloevo.conduct_heat()
            haloevo.hydrostatic_adjustment()
        except Exception:
            break

        if np.any(np.isnan(haloevo.rho)) or np.any(np.isnan(haloevo.r)):
            break

        while (next_ckpt_idx < len(checkpoints_sorted)
               and haloevo.t >= checkpoints_sorted[next_ckpt_idx]):
            haloevo.save_halo()
            list_files, _ = rec.glob_pickle_files()
            if len(list_files) > 0:
                data = rec.get_halo_state_pickled(file_halo=list_files[-1])
                scale_r_kpc = haloevo.scale_r.to('kpc').value
                scale_rho = haloevo.scale_rho.to('Msun/pc**3').value
                r_kpc = data['r'] * scale_r_kpc
                rho_arr = data['rho'] * scale_rho
                M_in = projected_enclosed_mass(r_kpc, rho_arr, 0.02,
                                               r_unit='kpc', rho_unit='Msun_pc3')
                M_out = projected_enclosed_mass(r_kpc, rho_arr, 0.09,
                                                r_unit='kpc', rho_unit='Msun_pc3')
                ratio = M_in / M_out if M_out > 0 else np.nan
                results.append({
                    't_dimless': float(data['t']),
                    't_gyr': float(data['t'] * haloevo.scale_t.to('Gyr').value),
                    'mass_ratio': float(ratio),
                })
                if stop_after_cross and len(results) >= 2:
                    d_prev = results[-2]['mass_ratio'] - MASS_RATIO_OBS
                    d_curr = results[-1]['mass_ratio'] - MASS_RATIO_OBS
                    if d_prev * d_curr < 0:
                        found_cross = True
                # Early bail: only bail if ratio is INCREASING (thermal expansion
                # dominating) AND well above obs. For M3 (weak dissipation),
                # ratio decreases slowly so we must not bail early.
                if len(results) >= 5 and ratio > MASS_RATIO_OBS + 0.15:
                    # Check trend: is ratio still decreasing?
                    r_recent = [results[i]['mass_ratio'] for i in range(-3, 0)]
                    if r_recent[-1] >= r_recent[0]:  # not decreasing
                        found_cross = True  # stop signal, mark as no-cross
            next_ckpt_idx += 1

        if found_cross:
            break
        if next_ckpt_idx >= len(checkpoints_sorted):
            break

        current_rho = haloevo.get_central_quantity(haloevo.rho)
        if current_rho > RHO_FACTOR_END * haloevo.rho_center:
            haloevo.save_halo()
            break

        if (step + 1) % SAVE_EVERY == 0:
            haloevo.save_halo()

    return results


def run_one_grid_point(model_label, model_key, sigma_m, r_s_kpc, out_dir,
                       rho_0=RHO_0_FIXED):
    label = f"{model_label}_sig{sigma_m:.4e}_rs{r_s_kpc:.4f}"

    init_ratio = compute_initial_ratio(r_s_kpc, rho_0)
    if init_ratio < MASS_RATIO_OBS - 3 * MASS_RATIO_ERR:
        return {
            'model_label': model_label, 'model_key': model_key,
            'sigma_m_cm2_g': sigma_m, 'r_s_kpc': r_s_kpc,
            'rho_0_msun_pc3': rho_0,
            't_cross_gyr': 0.0, 'viable': False,
            'n_checkpoints': 0, 'final_ratio': init_ratio, 'final_t_gyr': 0.0,
            'init_ratio': init_ratio, 'skipped': 'init_below_obs',
        }

    out_subdir = os.path.join(out_dir, '_tmp', label)
    os.makedirs(out_subdir, exist_ok=True)
    rec = HaloRecord(out_subdir)

    try:
        if model_key in ('M1_dark_photon_massive', 'M2_scalar_phi_massive'):
            bm = benchmark_models()
            p_model = bm[model_key]
            p_cal = calibrate_alpha_D(p_model, sigma_m)
            sigma_m_eff, rdiss_eff, _, _ = effective_sigma_m_and_rdiss(
                lambda v: sigma_T_born(np.atleast_1d(v), p_cal),
                lambda v: r_diss(np.atleast_1d(v), p_cal),
                np.logspace(1.5, 6.5, 40),
            )
            evo = DissipativeHalo(rec,
                                  sigma_m_eff_callable=sigma_m_eff,
                                  rdiss_eff_callable=rdiss_eff,
                                  flag_dissipation=True,
                                  profile='NFW', r_s=r_s_kpc, rho_s=rho_0,
                                  sigma_m_with_units=sigma_m, w_units=W_UNITS,
                                  n_shells=N_SHELLS, r_max=50.0, r_min=0.02,
                                  flag_hydrostatic_initial=True,
                                  flag_timestep_use_relaxation=True,
                                  flag_timestep_use_energy=True)
        elif model_key == 'const_rdiss_1p05':
            sig_fn = lambda T: np.full(np.atleast_1d(T).shape, sigma_m, dtype=float)
            rd_fn = lambda T: np.full(np.atleast_1d(T).shape, 1.05, dtype=float)
            evo = DissipativeHalo(rec,
                                  sigma_m_eff_callable=sig_fn,
                                  rdiss_eff_callable=rd_fn,
                                  flag_dissipation=True,
                                  profile='NFW', r_s=r_s_kpc, rho_s=rho_0,
                                  sigma_m_with_units=sigma_m, w_units=W_UNITS,
                                  n_shells=N_SHELLS, r_max=50.0, r_min=0.02,
                                  flag_hydrostatic_initial=True,
                                  flag_timestep_use_relaxation=True,
                                  flag_timestep_use_energy=True)
        else:
            raise ValueError(f"Unknown model_key: {model_key}")

        results = evolve_with_checkpoints(evo, T_CHECKPOINTS, label=label,
                                            stop_after_cross=True, wall_budget_s=30)

        t_cross_gyr = np.nan
        if len(results) >= 2:
            ts = np.array([r['t_gyr'] for r in results])
            rs = np.array([r['mass_ratio'] for r in results])
            diff = rs - MASS_RATIO_OBS
            for i in range(len(diff) - 1):
                if diff[i] * diff[i + 1] < 0:
                    lt = np.log10(ts[i]) + (np.log10(ts[i+1]) - np.log10(ts[i])) * (-diff[i]) / (diff[i+1] - diff[i])
                    t_cross_gyr = 10.0 ** lt
                    break

        if np.isnan(t_cross_gyr):
            if len(results) > 0:
                final_ratio = results[-1]['mass_ratio']
                if final_ratio > MASS_RATIO_OBS:
                    t_cross_gyr = np.inf
                else:
                    t_cross_gyr = 0.0

        result = {
            'model_label': model_label, 'model_key': model_key,
            'sigma_m_cm2_g': sigma_m, 'r_s_kpc': r_s_kpc,
            'rho_0_msun_pc3': rho_0,
            't_cross_gyr': t_cross_gyr,
            'viable': bool(not np.isnan(t_cross_gyr) and t_cross_gyr <= T_OBS_GYR),
            'n_checkpoints': len(results),
            'final_ratio': results[-1]['mass_ratio'] if results else np.nan,
            'final_t_gyr': results[-1]['t_gyr'] if results else np.nan,
            'init_ratio': init_ratio,
            'skipped': 'no',
        }
        if results:
            result['checkpoints'] = results

    except Exception as e:
        result = {
            'model_label': model_label, 'model_key': model_key,
            'sigma_m_cm2_g': sigma_m, 'r_s_kpc': r_s_kpc,
            'rho_0_msun_pc3': rho_0,
            't_cross_gyr': np.nan, 'viable': False,
            'n_checkpoints': 0, 'final_ratio': np.nan, 'final_t_gyr': np.nan,
            'init_ratio': init_ratio, 'skipped': f'error: {e}',
        }
    finally:
        import shutil
        shutil.rmtree(out_subdir, ignore_errors=True)

    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('out_dir', help='Output directory')
    parser.add_argument('--model', default='all', choices=['M1', 'M2', 'M3', 'all'])
    parser.add_argument('--n_sigma', type=int, default=20)
    parser.add_argument('--n_rs', type=int, default=15)
    parser.add_argument('--sigma_min', type=float, default=SIGMA_MIN)
    parser.add_argument('--sigma_max', type=float, default=SIGMA_MAX)
    parser.add_argument('--rs_min', type=float, default=RS_MIN)
    parser.add_argument('--rs_max', type=float, default=RS_MAX)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(os.path.join(args.out_dir, '_tmp'), exist_ok=True)

    sigma_grid = np.logspace(np.log10(args.sigma_min), np.log10(args.sigma_max), args.n_sigma)
    rs_grid = np.logspace(np.log10(args.rs_min), np.log10(args.rs_max), args.n_rs)
    print(f"Grid: {args.n_sigma} sigma x {args.n_rs} r_s = {args.n_sigma * args.n_rs} points per model")
    print(f"  sigma/m: [{sigma_grid[0]:.4e}, {sigma_grid[-1]:.4e}] cm²/g")
    print(f"  r_s:     [{rs_grid[0]:.4f}, {rs_grid[-1]:.4f}] kpc")
    print(f"  rho_0:   {RHO_0_FIXED} Msun/pc^3 (fixed)")

    # Print init_ratio at grid corners
    for rs in [rs_grid[0], rs_grid[-1]]:
        ir = compute_initial_ratio(rs, RHO_0_FIXED)
        print(f"  init_ratio at r_s={rs:.4f}: {ir:.4f}")

    models_to_run = list(MODELS.keys()) if args.model == 'all' else [args.model]
    total_runs = len(models_to_run) * args.n_sigma * args.n_rs
    print(f"Models: {models_to_run} -> {total_runs} total simulations")

    summary_csv = os.path.join(args.out_dir, 'grid_summary.csv')
    all_results = []

    t_start = time.time()
    run_count = 0

    for ml in models_to_run:
        mk = MODELS[ml]['key']
        print(f"\n=== Model {ml} ({mk}) ===")

        for i_sig, sig in enumerate(sigma_grid):
            for i_rs, rs in enumerate(rs_grid):
                run_count += 1
                elapsed = time.time() - t_start
                rate = run_count / max(elapsed, 1)
                eta = (total_runs - run_count) / max(rate, 0.01)
                print(f"[{run_count}/{total_runs}] {ml} sig={sig:.4e} rs={rs:.4f} "
                      f"({elapsed:.0f}s, {rate:.1f}/s, ETA {eta:.0f}s)", end=' ', flush=True)

                result = run_one_grid_point(ml, mk, sig, rs, args.out_dir)
                all_results.append(result)

                if result.get('skipped', 'no') != 'no':
                    print(f"SKIP ({result['skipped']})")
                elif result.get('viable'):
                    print(f"VIABLE t_cross={result['t_cross_gyr']:.4f} Gyr")
                elif result['t_cross_gyr'] == np.inf:
                    print(f"no-cross (ratio high)")
                else:
                    print(f"too-fast (t_cross={result['t_cross_gyr']:.4f})")

                if run_count % 10 == 0:
                    pd.DataFrame([{k: v for k, v in r.items() if k != 'checkpoints'}
                                  for r in all_results]).to_csv(summary_csv, index=False)

    df = pd.DataFrame([{k: v for k, v in r.items() if k != 'checkpoints'}
                       for r in all_results])
    df.to_csv(summary_csv, index=False)
    print(f"\nFinal summary -> {summary_csv}")
    print(f"Total viable points: {df.viable.sum()}/{len(df)}")
    for ml in models_to_run:
        sub = df[df.model_label == ml]
        print(f"  {ml}: {sub.viable.sum()}/{len(sub)} viable")

    with open(os.path.join(args.out_dir, 'grid_full.pkl'), 'wb') as f:
        pickle.dump(all_results, f)

    print(f"\nTotal time: {time.time() - t_start:.1f}s")


if __name__ == '__main__':
    main()
