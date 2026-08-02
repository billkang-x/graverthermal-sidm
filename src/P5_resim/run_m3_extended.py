"""
Re-run M3 grid scan points with extended evolution limits.

The original grid scan used MAX_STEPS=30000, wall_budget_s=30, and
checkpoints up to t=40 (dimless) ~ 0.053 Gyr. For M3 (const r_diss=1.05),
the cooling is weak enough that the mass ratio continues decreasing
slowly well beyond 0.053 Gyr. This script re-runs M3 points that:
  - Were not viable in the original scan
  - Had final_ratio > MASS_RATIO_OBS (ratio still above observed)
  - Were not skipped as init_below_obs

with extended limits: MAX_STEPS=200000, wall_budget_s=120,
checkpoints up to t=5120 (dimless) ~ 6.8 Gyr (covers t_obs=6.37).
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
sys.path.insert(0, os.path.join(_PROJ_ROOT, 'src', 'P3_rescaling'))

import matplotlib
matplotlib.use('Agg')

from astropy import units as ut
from astropy import constants as ct

from dissipative_halo import DissipativeHalo
from SourcePy.evolve import Halo
from SourcePy.record import HaloRecord
from rescale import (
    projected_enclosed_mass,
    M_20PC, M_90PC, MASS_RATIO_OBS, MASS_RATIO_ERR,
)

# Extended constants for M3
T_OBS_GYR = 6.37
RHO_0_FIXED = 10.0
N_SHELLS = 100
T_EPSILON = 1e-2
R_EPSILON = 1e-12
RHO_FACTOR_END = 1000.0
MAX_STEPS = 200000  # Extended: 200k steps
SAVE_EVERY = 2000
W_UNITS = 100.0

# Extended checkpoints: up to t=5120 dimless ~ 6.8 Gyr
T_CHECKPOINTS = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0,
                 80.0, 160.0, 320.0, 640.0, 1280.0, 2560.0, 5120.0]

WALL_BUDGET_S = 120  # 2 minutes per point


def evolve_extended(haloevo, checkpoints, label):
    """Extended evolution with more checkpoints and longer wall budget."""
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

    for step in range(MAX_STEPS):
        if time.time() - t_start_wall > WALL_BUDGET_S:
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
                if len(results) >= 2:
                    d_prev = results[-2]['mass_ratio'] - MASS_RATIO_OBS
                    d_curr = results[-1]['mass_ratio'] - MASS_RATIO_OBS
                    if d_prev * d_curr < 0:
                        found_cross = True
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


def run_m3_point(sigma_m, r_s_kpc, out_dir, rho_0=RHO_0_FIXED):
    """Run one M3 grid point with extended limits."""
    label = f"M3_ext_sig{sigma_m:.4e}_rs{r_s_kpc:.4f}"
    init_ratio = None

    # Compute init ratio
    r_arr = np.logspace(-2, 1, 200)
    x = r_arr / r_s_kpc
    rho_arr = rho_0 / (x * (1.0 + x)**2)
    M_in = projected_enclosed_mass(r_arr, rho_arr, 0.02, r_unit='kpc', rho_unit='Msun_pc3')
    M_out = projected_enclosed_mass(r_arr, rho_arr, 0.09, r_unit='kpc', rho_unit='Msun_pc3')
    init_ratio = M_in / M_out if M_out > 0 else np.nan

    out_subdir = os.path.join(out_dir, '_tmp', label)
    os.makedirs(out_subdir, exist_ok=True)
    rec = HaloRecord(out_subdir)

    try:
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

        results = evolve_extended(evo, T_CHECKPOINTS, label=label)

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
            'model_label': 'M3', 'model_key': 'const_rdiss_1p05',
            'sigma_m_cm2_g': sigma_m, 'r_s_kpc': r_s_kpc,
            'rho_0_msun_pc3': rho_0,
            't_cross_gyr': t_cross_gyr,
            'viable': bool(not np.isnan(t_cross_gyr) and t_cross_gyr <= T_OBS_GYR),
            'n_checkpoints': len(results),
            'final_ratio': results[-1]['mass_ratio'] if results else np.nan,
            'final_t_gyr': results[-1]['t_gyr'] if results else np.nan,
            'init_ratio': init_ratio,
            'skipped': 'no',
            'extended': True,
        }
        if results:
            result['checkpoints'] = results

    except Exception as e:
        result = {
            'model_label': 'M3', 'model_key': 'const_rdiss_1p05',
            'sigma_m_cm2_g': sigma_m, 'r_s_kpc': r_s_kpc,
            'rho_0_msun_pc3': rho_0,
            't_cross_gyr': np.nan, 'viable': False,
            'n_checkpoints': 0, 'final_ratio': np.nan, 'final_t_gyr': np.nan,
            'init_ratio': init_ratio, 'skipped': f'error: {e}',
            'extended': True,
        }
    finally:
        import shutil
        shutil.rmtree(out_subdir, ignore_errors=True)

    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('out_dir', help='Output directory')
    parser.add_argument('--max_points', type=int, default=60,
                        help='Max number of points to re-run (default 60)')
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(os.path.join(args.out_dir, '_tmp'), exist_ok=True)

    # Load original grid scan results
    df_orig = pd.read_csv('data/P5_grid_scan/grid_summary.csv')

    # Select M3 points to re-run: not viable, not skipped, ratio > obs
    m3_nv = df_orig[(df_orig.model_label == 'M3') &
                    (~df_orig.viable) &
                    (df_orig.skipped == 'no') &
                    (df_orig.final_ratio > MASS_RATIO_OBS)]
    print(f"M3 non-viable points with ratio > obs: {len(m3_nv)}")
    print(f"  sigma/m range: [{m3_nv.sigma_m_cm2_g.min():.4f}, {m3_nv.sigma_m_cm2_g.max():.4f}]")
    print(f"  r_s range:     [{m3_nv.r_s_kpc.min():.4f}, {m3_nv.r_s_kpc.max():.4f}]")

    # Prioritize: points closest to observed ratio (most likely to cross with more time)
    m3_nv = m3_nv.assign(ratio_diff=m3_nv.final_ratio - MASS_RATIO_OBS)
    m3_nv = m3_nv.sort_values('ratio_diff').head(args.max_points)

    print(f"\nRe-running {len(m3_nv)} M3 points with extended limits...")
    print(f"  MAX_STEPS={MAX_STEPS}, WALL_BUDGET={WALL_BUDGET_S}s, max_ckpt=5120 (dimless)")

    all_results = []
    summary_csv = os.path.join(args.out_dir, 'm3_extended_summary.csv')
    t_start = time.time()

    for i, (_, row) in enumerate(m3_nv.iterrows()):
        sig = row['sigma_m_cm2_g']
        rs = row['r_s_kpc']
        elapsed = time.time() - t_start
        print(f"[{i+1}/{len(m3_nv)}] sig={sig:.4f} rs={rs:.4f} "
              f"({elapsed:.0f}s)", end=' ', flush=True)

        result = run_m3_point(sig, rs, args.out_dir)
        all_results.append(result)

        if result.get('viable'):
            print(f"VIABLE t_cross={result['t_cross_gyr']:.4f} Gyr")
        elif result['t_cross_gyr'] == np.inf:
            print(f"no-cross (final={result['final_ratio']:.4f} at t={result['final_t_gyr']:.4f})")
        else:
            print(f"too-fast t_cross={result['t_cross_gyr']:.4f}")

        if (i + 1) % 5 == 0:
            pd.DataFrame([{k: v for k, v in r.items() if k != 'checkpoints'}
                          for r in all_results]).to_csv(summary_csv, index=False)

    df_out = pd.DataFrame([{k: v for k, v in r.items() if k != 'checkpoints'}
                           for r in all_results])
    df_out.to_csv(summary_csv, index=False)
    print(f"\nFinal summary -> {summary_csv}")
    print(f"Viable: {df_out.viable.sum()}/{len(df_out)}")
    print(f"Total time: {time.time() - t_start:.1f}s")

    with open(os.path.join(args.out_dir, 'm3_extended_full.pkl'), 'wb') as f:
        pickle.dump(all_results, f)


if __name__ == '__main__':
    main()
