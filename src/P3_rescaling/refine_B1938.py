"""
B1938+666 refinement: trace the full viable region in (σ/m, t_evo) using
the rescaling symmetry.

For velocity-INDEPENDENT models (elastic, const r_diss), the rescaling
symmetry of Appendix G means:

  The simulated halo at (r_s_sim, ρ_0_sim, σ_sim, t_sim) with mass ratio
  matching the observation at r_2D/r_s = x_0 maps to a physical halo with

      r_s_phys = λ(x_0) r_s_sim
      ρ_0_phys = μ(x_0) ρ_0_sim / λ³
      σ_phys   = (λ²/μ) σ_sim
      t_evo    = √(λ³/μ) t_sim

  Crucially, for a FIXED matching point (x_0) on a FIXED snapshot,
  varying σ_sim traces a curve σ_phys(σ_sim) = (λ²/μ) σ_sim that is
  LINEAR in σ_sim. So a single simulation traces the full family.

However, the snapshot TIME also enters: t_sim corresponds to the physical
evolution time t_evo = √(λ³/μ) t_sim. As we scan σ_sim, we get a family
of physical halos with the SAME mass ratio but DIFFERENT (σ_phys, t_evo).

For velocity-DEPENDENT models (M1, M2), the symmetry is broken — we cannot
analytically extend. But we can still scan α_D to generate model curves
through the (σ_phys, t_evo) plane, and overlay the B1938 viable region
extracted from the elastic case.

This module:
  1. Loads the elastic P2 run
  2. For each snapshot × each matching r_2D/rs, computes (λ, μ)
  3. For each matching point, traces the family {σ_phys, t_evo} over σ_sim
     in [0.01, 10000] cm²/g
  4. Filters by t_evo ≤ 6.37 Gyr
  5. Produces the viable region envelope in (σ_phys, t_evo)
"""
from __future__ import annotations

import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy import units as ut
from astropy import constants as ct

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJ_ROOT = os.path.normpath(os.path.join(_HERE, '..', '..'))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_PROJ_ROOT, 'external', 'gravothermalsidm'))

from rescale import (
    find_matching_radii, compute_rescaling,
    M_20PC, M_90PC, MASS_RATIO_OBS, MASS_RATIO_ERR,
    R_INNER_PC, R_OUTER_PC, R_RATIO, RS_SIM_KPC, T_ZOBS_GYR,
)
from SourcePy.record import HaloRecord

OUT_FIG = os.path.join(_PROJ_ROOT, 'figures')
OUT_DATA = os.path.join(_PROJ_ROOT, 'data')


def get_scales(halo_ini):
    r_s = halo_ini['r_s'] * ut.kpc
    rho_s = halo_ini['rho_s'] * ut.M_sun / ut.pc**3
    scale_r_kpc = r_s.to('kpc').value
    scale_rho = rho_s.to('Msun/pc**3').value
    scale_t = (1.0 / np.sqrt(4.0 * np.pi * rho_s * ct.G)).to('Gyr').value
    return scale_r_kpc, scale_rho, scale_t


def trace_viable_region(model_dir='elastic', sigma_sim_arr=None,
                        n_scan_rs=80, snapshot_stride=5):
    """Trace the full B1938 viable region using the elastic rescaling symmetry.

    For each (snapshot, r_2D/rs) match point, compute (λ, μ) from the
    simulation mass ratio. Then for each σ_sim value, the physical
    parameters are:
        σ_phys = (λ²/μ) σ_sim
        t_evo  = √(λ³/μ) t_sim

    We collect all (σ_phys, t_evo) pairs that satisfy t_evo ≤ T_ZOBS.

    Returns DataFrame.
    """
    if sigma_sim_arr is None:
        sigma_sim_arr = np.logspace(-2, 4, 50)  # 0.01 to 10000 cm²/g

    dir_data = os.path.join(_PROJ_ROOT, 'data', 'P2_runs', model_dir)
    if not os.path.isdir(dir_data):
        print(f"Skip: {dir_data} not found")
        return pd.DataFrame()

    halorec = HaloRecord(dir_data)
    list_files, _ = halorec.glob_pickle_files()
    halo_ini, _ = halorec.get_halo_initialization()
    scale_r_kpc, scale_rho, scale_t = get_scales(halo_ini)

    rows = []
    files_iter = list_files[::snapshot_stride]
    for i, f in enumerate(files_iter):
        try:
            data = halorec.get_halo_state_pickled(file_halo=f)
            if np.any(np.isnan(data.get('rho', [np.nan]))):
                continue
        except Exception:
            continue

        r_kpc = data['r'] * scale_r_kpc
        rho = data['rho'] * scale_rho
        t_sim_gyr = data['t'] * scale_t

        scan = find_matching_radii(r_kpc, rho, RS_SIM_KPC, n_scan=n_scan_rs)
        for r2D_rs, M_inner, M_outer, ratio in scan:
            if not (M_inner > 0 and M_outer > 0 and ratio > 0):
                continue
            n_sigma = abs(ratio - MASS_RATIO_OBS) / MASS_RATIO_ERR
            if n_sigma > 3:
                continue

            lam, mu = compute_rescaling(r2D_rs, M_inner, M_outer)
            if not (lam > 0 and mu > 0):
                continue

            # σ_phys = (λ²/μ) σ_sim;  t_evo = √(λ³/μ) t_sim
            sigma_phys_arr = (lam**2 / mu) * sigma_sim_arr
            # t_evo does not depend on σ_sim (under rescaling symmetry),
            # so it's a scalar for this (snapshot, r_2D/rs) match.
            t_evo_scalar = float(np.sqrt(lam**3 / mu) * t_sim_gyr)

            for s_sim, s_phys in zip(sigma_sim_arr, sigma_phys_arr):
                rows.append({
                    'model': model_dir,
                    'snapshot_idx': i,
                    'snapshot_time_gyr': t_sim_gyr,
                    'r2D_rs': r2D_rs,
                    'n_sigma': n_sigma,
                    'lambda': lam,
                    'mu': mu,
                    'sigma_sim': s_sim,
                    'sigma_phys_cm2_g': s_phys,
                    't_evo_gyr': t_evo_scalar,
                    't_ok': t_evo_scalar <= T_ZOBS_GYR,
                })

    return pd.DataFrame(rows)


def plot_viable_region(df, out_path):
    """Plot the viable region in (σ_phys, t_evo) colored by r_2D/rs."""
    fig, ax = plt.subplots(figsize=(9, 7))
    viable = df[df['t_ok']]
    excluded = df[~df['t_ok']]

    # Color by r_2D/rs
    if len(viable):
        sc = ax.scatter(viable['sigma_phys_cm2_g'], viable['t_evo_gyr'],
                        c=viable['r2D_rs'], s=3, cmap='viridis', alpha=0.5,
                        edgecolors='none')
        cb = fig.colorbar(sc, ax=ax, pad=0.02)
        cb.set_label(r'$r_{2D}/r_s$', fontsize=11)
    if len(excluded):
        ax.scatter(excluded['sigma_phys_cm2_g'], excluded['t_evo_gyr'],
                   s=2, c='lightgray', marker='x', alpha=0.2,
                   label=r'$t_{\rm evo} > 6.37$ Gyr (excluded)')

    ax.axhline(T_ZOBS_GYR, color='red', linestyle='--', lw=1.2,
               label=r'$t_{\rm zobs} = 6.37$ Gyr')
    # Shade Schmidt Table 2 benchmark points
    benchmarks = [(0.07, 0.41, 'Row #3'),
                  (1.0, 4.8, 'Row #2'),
                  (7.8, 28, 'Row #1 (excluded)')]
    for sm, te, lbl in benchmarks:
        ax.scatter([sm], [te], c='red', marker='*', s=120, zorder=10,
                   edgecolors='black', linewidths=0.8)
        ax.annotate(lbl, (sm, te), textcoords='offset points',
                     xytext=(8, 5), fontsize=9, color='red')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$\sigma_T / m_\chi$ [cm$^2$/g] (physical)', fontsize=12)
    ax.set_ylabel(r'$t_{\rm evo}$ [Gyr]', fontsize=12)
    ax.set_title('B1938+666 viable region (elastic rescaling symmetry)',
                 fontsize=12)
    ax.legend(loc='lower right', fontsize=10)
    ax.grid(True, which='both', alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    print(f"Saved {out_path}")


def plot_viable_region_envelope(df, out_path):
    """Plot the envelope (upper and lower t_evo as function of σ_phys)."""
    viable = df[df['t_ok']]
    if len(viable) == 0:
        print("No viable points")
        return

    # Bin in σ_phys and compute envelope
    log_s = np.log10(viable['sigma_phys_cm2_g'].values)
    log_t = np.log10(viable['t_evo_gyr'].values)
    bins = np.linspace(log_s.min(), log_s.max(), 40)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    t_upper = np.zeros(len(bins) - 1)
    t_lower = np.zeros(len(bins) - 1)
    for i in range(len(bins) - 1):
        mask = (log_s >= bins[i]) & (log_s < bins[i+1])
        if np.any(mask):
            t_upper[i] = 10**np.max(log_t[mask])
            t_lower[i] = 10**np.min(log_t[mask])
        else:
            t_upper[i] = np.nan
            t_lower[i] = np.nan

    fig, ax = plt.subplots(figsize=(9, 7))
    valid = np.isfinite(t_upper)
    ax.fill_between(10**bin_centers[valid], t_lower[valid], t_upper[valid],
                    color='tab:blue', alpha=0.3,
                    label='Viable region (elastic)')
    ax.plot(10**bin_centers[valid], t_upper[valid], color='tab:blue', lw=2,
            label='Upper envelope (slowest evolution)')
    ax.plot(10**bin_centers[valid], t_lower[valid], color='tab:blue', lw=2,
            linestyle='--', label='Lower envelope (fastest evolution)')

    # Benchmarks
    benchmarks = [(0.07, 0.41, 'Row #3'), (1.0, 4.8, 'Row #2')]
    for sm, te, lbl in benchmarks:
        ax.scatter([sm], [te], c='red', marker='*', s=150, zorder=10,
                   edgecolors='black', linewidths=0.8)
        ax.annotate(lbl, (sm, te), textcoords='offset points',
                     xytext=(8, 5), fontsize=9, color='red')

    ax.axhline(T_ZOBS_GYR, color='red', linestyle='--', lw=1.2,
               label=r'$t_{\rm zobs} = 6.37$ Gyr')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$\sigma_T / m_\chi$ [cm$^2$/g] (physical)', fontsize=12)
    ax.set_ylabel(r'$t_{\rm evo}$ [Gyr]', fontsize=12)
    ax.set_title('B1938+666 viable region envelope', fontsize=12)
    ax.legend(loc='upper left', fontsize=10)
    ax.grid(True, which='both', alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    print(f"Saved {out_path}")


def main():
    print("[B1938 refine] Tracing viable region via rescaling symmetry...")
    df = trace_viable_region(model_dir='elastic', n_scan_rs=80,
                              snapshot_stride=5)
    print(f"Total points: {len(df)} (viable: {int(df['t_ok'].sum())})")
    if len(df) == 0:
        return

    csv_path = os.path.join(OUT_DATA, 'P3_B1938_viable_region.csv')
    df.to_csv(csv_path, index=False)
    print(f"Saved {csv_path}")

    # Filter to physically relevant σ range
    df_phys = df[(df['sigma_phys_cm2_g'] >= 0.001) &
                 (df['sigma_phys_cm2_g'] <= 1e4)]
    print(f"In physical σ range: {len(df_phys)} "
          f"(viable: {int(df_phys['t_ok'].sum())})")

    plot_viable_region(df_phys,
                       os.path.join(OUT_FIG, 'P3_B1938_viable_scatter.png'))
    plot_viable_region_envelope(df_phys,
                       os.path.join(OUT_FIG, 'P3_B1938_viable_envelope.png'))

    # Summary statistics
    viable = df_phys[df_phys['t_ok']]
    print(f"\n=== Viable region summary ===")
    print(f"  σ_phys range: {viable['sigma_phys_cm2_g'].min():.3e} - "
          f"{viable['sigma_phys_cm2_g'].max():.3e} cm²/g")
    print(f"  t_evo range:  {viable['t_evo_gyr'].min():.4f} - "
          f"{viable['t_evo_gyr'].max():.4f} Gyr")
    print(f"  r_2D/rs range: {viable['r2D_rs'].min():.4f} - "
          f"{viable['r2D_rs'].max():.4f}")


if __name__ == '__main__':
    main()
