"""
Concentration scan: scan c_200 in {6, 8.5, 12, 16, 20} for M1 (dark photon)
to explore whether higher initial NFW concentration reaches the extended-
solution branch and produces stronger symmetry breaking.

For each c_200 we:
  1. Set r_s and rho_0 to match the target c_200 (keeping M_200 fixed)
  2. Run the M1 dark photon dissipative halo evolution
  3. Snapshot the projected enclosed mass ratio at several t_evo checkpoints
  4. Find which (t_evo, sigma/m) points match the observed B1938 ratio

The c_200 = 8.5 case reproduces the default P2 run; c_200 = 12 reproduces
the P2_highconc run. The new cases are c_200 = 6, 16, 20.

Outputs:
  data/P5_conc_scan/summary.csv : per-concentration summary
  data/P5_conc_scan/<c200>/snapshots/ : per-run snapshots
  figures/P5_conc_scan.png       : mass ratio vs t_evo for each c_200
"""
from __future__ import annotations

import os, sys, time, argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJ_ROOT = os.path.normpath(os.path.join(_HERE, '..', '..'))
sys.path.insert(0, os.path.join(_PROJ_ROOT, 'src', 'cross_sections'))
sys.path.insert(0, os.path.join(_PROJ_ROOT, 'src', 'fluid_runner'))
sys.path.insert(0, os.path.join(_PROJ_ROOT, 'external', 'gravothermalsidm'))
sys.path.insert(0, os.path.join(_PROJ_ROOT, 'src', 'P3_rescaling'))

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
# Concentration scan parameters
# ----------------------------------------------------------------------
# NFW concentration c_200 = r_200 / r_s.
# Fix M_200 = 1e12 Msun (typical subhalo). Then:
#   r_200 = (3 M_200 / (4 pi 200 rho_crit))^(1/3)
#   r_s = r_200 / c_200
#   rho_s = M_200 / (4 pi r_s^3 * (ln(1+c) - c/(1+c)))
#
# For the default run, c_200 ~ 8.5, r_s ~ 5 kpc, rho_s ~ 0.005 Msun/pc^3.

# Default P2 parameters (c_200 ~ 8.5): r_s=5.0 kpc, rho_0=5e-3
# Highconc P2 parameters (c_200 ~ 12): r_s=2.5 kpc, rho_0=1.47e-2
# We extend to c_200 = 6, 16, 20 with the same M_200 scaling.

def c200_to_nfw_params(c_200, M200_msun=1e12):
    """Compute (r_s_kpc, rho_s_msun_pc3) from c_200 and M_200.

    Uses the NFW profile: M(<r_200) = M_200, with
    rho_s = M_200 / (4*pi*r_s^3 * f(c)) where f(c) = ln(1+c) - c/(1+c).
    """
    # rho_crit(z=0) ~ 1.36e-7 Msun/pc^3 (with h=0.7)
    rho_crit = 1.36e-7  # Msun/pc^3
    rho_200 = 200.0 * rho_crit
    r_200_pc = (3.0 * M200_msun / (4.0 * np.pi * rho_200)) ** (1.0 / 3.0)
    r_s_pc = r_200_pc / c_200
    f_c = np.log(1.0 + c_200) - c_200 / (1.0 + c_200)
    rho_s = M200_msun / (4.0 * np.pi * r_s_pc ** 3 * f_c)
    return r_s_pc / 1000.0, rho_s  # r_s in kpc, rho_s in Msun/pc^3


# Target concentrations
CONCENTRATIONS = [6.0, 8.5, 12.0, 16.0, 20.0]

# Sigma/m to use for the M1 dark photon model (calibrated at v_ref=100 km/s)
# We pick a sigma that gives a reasonable evolution timescale.
SIGMA_M_TARGET = 0.1  # cm^2/g (typical SIDM value)

# Evolution checkpoints (dimensionless time)
# We evolve to a large t and snapshot at intervals
T_CHECKPOINTS_DIMLESS = [1.0, 5.0, 10.0, 20.0, 40.0, 60.0]

N_SHELLS = 100
MAX_STEPS = 200000
SAVE_EVERY = 500
RHO_FACTOR_END = 1000.0
T_EPSILON = 1e-2
R_EPSILON = 1e-12
W_UNITS = 100.0


def evolve_with_checkpoints(haloevo, checkpoints, label, out_subdir):
    """Evolve the halo, snapshotting at each checkpoint."""
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

    t0 = time.time()
    for step in range(MAX_STEPS):
        try:
            haloevo.conduct_heat()
            haloevo.hydrostatic_adjustment()
        except Exception as e:
            print(f"  [{label}] Error at step {step}: {e}")
            break

        if np.any(np.isnan(haloevo.rho)) or np.any(np.isnan(haloevo.r)):
            print(f"  [{label}] NaN at step {step}, t={haloevo.t:.4f}")
            break

        # Check if we've reached the next checkpoint
        while (next_ckpt_idx < len(checkpoints_sorted)
               and haloevo.t >= checkpoints_sorted[next_ckpt_idx]):
            haloevo.save_halo()
            # Extract mass ratio at this checkpoint
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
                    'M_inner': float(M_in),
                    'M_outer': float(M_out),
                    'mass_ratio': float(ratio),
                    'delta_pct': 100.0 * (ratio - MASS_RATIO_OBS) / MASS_RATIO_OBS,
                })
                print(f"  [{label}] ckpt {next_ckpt_idx}: t={data['t']:.2f} "
                      f"({data['t']*haloevo.scale_t.to('Gyr').value:.4e} Gyr) "
                      f"M(20)/M(90)={ratio:.4f} Δ={100.0*(ratio-MASS_RATIO_OBS)/MASS_RATIO_OBS:+.2f}%")
            next_ckpt_idx += 1

        if next_ckpt_idx >= len(checkpoints_sorted):
            print(f"  [{label}] All checkpoints done at step {step+1}")
            break

        current_rho = haloevo.get_central_quantity(haloevo.rho)
        if current_rho > RHO_FACTOR_END * haloevo.rho_center:
            haloevo.save_halo()
            print(f"  [{label}] Collapse trigger at t={haloevo.t:.4f}")
            break

        if (step + 1) % SAVE_EVERY == 0:
            haloevo.save_halo()

    elapsed = time.time() - t0
    print(f"  [{label}] done in {elapsed:.1f}s, {next_ckpt_idx}/{len(checkpoints_sorted)} checkpoints")
    return results, elapsed


def run_concentration(c_200, out_dir, max_seconds=600):
    """Run M1 dark photon evolution for a given c_200."""
    r_s, rho_s = c200_to_nfw_params(c_200)
    print(f"\n[c_200={c_200}] r_s={r_s:.4f} kpc, rho_s={rho_s:.4e} Msun/pc^3")

    out_subdir = os.path.join(out_dir, f"c{c_200:.0f}", 'snapshots')
    os.makedirs(out_subdir, exist_ok=True)
    rec = HaloRecord(out_subdir)

    bm = benchmark_models()
    p_model = bm['M1_dark_photon_massive']

    # Calibrate alpha_D so sigma_T_born(v_ref=100 km/s) = SIGMA_M_TARGET
    sig_now = sigma_T_born(np.array([W_UNITS]), p_model)[0]
    ratio = SIGMA_M_TARGET / sig_now
    import dataclasses
    p_cal = dataclasses.replace(p_model, alpha_D=p_model.alpha_D * np.sqrt(ratio))

    sigma_m_eff, rdiss_eff, _, _ = effective_sigma_m_and_rdiss(
        lambda v: sigma_T_born(np.atleast_1d(v), p_cal),
        lambda v: r_diss(np.atleast_1d(v), p_cal),
        np.logspace(1.5, 6.5, 40),
    )

    evo = DissipativeHalo(rec,
                          sigma_m_eff_callable=sigma_m_eff,
                          rdiss_eff_callable=rdiss_eff,
                          flag_dissipation=True,
                          profile='NFW', r_s=r_s, rho_s=rho_s,
                          sigma_m_with_units=SIGMA_M_TARGET, w_units=W_UNITS,
                          n_shells=N_SHELLS, r_max=50.0, r_min=0.02,
                          flag_hydrostatic_initial=True,
                          flag_timestep_use_relaxation=True,
                          flag_timestep_use_energy=True)

    results, elapsed = evolve_with_checkpoints(
        evo, T_CHECKPOINTS_DIMLESS, label=f"c{c_200:.0f}",
        out_subdir=out_subdir)

    for r in results:
        r['c_200'] = c_200
        r['r_s_kpc'] = r_s
        r['rho_s_msun_pc3'] = rho_s
        r['sigma_m_cm2_g'] = SIGMA_M_TARGET

    # Clean up snapshots to save disk
    import shutil
    shutil.rmtree(out_subdir, ignore_errors=True)

    return results, elapsed


def plot_conc_scan(all_results, out_path):
    """Mass ratio vs t_evo for each c_200."""
    fig, ax = plt.subplots(figsize=(8, 6))

    colors = plt.cm.viridis(np.linspace(0.2, 0.8, len(CONCENTRATIONS)))
    for i, c in enumerate(CONCENTRATIONS):
        sub = [r for r in all_results if r['c_200'] == c]
        if not sub:
            continue
        df = pd.DataFrame(sub).sort_values('t_dimless')
        ax.plot(df['t_gyr'], df['mass_ratio'], 'o-',
                color=colors[i], label=f'$c_{{200}} = {c:.0f}$', markersize=6)

    ax.axhline(MASS_RATIO_OBS, color='k', lw=1.5, ls='-',
               label=f'Observed ({MASS_RATIO_OBS:.3f})')
    ax.axhspan(MASS_RATIO_OBS - MASS_RATIO_ERR,
               MASS_RATIO_OBS + MASS_RATIO_ERR,
               color='gray', alpha=0.2, label=r'$\pm 1\sigma$')

    ax.set_xlabel(r'$t_{\rm evo}$ [Gyr]', fontsize=12)
    ax.set_ylabel(r'$M(r<20\,{\rm pc})/M(r<90\,{\rm pc})$', fontsize=12)
    ax.set_title(f'P5 concentration scan: M1 dark photon '
                 f'($\\sigma_T/m = {SIGMA_M_TARGET}$ cm$^2$/g)', fontsize=11)
    ax.set_xscale('log')
    ax.legend(loc='best', fontsize=9, framealpha=0.9)
    ax.grid(True, alpha=0.3, which='both')

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {out_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('out_dir', help='Output directory')
    parser.add_argument('--max_seconds', type=int, default=600,
                        help='Per-run wall-clock budget (default 600s)')
    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    all_results = []
    summary_csv = os.path.join(args.out_dir, 'summary.csv')

    for c in CONCENTRATIONS:
        try:
            results, elapsed = run_concentration(c, args.out_dir,
                                                 max_seconds=args.max_seconds)
            all_results.extend(results)
            print(f"  c={c}: {len(results)} checkpoints, {elapsed:.1f}s")
        except Exception as e:
            import traceback
            print(f"  ERROR for c={c}: {e}")
            traceback.print_exc()

        # Save progress
        if all_results:
            pd.DataFrame(all_results).to_csv(summary_csv, index=False)

    if all_results:
        df = pd.DataFrame(all_results)
        df.to_csv(summary_csv, index=False)
        print(f"\nFinal summary → {summary_csv}")
        print(df[['c_200', 't_gyr', 'mass_ratio', 'delta_pct']].to_string(index=False))

        out_plot = os.path.join(_PROJ_ROOT, 'figures', 'P5_conc_scan.png')
        plot_conc_scan(all_results, out_plot)


if __name__ == '__main__':
    main()
