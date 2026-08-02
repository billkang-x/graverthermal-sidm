"""
P2: Run dissipative halo evolution for elastic + dissipative benchmarks.

Uses parameters validated to reproduce Schmidt et al. 2026 Fig 2:
  - NFW: r_s=3.6 kpc, rho_0=7.09e-3 Msun/pc^3
  - sigma_T/m = 50 cm^2/g, w = 100 km/s
  - t_epsilon = 1e-2, n_shells = 100

Runs:
  1. Elastic control (sigma_m=50, r_diss=1.0)
  2. Const r_diss=1.05 (Schmidt benchmark)
  3. M1: Dark photon massive (calibrated, then Born-masked)
  4. M2: Scalar phi massive (calibrated, then Born-masked)
  5. M3: Massless control (calibrated, r_diss~1.05)

Each run uses a custom evolution loop with NaN checking and step-based saving.
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
from astropy import units as ut
from astropy import constants as ct

from dsidm_models import (
    benchmark_models,
    sigma_T_born,
    r_diss,
    DSIDMParameters,
    is_born_valid,
)
from emission_kernel import microscopic_cooling_sigma_m
from thermal_avg import effective_sigma_m_and_rdiss, effective_cooling_sigma_m_from_kernel
from dissipative_halo import DissipativeHalo
from SourcePy.evolve import Halo
from SourcePy.record import HaloRecord


# ----------------------------------------------------------------------
# Parameters (validated)
# ----------------------------------------------------------------------
RS_KPC = 3.6
RHO0 = 7.09e-3  # Msun/pc^3
SIGMA_M = 50.0  # cm^2/g
W_UNITS = 100.0  # km/s
N_SHELLS = 100
T_EPSILON = 1e-2
R_EPSILON = 1e-12
T_END = 150.0  # dimensionless; enough for core+collapse onset
RHO_FACTOR_END = 1000.0
MAX_STEPS = 500000
SAVE_EVERY = 500


# ----------------------------------------------------------------------
# Custom evolution loop with NaN checking
# ----------------------------------------------------------------------
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

        # NaN check
        if np.any(np.isnan(haloevo.rho)) or np.any(np.isnan(haloevo.r)):
            print(f"  [{label}] NaN at step {step}, t={haloevo.t:.4f}")
            break

        current_rho = haloevo.get_central_quantity(haloevo.rho)
        if haloevo.t >= t_end:
            haloevo.save_halo()
            print(f"  [{label}] Done: t={haloevo.t:.2f} ({(haloevo.t*haloevo.scale_t).to('Gyr').value:.2f} Gyr)")
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


# ----------------------------------------------------------------------
# Extract observables
# ----------------------------------------------------------------------
def extract_observables(halorec, haloevo, radii_kpc=None):
    if radii_kpc is None:
        radii_kpc = np.array([0.2, 0.8, 3.2])

    list_files, list_times = halorec.glob_pickle_files()
    if len(list_files) == 0:
        return None

    scale_r = haloevo.scale_r.to('kpc').value
    scale_rho = haloevo.scale_rho.to('Msun/pc**3').value
    scale_v = haloevo.scale_v.to('km/s').value
    scale_t = haloevo.scale_t.to('Gyr').value

    times_gyr, rho_center, v_center = [], [], []
    profiles_rho, profiles_v, profiles_r = [], [], []
    gamma_2D = {r: [] for r in radii_kpc}

    for f in list_files:
        try:
            data = halorec.get_halo_state_pickled(file_halo=f)
            if not data or np.any(np.isnan(data.get('rho', [np.nan]))):
                continue
        except:
            continue

        times_gyr.append(data['t'] * scale_t)
        rho_center.append(data['rho'][3] * scale_rho)
        v_center.append(np.sqrt(data['p'][3]/data['rho'][3]) * scale_v)
        profiles_rho.append(data['rho'])
        profiles_v.append(np.sqrt(data['p']/data['rho']))
        profiles_r.append(data['r'])

        # gamma_2D
        r_kpc = data['r'] * scale_r
        rho_phys = data['rho'] * scale_rho
        for r_target in radii_kpc:
            idx = np.searchsorted(r_kpc, r_target)
            if idx <= 0 or idx >= len(r_kpc):
                gamma_2D[r_target].append(np.nan)
                continue
            r1, r2 = r_kpc[idx-1], r_kpc[idx]
            rho1, rho2 = rho_phys[idx-1], rho_phys[idx]
            if rho1 <= 0 or rho2 <= 0 or r1 <= 0:
                gamma_2D[r_target].append(np.nan)
                continue
            gamma_2D[r_target].append(-np.log(rho2/rho1) / np.log(r2/r1))

    return {
        'times_gyr': np.array(times_gyr),
        'rho_center': np.array(rho_center),
        'v_center': np.array(v_center),
        'profiles_rho': np.array(profiles_rho),
        'profiles_v': np.array(profiles_v),
        'profiles_r': np.array(profiles_r),
        'gamma_2D': {r: np.array(gamma_2D[r]) for r in radii_kpc},
        'radii_kpc': radii_kpc,
        'scale_r': scale_r, 'scale_rho': scale_rho,
        'scale_v': scale_v, 'scale_t': scale_t,
    }


# ----------------------------------------------------------------------
# Calibrate DSIDM model to target sigma_T/m at v_ref
# ----------------------------------------------------------------------
def calibrate_model(p, target=50.0, v_ref=100.0):
    sig = sigma_T_born(np.array([v_ref]), p)[0]
    ratio = target / sig
    import dataclasses
    return dataclasses.replace(p, alpha_D=p.alpha_D * np.sqrt(ratio))


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    dir_output = os.path.join(_PROJ_ROOT, 'data', 'P2_runs')
    os.makedirs(dir_output, exist_ok=True)

    print("=" * 80)
    print("P2: Dissipative halo evolution")
    print(f"NFW: r_s={RS_KPC}, rho_0={RHO0}, sigma_m={SIGMA_M}, w={W_UNITS}")
    print("=" * 80)

    results = {}

    # Helper to check if a run is already done (has snapshots)
    def run_done(name):
        dir_data = os.path.join(dir_output, name)
        if not os.path.exists(dir_data):
            return False
        pickles = [f for f in os.listdir(dir_data) if f.endswith('.pickle')]
        return len(pickles) > 10  # need at least 10 snapshots

    # Helper to load results from existing run
    def load_existing(name):
        dir_data = os.path.join(dir_output, name)
        rec = HaloRecord(dir_data)
        # Create a dummy halo to get scales
        evo = Halo(rec, profile='NFW', r_s=RS_KPC, rho_s=RHO0,
                    sigma_m_with_units=SIGMA_M, w_units=W_UNITS,
                    n_shells=N_SHELLS, r_max=50.0, r_min=0.02,
                    flag_hydrostatic_initial=False)
        return extract_observables(rec, evo)

    # ---- 1. Elastic control ----
    print("\n--- Elastic control (sigma_m=50, r_diss=1.0) ---")
    name = 'elastic'
    if run_done(name):
        print(f"  [{name}] Already done, loading existing results...")
        results[name] = load_existing(name)
    else:
        dir_data = os.path.join(dir_output, name)
        if os.path.exists(dir_data):
            shutil.rmtree(dir_data)
        rec = HaloRecord(dir_data)
        evo = Halo(rec, profile='NFW', r_s=RS_KPC, rho_s=RHO0,
                    sigma_m_with_units=SIGMA_M, w_units=W_UNITS,
                    n_shells=N_SHELLS, r_max=50.0, r_min=0.02,
                    flag_hydrostatic_initial=True,
                    flag_timestep_use_relaxation=True,
                    flag_timestep_use_energy=True)
        evolve_safe(evo, label=name)
        results[name] = extract_observables(rec, evo)

    # ---- 2. Constant r_diss=1.05 ----
    print("\n--- Const r_diss=1.05 (sigma_m=50) ---")
    name = 'const_rdiss_1p05'
    if run_done(name):
        print(f"  [{name}] Already done, loading existing results...")
        results[name] = load_existing(name)
    else:
        dir_data = os.path.join(dir_output, name)
        if os.path.exists(dir_data):
            shutil.rmtree(dir_data)
        rec = HaloRecord(dir_data)
        sig_fn = lambda T: np.full(np.atleast_1d(T).shape, SIGMA_M, dtype=float)
        rd_fn = lambda T: np.full(np.atleast_1d(T).shape, 1.05, dtype=float)
        evo = DissipativeHalo(rec,
                              sigma_m_eff_callable=sig_fn,
                              rdiss_eff_callable=rd_fn,
                              cooling_sigma_m_eff_callable=lambda T: np.full(
                                  np.atleast_1d(T).shape, SIGMA_M * (1.05 - 1.0), dtype=float
                              ),
                              flag_dissipation=True,
                              profile='NFW', r_s=RS_KPC, rho_s=RHO0,
                              sigma_m_with_units=SIGMA_M, w_units=W_UNITS,
                              n_shells=N_SHELLS, r_max=50.0, r_min=0.02,
                              flag_hydrostatic_initial=True,
                              flag_timestep_use_relaxation=True,
                              flag_timestep_use_energy=True)
        evolve_safe(evo, label=name)
        results[name] = extract_observables(rec, evo)

    # ---- 3-5. Velocity-dependent models ----
    # Skip M3_massless_control: Rutherford divergence makes it numerically unstable.
    # The const_rdiss_1p05 above already serves as the "massless emission" control
    # (constant r_diss, rescaling symmetry holds).
    models = {k: v for k, v in benchmark_models().items()
              if k != 'M3_massless_control'}
    for mname, p in models.items():
        print(f"\n--- {mname} ---")
        safe_name = mname.replace('-', '_')
        p_cal = calibrate_model(p, target=SIGMA_M, v_ref=100.0)
        print(f"  alpha_D={p_cal.alpha_D:.4e}, m_med={p_cal.m_mediator:.4e}")
        if not is_born_valid(p_cal):
            print("  Skipping quantitative run: calibrated point is outside the Born mask.")
            results[mname] = None
            continue
        if run_done(safe_name):
            print(f"  [{safe_name}] Already done, loading existing results...")
            results[mname] = load_existing(safe_name)
            continue

        sigma_m_eff, rdiss_eff, _, _ = effective_sigma_m_and_rdiss(
            lambda v: sigma_T_born(np.atleast_1d(v), p_cal),
            lambda v: r_diss(np.atleast_1d(v), p_cal),
            np.logspace(1.5, 6.5, 40)
        )
        cooling_sigma_eff, _ = effective_cooling_sigma_m_from_kernel(
            lambda v: microscopic_cooling_sigma_m(v, p_cal),
            np.logspace(1.5, 6.5, 40),
        )

        dir_data = os.path.join(dir_output, safe_name)
        if os.path.exists(dir_data):
            shutil.rmtree(dir_data)
        rec = HaloRecord(dir_data)
        evo = DissipativeHalo(rec,
                              sigma_m_eff_callable=sigma_m_eff,
                              rdiss_eff_callable=rdiss_eff,
                              cooling_sigma_m_eff_callable=cooling_sigma_eff,
                              flag_dissipation=True,
                              profile='NFW', r_s=RS_KPC, rho_s=RHO0,
                              sigma_m_with_units=SIGMA_M, w_units=W_UNITS,
                              n_shells=N_SHELLS, r_max=50.0, r_min=0.02,
                              flag_hydrostatic_initial=True,
                              flag_timestep_use_relaxation=True,
                              flag_timestep_use_energy=True)
        evolve_safe(evo, label=safe_name)
        results[mname] = extract_observables(rec, evo)

    # ---- Save CSV ----
    print("\n" + "=" * 80)
    print("Saving results...")
    import csv
    csv_path = os.path.join(_PROJ_ROOT, 'data', 'P2_collapse_times.csv')
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['model', 't_core_gyr', 'rho_min', 't_final_gyr', 'rho_final', 'v_final'])
        for name, obs in results.items():
            if obs is None or len(obs['times_gyr']) == 0:
                continue
            idx_min = np.argmin(obs['rho_center'])
            w.writerow([name,
                        f'{obs["times_gyr"][idx_min]:.4f}',
                        f'{obs["rho_center"][idx_min]:.4e}',
                        f'{obs["times_gyr"][-1]:.4f}',
                        f'{obs["rho_center"][-1]:.4e}',
                        f'{obs["v_center"][-1]:.2f}'])
    print(f"  Saved: {csv_path}")

    # ---- Plot ----
    plot_results(results, _PROJ_ROOT)
    return results


# ----------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------
def plot_results(results, proj_root):
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    colors = {'elastic': 'b', 'const_rdiss_1p05': 'k',
              'M1_dark_photon_massive': 'r',
              'M2_scalar_phi_massive': 'g',
              'M3_massless_control': 'm'}
    labels = {'elastic': 'Elastic (r_diss=1.0)',
              'const_rdiss_1p05': 'Const r_diss=1.05',
              'M1_dark_photon_massive': 'M1: Dark photon m_V=1.1 keV',
              'M2_scalar_phi_massive': 'M2: Scalar m_phi=6.95 keV',
              'M3_massless_control': 'M3: Massless (control)'}

    # (a) Central density
    ax = axes[0, 0]
    for name, obs in results.items():
        if obs is None:
            continue
        ax.plot(obs['times_gyr'], obs['rho_center'],
                color=colors.get(name, 'c'), label=labels.get(name, name), lw=1.5)
    ax.set_xlabel('Time [Gyr]')
    ax.set_ylabel(r'$\rho_{\rm center}$ [$M_\odot$/pc$^3$]')
    ax.set_yscale('log')
    ax.set_title('Central density evolution')
    ax.legend(fontsize=8)
    ax.set_xlim(0, 10)
    ax.grid(True, alpha=0.3)

    # (b) Central velocity dispersion
    ax = axes[0, 1]
    for name, obs in results.items():
        if obs is None:
            continue
        ax.plot(obs['times_gyr'], obs['v_center'],
                color=colors.get(name, 'c'), label=labels.get(name, name), lw=1.5)
    ax.set_xlabel('Time [Gyr]')
    ax.set_ylabel(r'$\nu_{\rm center}$ [km/s]')
    ax.set_title('Central velocity dispersion')
    ax.legend(fontsize=8)
    ax.set_xlim(0, 10)
    ax.grid(True, alpha=0.3)

    # (c) gamma_2D at r=0.2 kpc
    ax = axes[1, 0]
    r_ref = 0.2
    for name, obs in results.items():
        if obs is None:
            continue
        ax.plot(obs['times_gyr'], obs['gamma_2D'][r_ref],
                color=colors.get(name, 'c'), label=labels.get(name, name), lw=1.5)
    ax.set_xlabel('Time [Gyr]')
    ax.set_ylabel(r'$\gamma(r=0.2\,{\rm kpc})$')
    ax.set_title('Inner density slope')
    ax.legend(fontsize=8)
    ax.set_xlim(0, 10)
    ax.grid(True, alpha=0.3)

    # (d) Final density profiles
    ax = axes[1, 1]
    for name, obs in results.items():
        if obs is None or len(obs['profiles_rho']) == 0:
            continue
        r_kpc = obs['profiles_r'][-1] * obs['scale_r']
        rho = obs['profiles_rho'][-1] * obs['scale_rho']
        idx = np.argsort(r_kpc)
        ax.plot(r_kpc[idx], rho[idx], color=colors.get(name, 'c'),
                label=labels.get(name, name), lw=1.5)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('r [kpc]')
    ax.set_ylabel(r'$\rho$ [$M_\odot$/pc$^3$]')
    ax.set_title('Final density profiles')
    ax.legend(fontsize=8)
    ax.set_xlim(0.01, 50)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig_path = os.path.join(proj_root, 'figures', 'P2_evolution.png')
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    print(f"  Saved: {fig_path}")
    plt.close()


if __name__ == '__main__':
    main()
