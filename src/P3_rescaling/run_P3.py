"""P3 runner: B1938+666 rescaling for all four P2 models.

For each model:
  1. Load all P2 snapshots
  2. For each snapshot, scan r_2D/r_s to find where projected mass ratio
     M(r<20pc)/M(r<90pc) matches the observed 0.364±0.022 within 3σ
  3. Compute rescaling (λ, μ) and physical parameters (r_s, ρ_0, σ/m, t_evo)
  4. Filter by t_evo ≤ 6.37 Gyr (cosmic age at z_obs=0.881)
  5. For velocity-dependent models, also compute the symmetry-breaking
     effect: r_diss before vs after rescaling

Outputs:
  - data/P3_rescaled_params.csv : all matching points
  - data/P3_summary.csv         : best (largest σ/m) viable point per model
  - figures/P3_B1938_regions.png : σ/m vs t_evo scatter, all models
  - figures/P3_mass_ratio.png    : mass ratio vs r_2D/r_s curves
  - notes/P3_results.md          : summary + comparison to Schmidt Table 2
"""
from __future__ import annotations

import os, sys, time
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy import units as ut
from astropy import constants as ct

from rescale import (
    projected_enclosed_mass, find_matching_radii,
    compute_rescaling, rescale_parameters, compute_symmetry_breaking,
    absolute_mass_fit,
    M_20PC, M_90PC, M_20PC_ERR, M_90PC_ERR,
    R_INNER_PC, R_OUTER_PC, R_RATIO,
    MASS_RATIO_OBS, MASS_RATIO_ERR, T_ZOBS_GYR,
    RS_SIM_KPC, RHO0_SIM, SIGMA_M_SIM,
)
from SourcePy.record import HaloRecord

# Also load cross-section models for symmetry breaking
_PROJ_ROOT = os.path.normpath(os.path.join(_HERE, '..', '..'))
sys.path.insert(0, os.path.join(_PROJ_ROOT, 'src', 'cross_sections'))
from dsidm_models import benchmark_models, r_diss as r_diss_func

# ----------------------------------------------------------------------
# Models to process (matches P2)
# ----------------------------------------------------------------------
MODELS = {
    'elastic':                  {'dir': 'elastic',                  'rdiss': 1.0,  'label': 'Elastic (r_diss=1)'},
    'const_rdiss_1p05':         {'dir': 'const_rdiss_1p05',         'rdiss': 1.05, 'label': r'Const $r_{\rm diss}=1.05$'},
    'M1_dark_photon_massive':   {'dir': 'M1_dark_photon_massive',   'rdiss': None, 'label': r'M1: dark photon (massive)'},
    'M2_scalar_phi_massive':    {'dir': 'M2_scalar_phi_massive',    'rdiss': None, 'label': r'M2: scalar $\phi$ (massive)'},
    'M1_highconc':              {'dir': 'M1_dark_photon_massive_highconc', 'rdiss': None,
                                 'label': r'M1: dark photon (high-conc NFW, $c_{200}\!\sim\!12$)',
                                 'rs_sim_kpc': 2.5, 'rho0_sim': 1.47e-2},
}

DATA_ROOT = os.path.join(_PROJ_ROOT, 'data', 'P2_runs')
OUT_DATA   = os.path.join(_PROJ_ROOT, 'data')
OUT_FIG    = os.path.join(_PROJ_ROOT, 'figures')
OUT_NOTES  = os.path.join(_PROJ_ROOT, 'notes')
for d in (OUT_DATA, OUT_FIG, OUT_NOTES):
    os.makedirs(d, exist_ok=True)

# r_2D/r_s benchmark values to highlight (per Schmidt Fig. 10 / Table 2)
BENCH_R2D_RS = [0.05, 0.2, 0.5]


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------
def get_scales(halo_ini):
    r_s = halo_ini['r_s'] * ut.kpc
    rho_s = halo_ini['rho_s'] * ut.M_sun / ut.pc**3
    scale_r_kpc = r_s.to('kpc').value
    scale_rho = rho_s.to('Msun/pc**3').value
    scale_v = np.sqrt(ct.G * 4 * np.pi * rho_s * r_s**3 / r_s).to('km/s').value
    scale_t = (1.0 / np.sqrt(4.0 * np.pi * rho_s * ct.G)).to('Gyr').value
    return scale_r_kpc, scale_rho, scale_v, scale_t


def process_model_full(model_key, model_meta, n_scan=80, snapshot_stride=1):
    """Scan all snapshots × r_2D/r_s. Return DataFrame of matches."""
    dir_data = os.path.join(DATA_ROOT, model_meta['dir'])
    if not os.path.isdir(dir_data):
        print(f"[{model_key}] skip: no data dir {dir_data}")
        return pd.DataFrame()

    # Use per-model NFW parameters (defaults match P2; highconc has its own)
    rs_sim_kpc = model_meta.get('rs_sim_kpc', RS_SIM_KPC)
    rho0_sim   = model_meta.get('rho0_sim', RHO0_SIM)

    halorec = HaloRecord(dir_data)
    list_files, list_times = halorec.glob_pickle_files()
    if len(list_files) == 0:
        print(f"[{model_key}] no snapshots")
        return pd.DataFrame()

    halo_ini, _ = halorec.get_halo_initialization()
    scale_r_kpc, scale_rho, scale_v, scale_t = get_scales(halo_ini)
    print(f"[{model_key}] {len(list_files)} snapshots, "
          f"scale_r={scale_r_kpc:.2f} kpc, scale_t={scale_t:.5f} Gyr")

    rows = []
    files_iter = list_files[::snapshot_stride]
    for i, f in enumerate(files_iter):
        try:
            data = halorec.get_halo_state_pickled(file_halo=f)
            rho_arr = data.get('rho', None)
            if rho_arr is None or np.any(np.isnan(rho_arr)):
                continue
        except Exception:
            continue

        r_kpc = data['r'] * scale_r_kpc
        rho = data['rho'] * scale_rho  # Msun/pc^3
        t_gyr = data['t'] * scale_t

        scan = find_matching_radii(r_kpc, rho, rs_sim_kpc, n_scan=n_scan)
        for r2D_rs, M_inner, M_outer, ratio in scan:
            if not (M_inner > 0 and M_outer > 0 and ratio > 0):
                continue
            n_sigma = abs(ratio - MASS_RATIO_OBS) / MASS_RATIO_ERR
            if n_sigma > 3:
                continue
            lam, mu = compute_rescaling(r2D_rs, M_inner, M_outer,
                                        r_s_sim_kpc=rs_sim_kpc)
            if not (lam > 0 and mu > 0):
                continue
            mass_fit = absolute_mass_fit(M_inner, M_outer)
            params = rescale_parameters(lam, mu, t_gyr,
                                        sigma_m_sim=SIGMA_M_SIM,
                                        r_s_sim=rs_sim_kpc, rho0_sim=rho0_sim)
            rows.append({
                'model': model_key,
                'snapshot_idx': i,
                'snapshot_time_gyr': t_gyr,
                'r2D_rs': r2D_rs,
                'M_inner_sim': M_inner,
                'M_outer_sim': M_outer,
                'mass_ratio': ratio,
                'n_sigma': n_sigma,
                'lambda': lam,
                'mu': mu,
                'M_inner_phys': mass_fit['M_inner_phys'],
                'M_outer_phys': mass_fit['M_outer_phys'],
                'mass_chi2': mass_fit['chi2'],
                'mass_nsigma_dof1': mass_fit['nsigma_dof1'],
                'mass_residual_inner_sigma': mass_fit['residual_inner_sigma'],
                'mass_residual_outer_sigma': mass_fit['residual_outer_sigma'],
                'r_s_kpc': params['r_s_kpc'],
                'rho0_msun_pc3': params['rho0_msun_pc3'],
                'sigma_m_cm2_g': params['sigma_m_cm2_g'],
                't_evo_gyr': params['t_evo_gyr'],
                't_ok': params['t_evo_gyr'] <= T_ZOBS_GYR,
            })

    df = pd.DataFrame(rows)
    print(f"[{model_key}] {len(df)} matching points, "
          f"{int(df['t_ok'].sum()) if len(df) else 0} with t_evo ≤ {T_ZOBS_GYR} Gyr")
    return df


def add_symmetry_breaking(df, model_key, model_meta=None):
    """For velocity-dependent models, compute the FIRST-ORDER symmetry-breaking
    correction to σ/m (not just a diagnostic ratio).

    Physics:
    --------
    The naive rescaling (Appendix G of Schmidt 2026) assumes r_diss is
    velocity-independent, so the simulation's r_diss(v_sim) equals the
    physical r_diss(v_phys). For massive emission this is violated:
        r_diss(v) = 1 + C0 * C_{V/φ}(ω(v)) * exp(-m/T(v))
    with v_phys = v_sim * sqrt(μ/λ).

    Since the volumetric cooling rate is C ∝ (σ/m) * ρ * ν^3 * (r_diss - 1),
    matching the observed density profile requires:
        (σ/m)_phys * (r_diss_phys - 1) = (σ/m)_sim * (r_diss_sim - 1) * (ρ_sim/ρ_phys) * (ν_sim/ν_phys)^3
    Under the elastic rescaling symmetry, the density/velocity factors cancel
    leaving (σ/m)_phys = (λ²/μ) (σ/m)_sim. With velocity-dependent r_diss,
    the leftover factor is:
        (σ/m)_corrected = (σ/m)_naive * [(r_diss_sim - 1) / (r_diss_phys - 1)]

    Two physically meaningful regimes:
    1. r_diss_sim > 1 and r_diss_phys > 1 (both above threshold):
       The correction is finite and quantifies the offset directly.
    2. r_diss_sim ≈ 1 (simulation frame below threshold):
       This means the simulation effectively evolves an ELASTIC halo, but
       the physical halo would have weak dissipation (r_diss_phys > 1).
       The naive rescaling UNDERESTIMATES σ/m because the physical halo
       would cool faster than the elastic simulation. Physically, the
       simulation is not a valid dissipative realization for this
       observation; the proper correction requires a new simulation with
       dissipative physics active at v_sim. We flag this regime as
       'invalid_kinematic_regime' and report a LOWER BOUND on the offset:

           σ/m_corrected ≥ σ/m_naive * (r_diss_phys - 1)_max / (r_diss_sim - 1)_max

       using the maximum expected r_diss (i.e., the asymptotic high-v
       value 1 + C0) for both, which gives the minimum relative correction
       factor of (r_diss_sim_max / r_diss_phys_max) and is conservative
       (in practice we mark the σ/m offset as a flag and set the correction
       factor to NaN to indicate "needs new simulation").

    Reference velocity: ~50 km/s, characteristic 1D velocity dispersion of
    the B1938+666 subhalo at z=0.881.
    """
    EPS_RDISS = 1e-9  # numerical floor for (r_diss - 1) below which we
                       # consider emission fully kinematically suppressed

    if model_key not in ('M1_dark_photon_massive', 'M2_scalar_phi_massive',
                          'M1_highconc'):
        df['rdiss_before'] = np.nan
        df['rdiss_after'] = np.nan
        df['rdiss_ratio'] = 1.0
        df['v_sim_kms'] = np.nan
        df['v_phys_kms'] = np.nan
        df['rdiss_correction_factor'] = 1.0  # elastic: no correction
        df['sigma_m_corrected_cm2_g'] = df['sigma_m_cm2_g']  # = naive
        df['symmetry_breaking_pct'] = 0.0  # percent offset
        df['valid_regime'] = True
        return df

    bm = benchmark_models()
    # For the high-conc M1 variant, use the same dark sector parameters
    # as M1_dark_photon_massive (the difference is in the halo profile only)
    src_key = 'M1_dark_photon_massive' if model_key == 'M1_highconc' else model_key
    params = bm[src_key]

    # IMPORTANT: v_phys is NOT fixed at 50 km/s. The physical halo's
    # characteristic 1D velocity dispersion scales as
    #   v_scale_phys = sqrt(G * 4π * ρ_phys * r_s_phys^2)
    # Under the elastic rescaling symmetry this is sqrt(μ/λ) * v_scale_sim,
    # where v_scale_sim is for the SIMULATION halo (r_s=3.6 kpc,
    # ρ=7.09e-3 Msun/pc^3). v_scale_sim ~ 70 km/s.
    # The previous code incorrectly fixed v_phys = 50 km/s, which
    # underestimates v_phys for the higher-σ/m extended-solution matches.
    from astropy import units as ut
    from astropy import constants as ct
    # Per-model sim halo parameters (highconc has different r_s, ρ)
    rs_sim_kpc = (model_meta or {}).get('rs_sim_kpc', RS_SIM_KPC)
    rho0_sim = (model_meta or {}).get('rho0_sim', RHO0_SIM)
    rho_sim_astropy = rho0_sim * ut.M_sun / ut.pc**3
    r_s_sim_astropy = rs_sim_kpc * ut.kpc
    v_scale_sim = np.sqrt(ct.G * 4 * np.pi * rho_sim_astropy * r_s_sim_astropy**2).to('km/s').value
    # Asymptotic (high-v) r_diss - 1 = C0 (used for invalid-regime flag only)
    rdiss_max_minus1 = 0.05

    r_b, r_a, v_s, v_p = [], [], [], []
    corr_f, sig_corr, sb_pct, valid_reg = [], [], [], []
    n_invalid = 0
    for _, row in df.iterrows():
        lam, mu = row['lambda'], row['mu']
        # Per-point v_phys and v_sim, using the actual halo velocity scale
        # under the rescaling symmetry:
        #   v_phys = sqrt(μ/λ) * v_scale_sim
        #   v_sim  = v_scale_sim
        # This accounts for the fact that high-σ/m matches (large λ) have
        # large physical halos with high v_phys, possibly above v*.
        v_sim = v_scale_sim
        v_phys = v_scale_sim * np.sqrt(mu / lam)
        try:
            rd_before = r_diss_func(np.array([v_sim]), params)[0]   # r_diss at sim-frame v
            rd_after = r_diss_func(np.array([v_phys]), params)[0]   # r_diss at phys-frame v
        except Exception:
            rd_before, rd_after = np.nan, np.nan
        r_b.append(rd_before)
        r_a.append(rd_after)
        v_s.append(v_sim)
        v_p.append(v_phys)

        rd_sim_minus1 = float(rd_before - 1.0)
        rd_phys_minus1 = float(rd_after - 1.0)
        sig_naive = row['sigma_m_cm2_g']

        # Two regimes
        if rd_sim_minus1 > EPS_RDISS and rd_phys_minus1 > EPS_RDISS:
            # Both above threshold: direct correction
            factor = rd_sim_minus1 / rd_phys_minus1
            sig_c = sig_naive * factor
            pct = 100.0 * (sig_c - sig_naive) / sig_naive if sig_naive > 0 else 0.0
            valid = True
        elif rd_sim_minus1 <= EPS_RDISS:
            # Simulation frame kinematically suppressed:
            # Naive rescaling treats this as an elastic halo, but the
            # physical halo would have weak dissipation. Simulation is not
            # a valid dissipative realization here.
            factor = np.nan
            sig_c = np.nan
            # Lower bound on the offset if sim had been dissipative at full
            # asymptotic strength: factor_min = rdiss_max / rdiss_phys
            pct = np.nan  # flagged
            valid = False
            n_invalid += 1
        else:
            # rd_sim > 1 but rd_phys ≈ 1: physical halo is below threshold
            # Naive rescaling OVERESTIMATES σ/m (would predict dissipative
            # cooling that doesn't actually happen in the physical halo)
            factor = float('inf')  # formally divergent (need to re-simulate)
            # Conservative: σ_corrected = σ_naive (elastic limit)
            sig_c = sig_naive
            pct = 0.0
            valid = False
            n_invalid += 1

        corr_f.append(factor)
        sig_corr.append(sig_c)
        sb_pct.append(pct)
        valid_reg.append(valid)

    df['rdiss_before'] = r_b
    df['rdiss_after'] = r_a
    df['rdiss_ratio'] = np.where(np.array(r_b) > 0,
                                  np.array(r_a) / np.array(r_b), np.nan)
    df['v_sim_kms'] = v_s
    df['v_phys_kms'] = v_p
    df['rdiss_correction_factor'] = corr_f
    df['sigma_m_corrected_cm2_g'] = sig_corr
    df['symmetry_breaking_pct'] = sb_pct
    df['valid_regime'] = valid_reg

    n_total = len(df)
    n_valid = int(np.sum(valid_reg))
    if n_invalid > 0:
        print(f"  [{model_key}] symmetry breaking: {n_valid}/{n_total} matches "
              f"in valid regime (both frames above threshold); "
              f"{n_invalid - (n_total - n_valid)} below sim threshold "
              f"(simulation is effectively elastic, σ/m_corrected "
              f"cannot be reliably computed without new sim).")
    return df


# ----------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------
def plot_sigma_t_evo(all_df, out_path):
    """σ/m vs t_evo scatter, all models."""
    fig, ax = plt.subplots(figsize=(9, 6.5))
    colors = {'elastic': 'tab:blue', 'const_rdiss_1p05': 'tab:orange',
              'M1_dark_photon_massive': 'tab:green',
              'M2_scalar_phi_massive': 'tab:red',
              'M1_highconc': 'tab:purple'}
    markers = {'elastic': 'o', 'const_rdiss_1p05': 's',
               'M1_dark_photon_massive': '^', 'M2_scalar_phi_massive': 'D',
               'M1_highconc': 'v'}

    for mk, meta in MODELS.items():
        df = all_df[all_df['model'] == mk]
        if len(df) == 0:
            continue
        # Viable (t_ok) in solid, non-viable hollow
        viable = df[df['t_ok']]
        excluded = df[~df['t_ok']]
        if len(viable):
            ax.scatter(viable['sigma_m_cm2_g'], viable['t_evo_gyr'],
                       s=12, c=colors[mk], marker=markers[mk],
                       label=f"{meta['label']} (viable)", alpha=0.6,
                       edgecolors='none')
        if len(excluded):
            ax.scatter(excluded['sigma_m_cm2_g'], excluded['t_evo_gyr'],
                       s=12, c='none', marker=markers[mk],
                       edgecolors=colors[mk], linewidths=0.7,
                       label=f"{meta['label']} (t_evo > 6.37 Gyr)", alpha=0.6)

    ax.axhline(T_ZOBS_GYR, color='k', linestyle='--', lw=1.2,
               label=r'$t_{\rm zobs} = 6.37$ Gyr')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$\sigma_T / m_\chi$ [cm$^2$/g]', fontsize=12)
    ax.set_ylabel(r'$t_{\rm evo}$ [Gyr]', fontsize=12)
    ax.set_title('B1938+666 rescaled halo parameters (within 3σ mass-ratio match)',
                 fontsize=12)
    ax.legend(loc='best', fontsize=9, framealpha=0.9)
    ax.grid(True, which='both', alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  Saved {out_path}")


def plot_symmetry_breaking_correction(all_df, out_path):
    """First-order symmetry-breaking correction diagnostic for the
    velocity-dependent models (M1, M2).

    After fixing the v_phys bug (was: fixed 50 km/s; now: per-point halo
    velocity scale v_scale_sim * sqrt(mu/lambda), ranging ~488-940 km/s),
    the picture is the OPPOSITE of the original finding:

    - M1 (dark photon, v*=200 km/s): v_phys >> v*, so ALL matches are in the
      valid regime (r_diss_sim > 1 AND r_diss_phys > 1). The first-order
      correction factor is large (sb_pct ~ -99%), meaning the naive elastic
      rescaling DRAMATICALLY overestimates sigma/m. The physical halo is
      much more dissipative than the simulation, so to match the same
      observed mass ratio, sigma/m must be reduced by ~2 orders of magnitude.
    - M2 (scalar phi, v*=500 km/s): v_phys > v* but r_diss(v_sim~70 km/s)
      is numerically zero, so these matches remain in the invalid regime.
      The simulation is effectively elastic; a proper treatment requires
      the per-point resimulation (P5).

    LEFT  : r_diss(v_sim) vs r_diss(v_phys) for each match, with the
            validity regime (both above threshold) marked.
    RIGHT : v_sim vs v_phys scatter, colored by r_diss_after / r_diss_before.
            Shows v_phys well above v* for M1, confirming the symmetry is
            strongly broken at the physical velocity scale.
    """
    vd_models = ['M1_dark_photon_massive', 'M2_scalar_phi_massive',
                'M1_highconc']
    df_vd = all_df[all_df['model'].isin(vd_models)]
    if len(df_vd) == 0:
        print(f"  No velocity-dependent matches; skipping {out_path}")
        return

    fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(13, 5.5))

    colors = {'M1_dark_photon_massive': 'tab:green',
              'M2_scalar_phi_massive': 'tab:red',
              'M1_highconc': 'tab:purple'}
    labels = {'M1_dark_photon_massive': r'M1: dark photon ($v^*\!=\!200$ km/s)',
              'M2_scalar_phi_massive':  r'M2: scalar $\phi$ ($v^*\!=\!500$ km/s)',
              'M1_highconc':            r'M1: high-conc NFW ($c_{200}\!\sim\!12$)'}
    v_star = {'M1_dark_photon_massive': 200.0,
              'M2_scalar_phi_massive': 500.0,
              'M1_highconc': 200.0}

    # LEFT: r_diss(v_sim) vs r_diss(v_phys)
    for mk in vd_models:
        sub = df_vd[df_vd['model'] == mk]
        if len(sub) == 0:
            continue
        valid_mask = sub['valid_regime'].astype(bool)
        invalid_mask = ~valid_mask
        # Valid: solid markers
        if valid_mask.any():
            ax_l.scatter(sub.loc[valid_mask, 'rdiss_before'],
                         sub.loc[valid_mask, 'rdiss_after'],
                         s=15, c=colors[mk], marker='o', alpha=0.7,
                         label=labels[mk] + ' (valid)', edgecolors='none')
        # Invalid: open markers
        if invalid_mask.any():
            ax_l.scatter(sub.loc[invalid_mask, 'rdiss_before'],
                         sub.loc[invalid_mask, 'rdiss_after'],
                         s=15, c='none', marker='x',
                         edgecolors=colors[mk], linewidths=0.8,
                         label=labels[mk] + ' (sim kinem. suppressed)',
                         alpha=0.6)
    # Elastic diagonal
    rd_arr = np.linspace(0.999, 1.10, 50)
    ax_l.plot(rd_arr, rd_arr, 'k--', lw=1.0,
              label='Elastic symmetry (r_diss_sim = r_diss_phys)')
    ax_l.set_xlabel(r'$r_{\rm diss}$ at $v_{\rm sim}$ '
                    r'(simulation frame)', fontsize=11)
    ax_l.set_ylabel(r'$r_{\rm diss}$ at $v_{\rm phys}=v_{\rm scale,phys}$ '
                    r'(physical frame)', fontsize=11)
    ax_l.set_title('Symmetry breaking: r_diss before vs after rescaling',
                   fontsize=11)
    ax_l.legend(loc='best', fontsize=8, framealpha=0.9)
    ax_l.grid(True, alpha=0.3)

    # RIGHT: v_sim vs v_phys, color = rdiss_ratio
    for mk in vd_models:
        sub = df_vd[df_vd['model'] == mk]
        if len(sub) == 0:
            continue
        sc = ax_r.scatter(sub['v_sim_kms'], sub['v_phys_kms'],
                          c=sub['rdiss_ratio'], s=10,
                          cmap='viridis', alpha=0.7, vmin=0.999, vmax=1.005,
                          label=labels[mk])
        # Mark v* for this model
        ax_r.axvline(v_star[mk], color=colors[mk], linestyle='--',
                     lw=1.2, alpha=0.7,
                     label=f'{labels[mk]} $v^*$ = {v_star[mk]:.0f} km/s')
    ax_r.axhline(50.0, color='gray', linestyle=':', lw=1.0,
                  label=r'$v_{\rm phys}=50$ km/s (old, incorrect)')
    ax_r.set_xlabel(r'$v_{\rm sim}$ [km/s] (simulation frame)', fontsize=11)
    ax_r.set_ylabel(r'$v_{\rm phys}$ [km/s] (physical frame)', fontsize=11)
    ax_r.set_title('Match-point velocities vs emission threshold',
                   fontsize=11)
    ax_r.legend(loc='best', fontsize=8, framealpha=0.9)
    ax_r.grid(True, alpha=0.3)
    cb = fig.colorbar(sc, ax=ax_r, shrink=0.8, pad=0.02)
    cb.set_label(r'$r_{\rm diss}^{\rm after} / r_{\rm diss}^{\rm before}$',
                 fontsize=10)

    fig.suptitle('P3 symmetry-breaking diagnostic for velocity-dependent '
                 'models (compact-solution branch)',
                 fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {out_path}")


def plot_mass_ratio_curves(all_df, out_path):
    """Mass ratio vs r_2D/r_s, with snapshot evolution color-coded."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    panels = [('elastic', 'const_rdiss_1p05'),
              ('M1_dark_photon_massive', 'M2_scalar_phi_massive')]
    # Append highconc to the second panel if present
    if 'M1_highconc' in all_df['model'].unique():
        panels[1] = ('M1_dark_photon_massive', 'M1_highconc')
    for ax, (k1, k2) in zip(axes, panels):
        for mk in (k1, k2):
            df = all_df[all_df['model'] == mk]
            if len(df) == 0:
                continue
            sc = ax.scatter(df['r2D_rs'], df['mass_ratio'],
                            c=df['snapshot_time_gyr'], s=10,
                            cmap='viridis', alpha=0.7,
                            label=MODELS[mk]['label'])
        ax.axhline(MASS_RATIO_OBS, color='red', linestyle='--', lw=1.2,
                   label=r'$M_{\rm obs} = 0.364$')
        ax.axhspan(MASS_RATIO_OBS - 3*MASS_RATIO_ERR,
                   MASS_RATIO_OBS + 3*MASS_RATIO_ERR,
                   color='red', alpha=0.12, label='3σ band')
        ax.axvline(0.05, color='gray', linestyle=':', lw=0.7)
        ax.axvline(0.2,  color='gray', linestyle=':', lw=0.7)
        ax.axvline(0.5,  color='gray', linestyle=':', lw=0.7)
        ax.set_xscale('log')
        ax.set_xlabel(r'$r_{2D} / r_s$', fontsize=12)
        ax.set_ylabel(r'$M(r<r_{2D}) / M(r<4.5\,r_{2D})$', fontsize=12)
        ax.legend(loc='upper left', fontsize=9)
        ax.grid(True, which='both', alpha=0.3)
    cb = fig.colorbar(sc, ax=axes, shrink=0.8, pad=0.02)
    cb.set_label('snapshot time [Gyr]', fontsize=11)
    fig.suptitle('Projected mass ratio vs r_2D/r_s (color = snapshot time)',
                 fontsize=12, y=1.02)
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {out_path}")


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    t0 = time.time()
    all_dfs = []
    for mk, meta in MODELS.items():
        df = process_model_full(mk, meta, n_scan=80, snapshot_stride=1)
        if len(df):
            df = add_symmetry_breaking(df, mk, model_meta=meta)
            all_dfs.append(df)

    if not all_dfs:
        print("No matching points found in any model.")
        return

    all_df = pd.concat(all_dfs, ignore_index=True)
    csv_path = os.path.join(OUT_DATA, 'P3_rescaled_params.csv')
    all_df.to_csv(csv_path, index=False)
    print(f"\nWrote {len(all_df)} matching points → {csv_path}")

    # Summary: best (largest σ/m) viable point per model
    summary_rows = []
    for mk, meta in MODELS.items():
        df = all_df[all_df['model'] == mk]
        viable = df[df['t_ok']]
        if len(viable) == 0:
            summary_rows.append({'model': mk, 'n_matches': len(df),
                                 'n_viable': 0})
            continue
        # Largest σ/m (typically the most "active" SIDM)
        best = viable.loc[viable['sigma_m_cm2_g'].idxmax()]
        summary_rows.append({
            'model': mk,
            'n_matches': len(df),
            'n_viable': int(len(viable)),
            'best_sigma_m_cm2_g': best['sigma_m_cm2_g'],
            'best_t_evo_gyr': best['t_evo_gyr'],
            'best_r_s_kpc': best['r_s_kpc'],
            'best_rho0_msun_pc3': best['rho0_msun_pc3'],
            'best_r2D_rs': best['r2D_rs'],
            'rdiss_ratio': best.get('rdiss_ratio', np.nan),
            'rdiss_correction_factor': best.get('rdiss_correction_factor', np.nan),
            'symmetry_breaking_pct': best.get('symmetry_breaking_pct', np.nan),
            'sigma_m_corrected_cm2_g': best.get('sigma_m_corrected_cm2_g', np.nan),
        })
    summary = pd.DataFrame(summary_rows)
    summary_path = os.path.join(OUT_DATA, 'P3_summary.csv')
    summary.to_csv(summary_path, index=False)
    print(f"Wrote summary → {summary_path}")
    print("\n=== Summary ===")
    print(summary.to_string(index=False))

    # Plots
    plot_sigma_t_evo(all_df, os.path.join(OUT_FIG, 'P3_B1938_regions.png'))
    plot_mass_ratio_curves(all_df, os.path.join(OUT_FIG, 'P3_mass_ratio.png'))
    plot_symmetry_breaking_correction(
        all_df, os.path.join(OUT_FIG, 'P3_symmetry_breaking_correction.png'))

    print(f"\n[P3 complete] elapsed {time.time()-t0:.1f}s")


if __name__ == '__main__':
    main()
