"""
B1938+666 rescaling for dissipative SIDM models.

Implements the rescaling symmetry of Schmidt et al. 2026 Appendix G:
  r → λ r,  m → μ m,  v → sqrt(μ/λ) v,  t → sqrt(λ³/μ) t

For each snapshot of a P2 evolution, computes:
  1. Projected 2D enclosed mass M(r<r_2D) and M(r<4.5*r_2D)
  2. The mass ratio M(r_2D)/M(4.5*r_2D) and compares to observed 0.364±0.022
  3. For matching snapshots, the rescaling parameters λ (Eq G.9) and μ (Eq G.11)
  4. Rescaled physical parameters: r_s, ρ_0, σ_T/m, t_evo

Observational constraints (Vegetti et al. 2026):
  M(r<20 pc) = (4.25 ± 0.21) × 10⁵ M☉
  M(r<90 pc) = (1.167 ± 0.039) × 10⁶ M☉
  z_obs = 0.881 → t(z_obs) = 6.37 Gyr → t_evo ≤ 6.37 Gyr

Key difference from Schmidt et al. 2026:
  For velocity-dependent models (M1, M2), the rescaling symmetry is BROKEN
  because r_diss(v) introduces a fixed velocity scale v* = sqrt(2m/μ).
  The rescaling v → sqrt(μ/λ) v changes v/v*, hence r_diss changes.
  We quantify this effect by comparing r_diss before and after rescaling.
"""

from __future__ import annotations

import os, sys
import numpy as np
from scipy.interpolate import interp1d

# Compatibility shim: np.trapz was removed in numpy 2.x; renamed to np.trapezoid.
if not hasattr(np, 'trapz'):
    np.trapz = np.trapezoid

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJ_ROOT = os.path.normpath(os.path.join(_HERE, '..', '..'))
sys.path.insert(0, os.path.join(_PROJ_ROOT, 'external', 'gravothermalsidm'))

from astropy import units as ut
from astropy import constants as ct
from SourcePy.record import HaloRecord


# ----------------------------------------------------------------------
# Observational constraints (Vegetti et al. 2026)
# ----------------------------------------------------------------------
M_20PC = 4.25e5      # M_sun
M_20PC_ERR = 0.21e5  # M_sun
M_90PC = 1.167e6     # M_sun
M_90PC_ERR = 0.039e6 # M_sun
R_INNER_PC = 20.0   # pc
R_OUTER_PC = 90.0   # pc
R_RATIO = R_OUTER_PC / R_INNER_PC  # = 4.5
MASS_RATIO_OBS = M_20PC / M_90PC  # = 0.364
MASS_RATIO_ERR = 0.022  # 1σ
T_ZOBS_GYR = 6.37   # Gyr (cosmic age at z=0.881)


def absolute_mass_fit(M_inner_sim, M_outer_sim):
    """Fit one mass-rescaling factor to both absolute B1938 masses.

    The rescaling maps simulated enclosed masses to ``M_phys = mu * M_sim``.
    This weighted least-squares fit is independent of the ratio preselection,
    which cannot constrain the absolute mass normalization by itself.
    """
    masses = np.asarray([M_inner_sim, M_outer_sim], dtype=float)
    observed = np.asarray([M_20PC, M_90PC], dtype=float)
    errors = np.asarray([M_20PC_ERR, M_90PC_ERR], dtype=float)
    if not np.all(np.isfinite(masses)) or np.any(masses <= 0):
        raise ValueError("simulated enclosed masses must be finite and positive")
    weights = 1.0 / errors**2
    denominator = float(np.sum(weights * masses**2))
    mu = float(np.sum(weights * masses * observed) / denominator)
    predicted = mu * masses
    residual = (predicted - observed) / errors
    return {
        'mu': mu,
        'M_inner_phys': float(predicted[0]),
        'M_outer_phys': float(predicted[1]),
        'chi2': float(np.sum(residual**2)),
        'dof': 1,
        'nsigma_dof1': float(np.sqrt(np.sum(residual**2))),
        'residual_inner_sigma': float(residual[0]),
        'residual_outer_sigma': float(residual[1]),
    }

# Simulation parameters (from P2)
RS_SIM_KPC = 3.6
RHO0_SIM = 7.09e-3  # Msun/pc^3
SIGMA_M_SIM = 50.0  # cm^2/g


# ----------------------------------------------------------------------
# Projected 2D enclosed mass
# ----------------------------------------------------------------------
def projected_enclosed_mass(r_3d, rho_3d, r_2d_target, r_unit='kpc',
                             rho_unit='Msun_pc3'):
    """Compute the 2D projected enclosed mass M_2D(< r_2d_target) from a 3D
    density profile rho(r).

    The fluid state stores a constant average density in each spherical shell.
    For a sphere of radius ``r``, the volume inside a projected cylinder of
    radius ``R`` is

        V_cyl(r, R) = 4*pi/3 * [r^3 - max(r^2 - R^2, 0)^(3/2)].

    Subtracting this volume at adjacent shell edges gives the exact projected
    volume of every piecewise-constant shell.  This avoids the endpoint
    singularity and resolution bias of a sampled Abel integral.

    Args:
        r_3d : array of 3D radii (same units as r_2d_target)
        rho_3d : array of 3D density (same units as mass/length^3 that will
                 give M in mass units when integrated against r*dr)
        r_2d_target : projected radius in the SAME units as r_3d
        r_unit, rho_unit : only used to verify the (r,rho) pair yields a mass.
            If r is in kpc but rho is in M_sun/pc^3, we internally convert r
            to pc so that the integral returns M_sun.

    Returns:
        M_2d in M_sun (assuming rho is in M_sun/pc^3 and r is in kpc, or any
        consistent mass/length pair).
    """
    # The fluid profile stores a shell-average density at each shell's outer
    # radius.  Project each constant-density shell analytically instead of
    # applying a singular Abel integral to sparse shell samples.
    idx = np.argsort(r_3d)
    r = np.asarray(r_3d, dtype=float)[idx]
    rho = np.asarray(rho_3d, dtype=float)[idx]
    if r.ndim != 1 or rho.ndim != 1 or len(r) != len(rho) or len(r) == 0:
        raise ValueError("r_3d and rho_3d must be non-empty 1D arrays of equal length")
    if np.any(~np.isfinite(r)) or np.any(~np.isfinite(rho)):
        raise ValueError("r_3d and rho_3d must be finite")
    if np.any(r <= 0) or np.any(np.diff(r) <= 0) or np.any(rho < 0):
        raise ValueError("radii must increase and densities must be non-negative")
    if r_2d_target <= 0:
        return 0.0

    # If r is in kpc but rho is in M_sun/pc^3, convert r to pc so the integral
    # has units compatible with shell volumes in pc^3.
    if r_unit == 'kpc' and rho_unit == 'Msun_pc3':
        r_pc = r * 1000.0
        r_2d_pc = r_2d_target * 1000.0
    else:
        r_pc = r
        r_2d_pc = r_2d_target

    r_inner = np.concatenate(([0.0], r_pc[:-1]))

    def cylinder_volume(radius):
        outside = np.clip(radius**2 - r_2d_pc**2, 0.0, None) ** 1.5
        return (4.0 * np.pi / 3.0) * (radius**3 - outside)

    projected_shell_volume = cylinder_volume(r_pc) - cylinder_volume(r_inner)
    return float(np.sum(rho * projected_shell_volume))


# ----------------------------------------------------------------------
# Find matching r_2D / r_s for a given snapshot
# ----------------------------------------------------------------------
def find_matching_radii(r_3d_kpc, rho_3d, r_s_kpc, n_scan=200):
    """Scan r_2D/r_s to find where M(r_2D)/M(4.5*r_2D) matches the observed ratio.

    Returns array of (r_2D/r_s, M_inner, M_outer, ratio, chi2) for all scanned radii.
    """
    r_2D_over_rs_arr = np.logspace(-2, 0.5, n_scan)  # 0.01 to ~3

    results = []
    for r2D_rs in r_2D_over_rs_arr:
        r_inner_kpc = r2D_rs * r_s_kpc
        r_outer_kpc = R_RATIO * r_inner_kpc

        # Compute projected enclosed masses
        M_inner = projected_enclosed_mass(r_3d_kpc, rho_3d, r_inner_kpc)
        M_outer = projected_enclosed_mass(r_3d_kpc, rho_3d, r_outer_kpc)

        if M_outer > 0:
            ratio = M_inner / M_outer
        else:
            ratio = np.nan

        results.append((r2D_rs, M_inner, M_outer, ratio))

    return np.array(results)


# ----------------------------------------------------------------------
# Compute rescaling parameters λ and μ (Appendix G)
# ----------------------------------------------------------------------
def compute_rescaling(r2D_rs, M_inner_sim, M_outer_sim,
                      r_s_sim_kpc=RS_SIM_KPC):
    """Compute λ and μ from Appendix G Eqs (G.9) and (G.11).

    λ = (20 pc / r_s_sim) × (r_2D/r_s)^{-1} = (20/3600) × (r_2D/r_s)^{-1}

    μ = [M_sim_inner * M_obs_inner / σ_inner² + M_sim_outer * M_obs_outer / σ_outer²]
        / [M_sim_inner² / σ_inner² + M_sim_outer² / σ_outer²]
    """
    # λ from Eq (G.9)
    r_inner_sim_kpc = r2D_rs * r_s_sim_kpc  # in kpc
    r_inner_sim_pc = r_inner_sim_kpc * 1000  # in pc
    lam = R_INNER_PC / r_inner_sim_pc

    # μ from Eq (G.11)
    mu = absolute_mass_fit(M_inner_sim, M_outer_sim)['mu']

    return lam, mu


# ----------------------------------------------------------------------
# Apply rescaling to get physical parameters
# ----------------------------------------------------------------------
def rescale_parameters(lam, mu, t_sim_gyr, sigma_m_sim=SIGMA_M_SIM,
                       r_s_sim=RS_SIM_KPC, rho0_sim=RHO0_SIM):
    """Compute rescaled physical parameters.

    From Appendix G:
      r_s_phys = λ × r_s_sim
      ρ_0_phys = μ × ρ_0_sim / λ³
      σ_T/m_phys = λ²/μ × σ_T/m_sim
      t_evo_phys = sqrt(λ³/μ) × t_sim
    """
    r_s_phys = lam * r_s_sim  # kpc
    rho0_phys = mu * rho0_sim / lam**3  # Msun/pc^3
    sigma_m_phys = lam**2 / mu * sigma_m_sim  # cm^2/g
    t_evo_phys = np.sqrt(lam**3 / mu) * t_sim_gyr  # Gyr

    return {
        'r_s_kpc': r_s_phys,
        'rho0_msun_pc3': rho0_phys,
        'sigma_m_cm2_g': sigma_m_phys,
        't_evo_gyr': t_evo_phys,
    }


# ----------------------------------------------------------------------
# Process all snapshots of a model
# ----------------------------------------------------------------------
def process_model(model_name, dir_data, scale_r_kpc, scale_rho,
                  scale_t_gyr, scale_v_kms, sigma_m_sim=SIGMA_M_SIM):
    """Process all snapshots of a model, find matches, compute rescaled params.

    Returns list of dicts with matching results.
    """
    halorec = HaloRecord(dir_data)
    list_files, list_times = halorec.glob_pickle_files()
    if len(list_files) == 0:
        print(f"  [{model_name}] No snapshots found")
        return []

    print(f"  [{model_name}] Processing {len(list_files)} snapshots...")

    matches = []
    n_within_3sigma = 0

    for i, f in enumerate(list_files):
        try:
            data = halorec.get_halo_state_pickled(file_halo=f)
            if not data or np.any(np.isnan(data.get('rho', [np.nan]))):
                continue
        except:
            continue

        # Dimensionful quantities
        r_kpc = data['r'] * scale_r_kpc
        rho = data['rho'] * scale_rho  # Msun/pc^3
        t_gyr = data['t'] * scale_t_gyr

        # Scan r_2D/r_s
        scan = find_matching_radii(r_kpc, rho, RS_SIM_KPC, n_scan=100)

        # Find where mass ratio matches observation
        for row in scan:
            r2D_rs, M_inner, M_outer, ratio = row
            if np.isnan(ratio) or M_inner <= 0 or M_outer <= 0:
                continue

            # Check if within 3σ of observed ratio
            n_sigma = abs(ratio - MASS_RATIO_OBS) / MASS_RATIO_ERR
            if n_sigma > 3:
                continue

            n_within_3sigma += 1

            # Compute rescaling and retain the joint absolute-mass fit.
            lam, mu = compute_rescaling(r2D_rs, M_inner, M_outer)
            if np.isnan(mu) or mu <= 0 or lam <= 0:
                continue
            mass_fit = absolute_mass_fit(M_inner, M_outer)

            params = rescale_parameters(lam, mu, t_gyr, sigma_m_sim)

            # Check t_evo constraint
            t_ok = params['t_evo_gyr'] <= T_ZOBS_GYR

            matches.append({
                'model': model_name,
                'snapshot_idx': i,
                'snapshot_time_gyr': t_gyr,
                'r2D_rs': r2D_rs,
                'mass_ratio': ratio,
                'n_sigma': n_sigma,
                'lambda': lam,
                'mu': mu,
                'mass_chi2': mass_fit['chi2'],
                'mass_nsigma_dof1': mass_fit['nsigma_dof1'],
                'mass_residual_inner_sigma': mass_fit['residual_inner_sigma'],
                'mass_residual_outer_sigma': mass_fit['residual_outer_sigma'],
                'r_s_kpc': params['r_s_kpc'],
                'rho0_msun_pc3': params['rho0_msun_pc3'],
                'sigma_m_cm2_g': params['sigma_m_cm2_g'],
                't_evo_gyr': params['t_evo_gyr'],
                't_ok': t_ok,
            })

    print(f"  [{model_name}] Found {n_within_3sigma} matching points "
          f"({len([m for m in matches if m['t_ok']])} with t_evo ≤ {T_ZOBS_GYR} Gyr)")
    return matches


# ----------------------------------------------------------------------
# Compute symmetry-breaking effect for velocity-dependent models
# ----------------------------------------------------------------------
def compute_symmetry_breaking(lam, mu, model_params, scale_v_kms):
    """For velocity-dependent models, compute how r_diss changes under rescaling.

    The rescaling v → sqrt(μ/λ) v changes the velocity scale.
    For massive emission, r_diss(v) has a characteristic v* = sqrt(2m/μ_red).
    After rescaling, the effective v in the halo changes by factor sqrt(μ/λ),
    so v/v* changes, and r_diss(v) changes.

    Returns the ratio r_diss_after / r_diss_before at the characteristic velocity.
    """
    # The velocity rescaling factor
    v_factor = np.sqrt(mu / lam)

    # For the model, compute r_diss at v_ref and at v_ref * v_factor
    # This requires importing the model
    sys.path.insert(0, os.path.join(_PROJ_ROOT, 'src', 'cross_sections'))
    from dsidm_models import sigma_T_born, r_diss

    # Reference velocity (dwarf scale)
    v_ref = 100.0  # km/s

    if model_params.m_mediator == 0 or model_params.m_mediator is None:
        # Massless: r_diss is constant, no symmetry breaking
        return 1.0, 1.0, v_factor
    else:
        # Massive emission: r_diss changes
        rd_before = r_diss(np.array([v_ref]), model_params)[0]
        rd_after = r_diss(np.array([v_ref * v_factor]), model_params)[0]
        ratio = rd_after / rd_before if rd_before > 0 else np.nan
        return rd_before, rd_after, v_factor
