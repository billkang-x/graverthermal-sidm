"""
Compute the effective σ_T/m for each dSIDM model at the velocity scales of
the observational constraints.

For velocity-dependent σ_T(v), the "effective" σ/m at a given halo velocity
dispersion T = <v^2> (3D) is the thermally-averaged transfer cross section:

    <σ_T v> / <v>  (velocity-weighted effective σ/m)

For the Born (Yukawa / Rutherford) models used in this project,
σ_T(v) ∝ 1/v⁴ at high v (massive mediator) or ∝ 1/v² log(v) at low v
(massless mediator, regulated). Hence <σ_T v> falls steeply with T.

We use the Maxwell-Boltzmann thermal averaging from src/fluid_runner/thermal_avg.py
to compute <σ_T v> / m_chi at the velocity scales of the four constraints.
"""
from __future__ import annotations

import os, sys
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJ_ROOT = os.path.normpath(os.path.join(_HERE, '..', '..'))
sys.path.insert(0, os.path.join(_PROJ_ROOT, 'src', 'cross_sections'))
sys.path.insert(0, os.path.join(_PROJ_ROOT, 'src', 'fluid_runner'))

from dsidm_models import (
    benchmark_models, sigma_T_born, r_diss as r_diss_func,
    DSIDMParameters, born_expansion_parameter, is_born_valid,
)
from thermal_avg import make_thermal_interpolators, mean_relative_speed


# Reference velocity dispersions (3D T = <v^2>) of each constraint, in (km/s)^2
# These are 1D velocity dispersions squared; the thermal_avg module uses
# T_km2_s2 = v_1D^2 directly.
CONSTRAINT_T = {
    'bullet':         3000.**2,    # ~9e6 (km/s)^2
    'cluster_cores':  1200.**2,    # ~1.4e6
    'dwarf_cores':    50.**2,      # ~2500
    'b1938':          80.**2,      # ~6400
}


def compute_effective_sigma_m(model_params: DSIDMParameters,
                              T_values_km2_s2: np.ndarray,
                              alpha_rescale: float = 1.0,
                              enforce_born: bool = True):
    """Compute <σ_T v>/m and <r_diss σ_T v>/m at each T value.

    alpha_rescale is an overall multiplicative factor on alpha_D, used to
    scan a family of models that share the same velocity dependence but
    differ in absolute normalization.
    """
    # Build a rescaled copy of the parameters
    if alpha_rescale != 1.0:
        p = DSIDMParameters(
            model=model_params.model,
            m_chi=model_params.m_chi,
            m_mediator=model_params.m_mediator,
            # alpha_rescale is a multiplicative factor on alpha_D. Squaring
            # it here would double-count the scan normalization.
            alpha_D=model_params.alpha_D * alpha_rescale,
            m_mediator_heavy=model_params.m_mediator_heavy,
            mediation=model_params.mediation,
            emission_type=model_params.emission_type,
        )
    else:
        p = model_params

    if enforce_born and not is_born_valid(p):
        shape = np.atleast_1d(T_values_km2_s2).shape
        return np.full(shape, np.nan), np.full(shape, np.nan)

    # sigma_v callable in km/s -> cm²/g
    sigma_v = lambda v: sigma_T_born(v, p)
    rdiss_v = lambda v: r_diss_func(v, p)

    # Thermal averaging over a v_axis up to ~5 sigma
    v_axis = np.logspace(-1, 4.5, 300)  # 0.1 to ~30000 km/s
    # make_thermal_interpolators returns (sig_avg_callable, rsig_avg_callable,
    # sig_arr, rsig_arr); we want the precomputed arrays sig_arr, rsig_arr
    # at the input T values.
    T_axis = np.atleast_1d(np.asarray(T_values_km2_s2, dtype=float))
    from thermal_avg import thermal_avg_sigma_v, thermal_avg_rsig_v
    sig_arr = np.array([thermal_avg_sigma_v(sigma_v, T) for T in T_axis])
    rsig_arr = np.array([thermal_avg_rsig_v(sigma_v, rdiss_v, T) for T in T_axis])

    # <σ_T v> is returned; to get an "effective" σ_T = <σ_T v> / <v>
    # we divide by the mean relative speed <v> = sqrt(8 T / pi)
    # (mean of 3D MB relative speed).
    v_mean = mean_relative_speed(T_axis)  # km/s
    sigma_eff = sig_arr / np.clip(v_mean, 1e-10, None)  # cm²/g
    rdiss_eff = np.where(sig_arr > 0, rsig_arr / np.clip(sig_arr, 1e-300, None), 1.0)
    return sigma_eff, rdiss_eff


def scan_alpha_for_model(model_key: str, model_params: DSIDMParameters,
                         alpha_log_range=(-3, 3), n_alpha=25):
    """Scan alpha_D over a logarithmic range and return effective σ/m at all
    constraint velocity scales for each alpha.

    Returns: dict with arrays 'alpha', 'sigma_eff' [n_alpha, n_constraints],
             'rdiss_eff' [n_alpha, n_constraints], and 'T' (constraint keys).
    """
    alpha_arr = np.logspace(alpha_log_range[0], alpha_log_range[1], n_alpha)
    T_arr = np.array([CONSTRAINT_T[k] for k in CONSTRAINT_T])
    constraint_keys = list(CONSTRAINT_T.keys())

    sigma_all = np.zeros((n_alpha, len(constraint_keys)))
    rdiss_all = np.zeros((n_alpha, len(constraint_keys)))
    born_parameter = np.zeros(n_alpha)

    for i, a in enumerate(alpha_arr):
        try:
            sig, rd = compute_effective_sigma_m(model_params, T_arr, alpha_rescale=a)
            sigma_all[i] = sig
            rdiss_all[i] = rd
            p_scaled = DSIDMParameters(
                model=model_params.model,
                m_chi=model_params.m_chi,
                m_mediator=model_params.m_mediator,
                alpha_D=model_params.alpha_D * a,
                m_mediator_heavy=model_params.m_mediator_heavy,
                mediation=model_params.mediation,
                emission_type=model_params.emission_type,
            )
            born_parameter[i] = born_expansion_parameter(p_scaled)
        except Exception as e:
            print(f"  [{model_key}] alpha={a:.2e} failed: {e}")
            sigma_all[i] = np.nan
            rdiss_all[i] = np.nan

    return {
        'alpha': alpha_arr,
        'T_keys': constraint_keys,
        'T_values': T_arr,
        'sigma_eff': sigma_all,
        'rdiss_eff': rdiss_all,
        'born_parameter': born_parameter,
    }


if __name__ == '__main__':
    bm = benchmark_models()
    T_arr = np.array([CONSTRAINT_T[k] for k in CONSTRAINT_T])
    keys = list(CONSTRAINT_T.keys())

    print(f"Constraints: {keys}")
    print(f"T values: {T_arr}")
    print()

    for name, p in bm.items():
        print(f"=== {name} ===")
        sig, rd = compute_effective_sigma_m(p, T_arr)
        for k, s, r in zip(keys, sig, rd):
            print(f"  {k:<15} T={CONSTRAINT_T[k]:.2e}  "
                  f"σ/m={s:.3e} cm²/g   r_diss={r:.4f}")
        print()
