"""Born-limit energy-loss kernel for identical fermion scattering.

This module implements the long-range identical-particle energy-differential
cross sections in Lankester-Broche & Pradler (2026), Eqs. (3.43) and (3.47),
and converts their energy integral into the effective cooling cross section
needed by the gravothermal fluid solver.

The implementation is deliberately limited to the two fermionic channels used
by this project: ``chi-V`` and ``chi-phi``.  Scalar and short-range channels
require their separate matrix elements and are rejected explicitly.
"""

from __future__ import annotations

import numpy as np

from dsidm_models import DSIDMParameters, GEV2_TO_CM2G, KM_S_TO_NAT


def _K(x, a: float, b: float):
    return a + b * x**2 / (2.0 - x) ** 2


def _L(x, a: float, b: float, c: float):
    x = np.asarray(x, dtype=float)
    s = np.sqrt(np.clip(1.0 - x, 0.0, None))
    numerator = a * (2.0 - x) ** 4 + b * (2.0 - x) ** 2 * x**2 + c * x**4
    denominator = (2.0 - x) ** 3 * np.clip(s, 1e-14, None)
    # This form avoids cancellation for x close to zero while preserving the
    # logarithmic endpoint structure of the published expression.
    log_ratio = np.log1p(2.0 * s / np.clip(1.0 - s, 1e-14, None))
    return numerator / denominator * log_ratio


def omega_dsigma_domega_long_range_identical(
    x,
    x_min: float,
    p: DSIDMParameters,
):
    """Return ``omega * d sigma/d omega`` in natural units (GeV^-2).

    ``x = omega / E_cm = 2 omega/(mu v_i^2)`` and
    ``x_min = m_emitted/E_cm``.  The formula includes the massive-emission
    correction ``kappa = m_emitted/omega = x_min/x``.
    """
    if p.mediation != "long":
        raise ValueError("the implemented emission kernel is long-range only")
    if p.model not in {"chi-V", "chi-phi", "chi_tilde-phi"}:
        raise ValueError("kernel supports only identical fermion chi-V/chi-phi channels")

    x = np.asarray(x, dtype=float)
    x = np.clip(x, max(float(x_min), 1e-14), 1.0 - 1e-12)
    x_min = float(np.clip(x_min, 0.0, 1.0))
    kappa = x_min / x
    phase_x = np.sqrt(np.clip(1.0 - x, 0.0, None))
    phase_kappa = np.sqrt(np.clip(1.0 - kappa**2, 0.0, None))

    if p.model == "chi-V":
        coupling_power = p.g_eff**6
        bracket = (
            _K(x, 17.0, -3.0)
            + 0.5 * _L(x, 12.0, -7.0, -3.0)
            + (kappa**2 / 4.0)
            * (_K(x, 52.0, -8.0) + _L(x, 16.0, -11.0, -4.0))
        )
        mass_factor = phase_kappa**3
    elif p.model in {"chi-phi", "chi_tilde-phi"}:
        coupling_power = p.g_eff**6
        bracket = (
            _K(x, 18.0, -2.0)
            + _L(x, 4.0, -4.0, -1.0)
            + kappa**2 * (_K(x, -16.0, 4.0) + _L(x, -8.0, 3.0, 2.0))
            + (kappa**4 / 4.0)
            * (_K(x, 52.0, -8.0) + _L(x, 16.0, -11.0, -4.0))
        )
        mass_factor = phase_kappa
    else:  # pragma: no cover - guarded above
        raise AssertionError("unreachable channel")

    prefactor = coupling_power * phase_x / (240.0 * np.pi**3 * p.m_chi**2)
    result = prefactor * mass_factor * bracket
    # Roundoff near x=1 can generate tiny negative values although the
    # differential rate is non-negative.  Large negative values indicate a
    # formula or parameter error and are not silently hidden.
    if np.any(result < -1e-10 * np.maximum(np.abs(prefactor), 1e-300)):
        raise FloatingPointError("negative energy-differential cross section")
    return np.clip(result, 0.0, None)


def radiated_energy_cross_section(
    v_km_s: float,
    p: DSIDMParameters,
    n_quadrature: int = 192,
) -> float:
    """Return ``integral d omega omega d sigma/d omega`` in GeV^-1."""
    v_nat = float(v_km_s) * KM_S_TO_NAT
    if v_nat <= 0:
        return 0.0
    e_cm = 0.5 * p.mu * v_nat**2
    emitted_mass = 0.0 if p.emission_type == "massless" else p.m_mediator
    if emitted_mass is None or emitted_mass < 0:
        raise ValueError("emitted mass must be non-negative")
    if e_cm <= emitted_mass:
        return 0.0

    x_min = emitted_mass / e_cm
    nodes, weights = np.polynomial.legendre.leggauss(n_quadrature)
    x = x_min + 0.5 * (nodes + 1.0) * (1.0 - x_min)
    kernel = omega_dsigma_domega_long_range_identical(x, x_min, p)
    return float(e_cm * 0.5 * (1.0 - x_min) * np.dot(weights, kernel))


def microscopic_cooling_sigma_m(v_km_s, p: DSIDMParameters):
    """Return the per-speed cooling combination in ``cm^2/g``.

    Matching the binary energy-loss rate to the Schmidt et al. fluid form gives

        sigma_cool/m = 2 Q(v) / (m_chi^2 v^2),

    where ``Q(v) = integral d omega omega d sigma/d omega``.  This quantity is
    already the product of the elastic transfer cross section and ``r_diss-1``;
    it must therefore be thermally averaged with the energy-weighted v^3
    moment, not reconstructed from two separate collision averages.
    """
    velocities = np.atleast_1d(np.asarray(v_km_s, dtype=float))
    output = np.zeros_like(velocities)
    for i, velocity in enumerate(velocities):
        v_nat = velocity * KM_S_TO_NAT
        if v_nat <= 0:
            continue
        q_energy = radiated_energy_cross_section(float(velocity), p)
        output[i] = 2.0 * q_energy / (p.m_chi**2 * v_nat**2) * GEV2_TO_CM2G
    return output

