"""Pure functions for the dissipative fluid cooling closure.

The source fluid equation gives the volumetric loss rate

    C_vol = (8/sqrt(pi)) (sigma_T/m_chi) (r_diss - 1) rho^2 nu^3.

After division by rho, ``specific_cooling_rate`` returns the specific
internal-energy loss rate used by :mod:`dissipative_halo`.  The
``specific_cooling_rate_from_moment`` companion accepts the exact
energy-weighted emission moment after its relative-speed Maxwell average.
"""

from __future__ import annotations

import numpy as np


COOLING_PREFACTOR = 8.0 / np.sqrt(np.pi)


def specific_cooling_rate(
    sigma_m,
    rho,
    nu,
    rdiss,
    dissipation_prefactor: float = 1.0,
):
    """Return ``C_vol/rho`` in the code's specific-energy units per time."""
    sigma_m = np.asarray(sigma_m, dtype=float)
    rho = np.asarray(rho, dtype=float)
    nu = np.asarray(nu, dtype=float)
    rdiss = np.asarray(rdiss, dtype=float)
    if not np.isfinite(dissipation_prefactor) or dissipation_prefactor < 0:
        raise ValueError("dissipation_prefactor must be finite and non-negative")
    if np.any(sigma_m < 0) or np.any(rho < 0) or np.any(nu < 0):
        raise ValueError("sigma_m, rho, and nu must be non-negative")
    return (
        dissipation_prefactor
        * COOLING_PREFACTOR
        * sigma_m
        * rho
        * nu**3
        * np.clip(rdiss - 1.0, 0.0, None)
    )


def specific_cooling_rate_from_moment(
    sigma_cooling_m,
    rho,
    nu,
    dissipation_prefactor: float = 1.0,
):
    """Return the specific rate from the energy-weighted effective cross section.

    ``sigma_cooling_m`` is defined by

        (1/4) <sigma_m(v) [r_diss(v)-1] v^3>
        = (8/sqrt(pi)) sigma_cooling_m(T) T^(3/2).

    This form preserves the exact velocity weighting while keeping the
    gravothermal solver's code-unit interface unchanged.
    """
    sigma_cooling_m = np.asarray(sigma_cooling_m, dtype=float)
    rho = np.asarray(rho, dtype=float)
    nu = np.asarray(nu, dtype=float)
    if not np.isfinite(dissipation_prefactor) or dissipation_prefactor < 0:
        raise ValueError("dissipation_prefactor must be finite and non-negative")
    if np.any(sigma_cooling_m < 0) or np.any(rho < 0) or np.any(nu < 0):
        raise ValueError("sigma_cooling_m, rho, and nu must be non-negative")
    return (
        dissipation_prefactor
        * COOLING_PREFACTOR
        * sigma_cooling_m
        * rho
        * nu**3
    )
