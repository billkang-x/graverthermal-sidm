"""Host-orbit and tidal-radius priors for the B1938+666 perturber.

Powell et al. (2025) report ``r_t = 53 +/- 1 pc`` when the unknown 3D
galactocentric radius is set equal to its projected value.  Their freely
truncated pseudo-Jaffe fit instead gives ``r_t = 149 +/- 18 pc`` and
``m_tot = (2.82 +/- 0.26) 1e6 Msun``.  This module propagates the missing
line-of-sight position and an explicit pericentre fraction.  It intentionally
keeps the orbital families labelled as sensitivity priors; no velocity or
author posterior chain is publicly available.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np


@dataclass(frozen=True)
class PublishedTidalReference:
    """Published normalization and digitized macro-model geometry."""

    projected_radius_pc: float = 1520.0
    tidal_radius_pc: float = 53.0
    total_mass_msun: float = 1.54e6
    host_density_slope: float = 1.76

    def to_dict(self) -> dict[str, float]:
        return asdict(self)


@dataclass(frozen=True)
class FreePseudoJaffeSummary:
    """Gaussian summary of the public free-truncation result."""

    total_mass_msun: float = 2.82e6
    total_mass_sigma_msun: float = 0.26e6
    tidal_radius_pc: float = 149.0
    tidal_radius_sigma_pc: float = 18.0


ORBIT_PRIORS = {
    "circular_upper_envelope": None,
    "phase_mixed": (4.0, 2.0),
    "radial_sensitivity": (2.0, 4.0),
}


def tidal_radius_power_law(
    pericentre_pc: np.ndarray | float,
    total_mass_msun: np.ndarray | float,
    *,
    reference: PublishedTidalReference = PublishedTidalReference(),
    host_density_slope: np.ndarray | float | None = None,
) -> np.ndarray:
    """Scale the published tidal radius through a power-law host.

    Powell et al. use ``r_t proportional to R [m/M(<R)]^(1/3)``.  For a
    host with 3D density slope ``gamma``, ``M(<R) proportional to R^(3-gamma)``
    and hence ``r_t proportional to m^(1/3) R^(gamma/3)``.
    """
    pericentre = np.asarray(pericentre_pc, dtype=float)
    total_mass = np.asarray(total_mass_msun, dtype=float)
    gamma = np.asarray(
        reference.host_density_slope
        if host_density_slope is None else host_density_slope,
        dtype=float,
    )
    if (
        np.any(~np.isfinite(pericentre)) or np.any(pericentre <= 0)
        or np.any(~np.isfinite(total_mass)) or np.any(total_mass <= 0)
        or np.any(~np.isfinite(gamma)) or np.any(gamma <= 0)
    ):
        raise ValueError("pericentre, mass, and host slope must be positive")
    return (
        reference.tidal_radius_pc
        * (total_mass / reference.total_mass_msun) ** (1.0 / 3.0)
        * (pericentre / reference.projected_radius_pc) ** (gamma / 3.0)
    )


def minimum_current_radius_pc(
    target_tidal_radius_pc: np.ndarray | float,
    total_mass_msun: np.ndarray | float,
    *,
    reference: PublishedTidalReference = PublishedTidalReference(),
    host_density_slope: np.ndarray | float | None = None,
) -> np.ndarray:
    """Orbit-independent radius lower bound from ``r_peri <= r_current``."""
    target = np.asarray(target_tidal_radius_pc, dtype=float)
    mass = np.asarray(total_mass_msun, dtype=float)
    gamma = np.asarray(
        reference.host_density_slope
        if host_density_slope is None else host_density_slope,
        dtype=float,
    )
    if np.any(target <= 0) or np.any(~np.isfinite(target)):
        raise ValueError("target tidal radius must be finite and positive")
    normalization = reference.tidal_radius_pc * (
        mass / reference.total_mass_msun
    ) ** (1.0 / 3.0)
    return reference.projected_radius_pc * (target / normalization) ** (3.0 / gamma)


def nfw_tracer_density(
    radius_pc: np.ndarray | float,
    scale_radius_pc: float,
) -> np.ndarray:
    """Unnormalized NFW-like number density used as a transparent sensitivity prior."""
    radius = np.asarray(radius_pc, dtype=float)
    if scale_radius_pc <= 0 or np.any(radius <= 0) or np.any(~np.isfinite(radius)):
        raise ValueError("radii and scale radius must be finite and positive")
    x = radius / scale_radius_pc
    return 1.0 / (x * (1.0 + x) ** 2)


def weighted_quantile(
    values: np.ndarray,
    quantiles: np.ndarray | tuple[float, ...],
    weights: np.ndarray,
) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    quantiles = np.asarray(quantiles, dtype=float)
    weights = np.asarray(weights, dtype=float)
    if values.ndim != 1 or weights.shape != values.shape:
        raise ValueError("values and weights must be one-dimensional and aligned")
    if np.any(~np.isfinite(values)) or np.any(~np.isfinite(weights)):
        raise ValueError("values and weights must be finite")
    if np.any(weights < 0) or np.sum(weights) <= 0:
        raise ValueError("weights must be non-negative with positive sum")
    if np.any((quantiles < 0) | (quantiles > 1)):
        raise ValueError("quantiles must lie in [0, 1]")
    order = np.argsort(values)
    sorted_values = values[order]
    sorted_weights = weights[order]
    cumulative = np.cumsum(sorted_weights) - 0.5 * sorted_weights
    cumulative /= np.sum(sorted_weights)
    return np.interp(quantiles, cumulative, sorted_values)


def importance_sample_orbit_posterior(
    sample_count: int,
    orbit_prior: str,
    *,
    seed: int = 1938666,
    reference: PublishedTidalReference = PublishedTidalReference(),
    imaging: FreePseudoJaffeSummary = FreePseudoJaffeSummary(),
    host_scale_radius_pc: float = 30_000.0,
    host_max_radius_pc: float = 200_000.0,
    host_slope_sigma: float = 0.02,
) -> dict[str, np.ndarray | float | str]:
    """Importance sample a labelled orbit sensitivity posterior.

    The proposal is uniform in line-of-sight coordinate.  It is reweighted to
    an NFW-like tracer density conditional on the observed projected radius,
    and to the Gaussian public measurement of the free truncation radius.
    The public mass summary is sampled as an independent Gaussian because its
    covariance with radius has not been released.
    """
    if sample_count < 10_000:
        raise ValueError("sample_count must be at least 10000")
    if orbit_prior not in ORBIT_PRIORS:
        raise ValueError(f"unknown orbit prior: {orbit_prior}")
    if host_max_radius_pc <= reference.projected_radius_pc:
        raise ValueError("host maximum radius must exceed projected radius")
    rng = np.random.default_rng(seed)
    z_max = np.sqrt(
        host_max_radius_pc**2 - reference.projected_radius_pc**2
    )
    line_of_sight = rng.uniform(-z_max, z_max, sample_count)
    current_radius = np.sqrt(
        reference.projected_radius_pc**2 + line_of_sight**2
    )
    beta_parameters = ORBIT_PRIORS[orbit_prior]
    if beta_parameters is None:
        pericentre_fraction = np.ones(sample_count)
    else:
        pericentre_fraction = rng.beta(*beta_parameters, sample_count)
    pericentre = current_radius * pericentre_fraction
    total_mass = rng.normal(
        imaging.total_mass_msun, imaging.total_mass_sigma_msun, sample_count
    )
    total_mass = np.maximum(total_mass, 1.0)
    host_slope = rng.normal(
        reference.host_density_slope, host_slope_sigma, sample_count
    )
    predicted_tidal_radius = tidal_radius_power_law(
        pericentre, total_mass, reference=reference,
        host_density_slope=host_slope,
    )
    radial_weight = nfw_tracer_density(current_radius, host_scale_radius_pc)
    tidal_log_weight = -0.5 * (
        (predicted_tidal_radius - imaging.tidal_radius_pc)
        / imaging.tidal_radius_sigma_pc
    ) ** 2
    log_weight = np.log(radial_weight) + tidal_log_weight
    log_weight -= float(np.max(log_weight))
    weight = np.exp(log_weight)
    weight /= np.sum(weight)
    effective_sample_size = 1.0 / float(np.sum(weight**2))
    return {
        "orbit_prior": orbit_prior,
        "host_scale_radius_pc": float(host_scale_radius_pc),
        "host_max_radius_pc": float(host_max_radius_pc),
        "line_of_sight_pc": line_of_sight,
        "current_radius_pc": current_radius,
        "pericentre_fraction": pericentre_fraction,
        "pericentre_pc": pericentre,
        "total_mass_msun": total_mass,
        "host_density_slope": host_slope,
        "predicted_tidal_radius_pc": predicted_tidal_radius,
        "weight": weight,
        "effective_sample_size": effective_sample_size,
    }
