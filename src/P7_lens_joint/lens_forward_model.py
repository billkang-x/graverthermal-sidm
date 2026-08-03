"""Independent visibility-plane lensing operators for B1938+666.

This implements the public mathematical model, not the private PRONTO
configuration.  In particular, it provides the elliptical power-law macro
lens, the power-law m=3 and m=4 convergence multipoles used by Powell et al.,
external shear, spherical perturbers, a bilinear pixel-source operator, and a
band-integrated NUFFT measurement operator.
"""

from __future__ import annotations

from dataclasses import dataclass

import finufft
import numpy as np
from lenstronomy.LensModel.Profiles.epl import EPL
from scipy import sparse
from scipy.linalg import cho_factor, cho_solve

from visibility_likelihood import subchannel_frequencies_hz


ARCSEC_TO_RADIAN = np.deg2rad(1.0 / 3600.0)
PUBLIC_REPRODUCTION_GAPS = (
    "exact custom RFI intervals",
    "registration between PRONTO model coordinates and the UVFITS phase centre",
    "adaptive source grid and source regularization matrix",
    "numerical macro-model prior and posterior chain",
    "PRONTO fast-chi2 and log-determinant preconditioner",
    "MultiNest settings and posterior samples",
)


@dataclass(frozen=True)
class MacroLensParameters:
    theta_e_arcsec: float
    gamma: float
    axis_ratio: float
    position_angle_radian: float
    center_x_arcsec: float
    center_y_arcsec: float
    shear: float = 0.0
    shear_angle_radian: float = 0.0
    a3: float = 0.0
    b3: float = 0.0
    a4: float = 0.0
    b4: float = 0.0


@dataclass(frozen=True)
class ImageGrid:
    nx: int
    ny: int
    pixel_scale_arcsec: float
    center_x_arcsec: float = 0.0
    center_y_arcsec: float = 0.0

    def __post_init__(self) -> None:
        if self.nx <= 0 or self.ny <= 0 or self.pixel_scale_arcsec <= 0:
            raise ValueError("grid dimensions and pixel scale must be positive")

    @property
    def pixel_area_arcsec2(self) -> float:
        return self.pixel_scale_arcsec**2

    def coordinates(self) -> tuple[np.ndarray, np.ndarray]:
        x = (
            np.arange(self.nx, dtype=float) - self.nx // 2
        ) * self.pixel_scale_arcsec + self.center_x_arcsec
        y = (
            np.arange(self.ny, dtype=float) - self.ny // 2
        ) * self.pixel_scale_arcsec + self.center_y_arcsec
        xx, yy = np.meshgrid(x, y)
        return xx, yy


@dataclass(frozen=True)
class SourceGrid:
    nx: int
    ny: int
    pixel_scale_arcsec: float
    center_x_arcsec: float = 0.0
    center_y_arcsec: float = 0.0

    def __post_init__(self) -> None:
        if self.nx < 2 or self.ny < 2 or self.pixel_scale_arcsec <= 0:
            raise ValueError("source grid needs at least 2x2 positive-scale pixels")

    @property
    def x_min_arcsec(self) -> float:
        return self.center_x_arcsec - (self.nx // 2) * self.pixel_scale_arcsec

    @property
    def y_min_arcsec(self) -> float:
        return self.center_y_arcsec - (self.ny // 2) * self.pixel_scale_arcsec


def _ellipticity_components(axis_ratio: float, angle: float) -> tuple[float, float]:
    if not 0 < axis_ratio <= 1:
        raise ValueError("axis ratio must lie in (0, 1]")
    modulus = (1.0 - axis_ratio) / (1.0 + axis_ratio)
    return modulus * np.cos(2.0 * angle), modulus * np.sin(2.0 * angle)


def power_law_multipole_deflection(
    x_arcsec: np.ndarray,
    y_arcsec: np.ndarray,
    *,
    gamma: float,
    order: int,
    sine_coefficient: float,
    cosine_coefficient: float,
    center_x_arcsec: float = 0.0,
    center_y_arcsec: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Deflection for Powell et al. equation (8)."""
    if order < 1 or not 1.0 < gamma < 3.0:
        raise ValueError("multipole order and host slope are not physical")
    x, y = np.broadcast_arrays(
        np.asarray(x_arcsec, dtype=float) - center_x_arcsec,
        np.asarray(y_arcsec, dtype=float) - center_y_arcsec,
    )
    radius = np.hypot(x, y)
    safe_radius = np.maximum(radius, 1.0e-12)
    angle = np.arctan2(y, x)
    radial_power = 3.0 - gamma
    denominator = radial_power**2 - order**2
    if abs(denominator) < 1.0e-12:
        raise ValueError("resonant multipole slope is unsupported")
    angular = (
        sine_coefficient * np.sin(order * angle)
        + cosine_coefficient * np.cos(order * angle)
    )
    angular_derivative = order * (
        sine_coefficient * np.cos(order * angle)
        - cosine_coefficient * np.sin(order * angle)
    )
    common = 2.0 * safe_radius ** (radial_power - 1.0) / denominator
    alpha_r = common * radial_power * angular
    alpha_angle = common * angular_derivative
    alpha_x = alpha_r * np.cos(angle) - alpha_angle * np.sin(angle)
    alpha_y = alpha_r * np.sin(angle) + alpha_angle * np.cos(angle)
    alpha_x = np.where(radius > 0, alpha_x, 0.0)
    alpha_y = np.where(radius > 0, alpha_y, 0.0)
    return alpha_x, alpha_y


def macro_deflection(
    x_arcsec: np.ndarray,
    y_arcsec: np.ndarray,
    parameters: MacroLensParameters,
) -> tuple[np.ndarray, np.ndarray]:
    """EPL, m=3/4 convergence multipoles, and external shear."""
    e1, e2 = _ellipticity_components(
        parameters.axis_ratio, parameters.position_angle_radian
    )
    epl = EPL()
    alpha_x, alpha_y = epl.derivatives(
        np.asarray(x_arcsec, dtype=float), np.asarray(y_arcsec, dtype=float),
        parameters.theta_e_arcsec, parameters.gamma, e1, e2,
        parameters.center_x_arcsec, parameters.center_y_arcsec,
    )
    for order, sine_coefficient, cosine_coefficient in (
        (3, parameters.a3, parameters.b3),
        (4, parameters.a4, parameters.b4),
    ):
        multipole_x, multipole_y = power_law_multipole_deflection(
            x_arcsec, y_arcsec, gamma=parameters.gamma, order=order,
            sine_coefficient=sine_coefficient,
            cosine_coefficient=cosine_coefficient,
            center_x_arcsec=parameters.center_x_arcsec,
            center_y_arcsec=parameters.center_y_arcsec,
        )
        alpha_x = alpha_x + multipole_x
        alpha_y = alpha_y + multipole_y
    shifted_x = np.asarray(x_arcsec, dtype=float) - parameters.center_x_arcsec
    shifted_y = np.asarray(y_arcsec, dtype=float) - parameters.center_y_arcsec
    shear_1 = parameters.shear * np.cos(2.0 * parameters.shear_angle_radian)
    shear_2 = parameters.shear * np.sin(2.0 * parameters.shear_angle_radian)
    alpha_x = alpha_x + shear_1 * shifted_x + shear_2 * shifted_y
    alpha_y = alpha_y + shear_2 * shifted_x - shear_1 * shifted_y
    return np.asarray(alpha_x), np.asarray(alpha_y)


def pseudo_jaffe_projected_mass(
    radius_pc: np.ndarray | float,
    total_mass_msun: float,
    truncation_radius_pc: float,
) -> np.ndarray:
    """Cylindrical mass for the singular pseudo-Jaffe profile in the paper."""
    radius = np.asarray(radius_pc, dtype=float)
    if total_mass_msun <= 0 or truncation_radius_pc <= 0 or np.any(radius < 0):
        raise ValueError("mass and truncation radius must be positive")
    x = radius / truncation_radius_pc
    return total_mass_msun * (1.0 + x - np.sqrt(1.0 + x**2))


def spherical_mass_deflection(
    x_arcsec: np.ndarray,
    y_arcsec: np.ndarray,
    *,
    center_x_arcsec: float,
    center_y_arcsec: float,
    projected_mass_function,
    pc_per_arcsec: float,
    sigma_critical_msun_arcsec2: float = 1.50e11,
) -> tuple[np.ndarray, np.ndarray]:
    """Deflection from an arbitrary cylindrical enclosed-mass function."""
    if pc_per_arcsec <= 0 or sigma_critical_msun_arcsec2 <= 0:
        raise ValueError("distance scale and critical density must be positive")
    dx = np.asarray(x_arcsec, dtype=float) - center_x_arcsec
    dy = np.asarray(y_arcsec, dtype=float) - center_y_arcsec
    angular_radius = np.hypot(dx, dy)
    physical_radius = angular_radius * pc_per_arcsec
    enclosed_mass = np.asarray(projected_mass_function(physical_radius), dtype=float)
    safe_radius = np.maximum(angular_radius, 1.0e-15)
    alpha_radius = enclosed_mass / (
        np.pi * sigma_critical_msun_arcsec2 * safe_radius
    )
    alpha_x = alpha_radius * dx / safe_radius
    alpha_y = alpha_radius * dy / safe_radius
    return (
        np.where(angular_radius > 0, alpha_x, 0.0),
        np.where(angular_radius > 0, alpha_y, 0.0),
    )


def ray_shoot(
    x_arcsec: np.ndarray,
    y_arcsec: np.ndarray,
    macro: MacroLensParameters,
    *,
    perturber_deflections: tuple[tuple[np.ndarray, np.ndarray], ...] = (),
) -> tuple[np.ndarray, np.ndarray]:
    alpha_x, alpha_y = macro_deflection(x_arcsec, y_arcsec, macro)
    for perturber_x, perturber_y in perturber_deflections:
        alpha_x = alpha_x + perturber_x
        alpha_y = alpha_y + perturber_y
    return np.asarray(x_arcsec) - alpha_x, np.asarray(y_arcsec) - alpha_y


def bilinear_source_operator(
    beta_x_arcsec: np.ndarray,
    beta_y_arcsec: np.ndarray,
    source_grid: SourceGrid,
    *,
    image_pixel_area_arcsec2: float = 1.0,
) -> sparse.csr_matrix:
    """Map source surface brightness to image-pixel flux."""
    beta_x, beta_y = np.broadcast_arrays(
        np.asarray(beta_x_arcsec, dtype=float), np.asarray(beta_y_arcsec, dtype=float)
    )
    flat_x = beta_x.ravel()
    flat_y = beta_y.ravel()
    fractional_x = (flat_x - source_grid.x_min_arcsec) / source_grid.pixel_scale_arcsec
    fractional_y = (flat_y - source_grid.y_min_arcsec) / source_grid.pixel_scale_arcsec
    lower_x = np.floor(fractional_x).astype(int)
    lower_y = np.floor(fractional_y).astype(int)
    valid = (
        (lower_x >= 0) & (lower_x < source_grid.nx - 1)
        & (lower_y >= 0) & (lower_y < source_grid.ny - 1)
    )
    rows = np.flatnonzero(valid)
    dx = fractional_x[valid] - lower_x[valid]
    dy = fractional_y[valid] - lower_y[valid]
    columns = []
    weights = []
    for offset_x, offset_y, weight in (
        (0, 0, (1.0 - dx) * (1.0 - dy)),
        (1, 0, dx * (1.0 - dy)),
        (0, 1, (1.0 - dx) * dy),
        (1, 1, dx * dy),
    ):
        columns.append(
            (lower_y[valid] + offset_y) * source_grid.nx
            + lower_x[valid] + offset_x
        )
        weights.append(weight * image_pixel_area_arcsec2)
    return sparse.coo_matrix(
        (
            np.concatenate(weights),
            (np.tile(rows, 4), np.concatenate(columns)),
        ),
        shape=(flat_x.size, source_grid.nx * source_grid.ny),
    ).tocsr()


def gradient_precision(
    source_grid: SourceGrid,
    regularization_strength: float,
    *,
    ridge_fraction: float = 1.0e-8,
) -> sparse.csr_matrix:
    """Positive-definite first-difference source-prior precision."""
    if regularization_strength <= 0 or ridge_fraction <= 0:
        raise ValueError("regularization and ridge must be positive")
    row_count = source_grid.ny * (source_grid.nx - 1)
    horizontal = sparse.lil_matrix(
        (row_count, source_grid.nx * source_grid.ny), dtype=float
    )
    row = 0
    for iy in range(source_grid.ny):
        for ix in range(source_grid.nx - 1):
            left = iy * source_grid.nx + ix
            horizontal[row, left] = -1.0
            horizontal[row, left + 1] = 1.0
            row += 1
    row_count = (source_grid.ny - 1) * source_grid.nx
    vertical = sparse.lil_matrix(
        (row_count, source_grid.nx * source_grid.ny), dtype=float
    )
    row = 0
    for iy in range(source_grid.ny - 1):
        for ix in range(source_grid.nx):
            lower = iy * source_grid.nx + ix
            vertical[row, lower] = -1.0
            vertical[row, lower + source_grid.nx] = 1.0
            row += 1
    difference = sparse.vstack((horizontal.tocsr(), vertical.tocsr()))
    identity = sparse.eye(source_grid.nx * source_grid.ny, format="csr")
    return regularization_strength * (
        difference.T @ difference + ridge_fraction * identity
    )


class BandIntegratedNufft:
    """NUFFT from a regular image-pixel flux grid to averaged visibilities."""

    def __init__(
        self,
        image_grid: ImageGrid,
        *,
        bandwidth_hz: float = 8.0e6,
        channel_count: int = 32,
        tolerance: float = 1.0e-8,
    ) -> None:
        if tolerance <= 0:
            raise ValueError("NUFFT tolerance must be positive")
        self.image_grid = image_grid
        self.bandwidth_hz = bandwidth_hz
        self.channel_count = channel_count
        self.tolerance = tolerance

    def _coordinates(
        self, uu_seconds: np.ndarray, vv_seconds: np.ndarray, frequency_hz: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        pixel_radians = self.image_grid.pixel_scale_arcsec * ARCSEC_TO_RADIAN
        return (
            2.0 * np.pi * uu_seconds * frequency_hz * pixel_radians,
            2.0 * np.pi * vv_seconds * frequency_hz * pixel_radians,
        )

    def forward(
        self,
        image_flux_jy: np.ndarray,
        uu_seconds: np.ndarray,
        vv_seconds: np.ndarray,
        centre_hz: np.ndarray,
    ) -> np.ndarray:
        image = np.asarray(image_flux_jy, dtype=np.complex128)
        expected = (self.image_grid.ny, self.image_grid.nx)
        if image.shape != expected:
            raise ValueError(f"image must have shape {expected}")
        uu, vv, centre = np.broadcast_arrays(
            np.asarray(uu_seconds, dtype=float),
            np.asarray(vv_seconds, dtype=float),
            np.asarray(centre_hz, dtype=float),
        )
        frequencies = subchannel_frequencies_hz(
            centre, bandwidth_hz=self.bandwidth_hz,
            channel_count=self.channel_count,
        )
        output = np.zeros(uu.shape, dtype=np.complex128)
        centre_x_rad = self.image_grid.center_x_arcsec * ARCSEC_TO_RADIAN
        centre_y_rad = self.image_grid.center_y_arcsec * ARCSEC_TO_RADIAN
        modes = np.ascontiguousarray(image.T)
        for channel in range(self.channel_count):
            frequency = frequencies[..., channel]
            x, y = self._coordinates(uu, vv, frequency)
            shifted = finufft.nufft2d2(
                x.ravel(), y.ravel(), modes, isign=-1, eps=self.tolerance
            ).reshape(uu.shape)
            phase = np.exp(
                -2j * np.pi * frequency
                * (uu * centre_x_rad + vv * centre_y_rad)
            )
            output += phase * shifted
        return output / self.channel_count

    def adjoint(
        self,
        visibility_values: np.ndarray,
        uu_seconds: np.ndarray,
        vv_seconds: np.ndarray,
        centre_hz: np.ndarray,
    ) -> np.ndarray:
        values, uu, vv, centre = np.broadcast_arrays(
            np.asarray(visibility_values, dtype=np.complex128),
            np.asarray(uu_seconds, dtype=float),
            np.asarray(vv_seconds, dtype=float),
            np.asarray(centre_hz, dtype=float),
        )
        frequencies = subchannel_frequencies_hz(
            centre, bandwidth_hz=self.bandwidth_hz,
            channel_count=self.channel_count,
        )
        output = np.zeros(
            (self.image_grid.nx, self.image_grid.ny), dtype=np.complex128
        )
        centre_x_rad = self.image_grid.center_x_arcsec * ARCSEC_TO_RADIAN
        centre_y_rad = self.image_grid.center_y_arcsec * ARCSEC_TO_RADIAN
        for channel in range(self.channel_count):
            frequency = frequencies[..., channel]
            x, y = self._coordinates(uu, vv, frequency)
            phase_conjugate = np.exp(
                2j * np.pi * frequency
                * (uu * centre_x_rad + vv * centre_y_rad)
            )
            output += finufft.nufft2d1(
                x.ravel(), y.ravel(), (values * phase_conjugate).ravel(),
                (self.image_grid.nx, self.image_grid.ny),
                isign=1, eps=self.tolerance,
            )
        return (output / self.channel_count).T


@dataclass(frozen=True)
class MarginalizedSourceResult:
    log_likelihood: float
    source_map: np.ndarray
    chi2: float
    prior_penalty: float
    logdet_hessian: float
    logdet_prior: float


def marginalized_source_log_likelihood_explicit(
    observed: np.ndarray,
    design_matrix: np.ndarray,
    component_sigma: np.ndarray,
    prior_precision: np.ndarray | sparse.spmatrix,
) -> MarginalizedSourceResult:
    """Exact Gaussian source marginalization for a finite design matrix."""
    data = np.asarray(observed, dtype=np.complex128)
    design = np.asarray(design_matrix, dtype=np.complex128)
    sigma = np.asarray(component_sigma, dtype=float)
    prior = (
        prior_precision.toarray()
        if sparse.issparse(prior_precision) else np.asarray(prior_precision, dtype=float)
    )
    if design.ndim != 2 or data.shape != (design.shape[0],):
        raise ValueError("data and design-matrix dimensions do not align")
    if sigma.shape != data.shape or np.any(sigma <= 0) or not np.isfinite(sigma).all():
        raise ValueError("component sigma must be positive and match the data")
    if prior.shape != (design.shape[1], design.shape[1]):
        raise ValueError("prior precision has the wrong shape")
    inverse_variance = sigma**-2
    weighted_design = design * inverse_variance[:, None]
    hessian = np.real(design.conj().T @ weighted_design) + prior
    right_hand_side = np.real(design.conj().T @ (inverse_variance * data))
    factor = cho_factor(hessian, lower=True, check_finite=True)
    source_map = cho_solve(factor, right_hand_side)
    residual = data - design @ source_map
    chi2 = float(np.sum(np.abs(residual) ** 2 * inverse_variance))
    prior_penalty = float(source_map @ prior @ source_map)
    sign_hessian, logdet_hessian = np.linalg.slogdet(hessian)
    sign_prior, logdet_prior = np.linalg.slogdet(prior)
    if sign_hessian <= 0 or sign_prior <= 0:
        raise ValueError("prior and posterior precision must be positive definite")
    normalization = float(np.sum(np.log(2.0 * np.pi * sigma**2)))
    log_likelihood = -normalization - 0.5 * (
        chi2 + prior_penalty + logdet_hessian - logdet_prior
    )
    return MarginalizedSourceResult(
        log_likelihood=float(log_likelihood),
        source_map=source_map,
        chi2=chi2,
        prior_penalty=prior_penalty,
        logdet_hessian=float(logdet_hessian),
        logdet_prior=float(logdet_prior),
    )
