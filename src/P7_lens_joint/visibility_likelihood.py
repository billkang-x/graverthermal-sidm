"""Noise, flagging, and complex likelihood utilities for GM068 visibilities.

The public pipeline UVFITS contains one channel per 8 MHz IF because AIPS
``SPLIT`` averaged the original 32 channels.  Sky models must therefore be
integrated across those 32 channel centres rather than evaluated only at the
reported IF centre.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


PUBLISHED_TRIANGLE_CODES = frozenset({258, 259, 515})
PUBLISHED_RFI_BASELINE_CODES = frozenset({3350, 4113, 4374, 4886})


def subchannel_frequencies_hz(
    centre_hz: np.ndarray | float,
    *,
    bandwidth_hz: float = 8.0e6,
    channel_count: int = 32,
) -> np.ndarray:
    """Return uniformly spaced channel centres for an averaged IF."""
    if bandwidth_hz <= 0:
        raise ValueError("bandwidth_hz must be positive")
    if channel_count <= 0:
        raise ValueError("channel_count must be positive")
    centre = np.asarray(centre_hz, dtype=float)
    if not np.isfinite(centre).all():
        raise ValueError("centre_hz must be finite")
    spacing = bandwidth_hz / channel_count
    offsets = (
        np.arange(channel_count, dtype=float) - 0.5 * (channel_count - 1)
    ) * spacing
    return centre[..., None] + offsets


def band_averaged_point_visibility(
    uu_seconds: np.ndarray | float,
    vv_seconds: np.ndarray | float,
    l_radians: np.ndarray | float,
    m_radians: np.ndarray | float,
    centre_hz: np.ndarray | float,
    *,
    flux_jy: np.ndarray | float = 1.0,
    bandwidth_hz: float = 8.0e6,
    channel_count: int = 32,
) -> np.ndarray:
    """Visibility of a point source, integrated over an averaged IF.

    ``uu_seconds`` and ``vv_seconds`` are the UVFITS random parameters.  The
    phase convention is ``exp[-2 pi i (u l + v m)]``.
    """
    uu, vv, ll, mm, centre, flux = np.broadcast_arrays(
        np.asarray(uu_seconds, dtype=float),
        np.asarray(vv_seconds, dtype=float),
        np.asarray(l_radians, dtype=float),
        np.asarray(m_radians, dtype=float),
        np.asarray(centre_hz, dtype=float),
        np.asarray(flux_jy, dtype=float),
    )
    frequencies = subchannel_frequencies_hz(
        centre, bandwidth_hz=bandwidth_hz, channel_count=channel_count
    )
    delay_seconds = uu * ll + vv * mm
    phase = -2j * np.pi * delay_seconds[..., None] * frequencies
    return flux * np.mean(np.exp(phase), axis=-1)


def complex_gaussian_log_likelihood(
    observed: np.ndarray,
    model: np.ndarray,
    component_sigma: np.ndarray | float,
    *,
    include_normalization: bool = True,
) -> float:
    """Log likelihood for independent Gaussian real and imaginary parts."""
    observed, model, sigma = np.broadcast_arrays(
        np.asarray(observed, dtype=np.complex128),
        np.asarray(model, dtype=np.complex128),
        np.asarray(component_sigma, dtype=float),
    )
    if not np.isfinite(observed.real).all() or not np.isfinite(observed.imag).all():
        raise ValueError("observed visibilities must be finite")
    if not np.isfinite(model.real).all() or not np.isfinite(model.imag).all():
        raise ValueError("model visibilities must be finite")
    if not np.isfinite(sigma).all() or np.any(sigma <= 0):
        raise ValueError("component_sigma must be finite and positive")
    residual_power = np.abs(observed - model) ** 2
    log_like = -0.5 * np.sum(residual_power / sigma**2, dtype=np.float64)
    if include_normalization:
        log_like -= np.sum(np.log(2.0 * np.pi * sigma**2), dtype=np.float64)
    return float(log_like)


def classify_noise_bins(
    noise_bins: pd.DataFrame,
    *,
    robust_z_threshold: float = 6.0,
    minimum_difference_count: int = 64,
) -> pd.DataFrame:
    """Classify published baseline removals and data-driven noisy intervals.

    Exact time intervals for four RFI-sensitive baselines were not published.
    For those baselines only, high-noise intervals are identified from the
    median absolute deviation of log RMS within that baseline.  The full
    EF-JB-WB triangle is removed exactly as described by Powell et al. (2025).
    """
    required = {
        "baseline_code", "time_bin", "complex_difference_count",
        "adjacent_sigma_jy",
    }
    missing = required.difference(noise_bins.columns)
    if missing:
        raise ValueError(f"noise-bin table is missing columns: {sorted(missing)}")
    if robust_z_threshold <= 0 or minimum_difference_count <= 0:
        raise ValueError("flagging controls must be positive")

    frame = noise_bins.copy()
    sigma = frame["adjacent_sigma_jy"].to_numpy(dtype=float)
    frame["published_triangle_flag"] = frame["baseline_code"].isin(
        PUBLISHED_TRIANGLE_CODES
    )
    frame["insufficient_pairs_flag"] = (
        frame["complex_difference_count"].to_numpy(dtype=int)
        < minimum_difference_count
    ) | ~np.isfinite(sigma) | (sigma <= 0)
    frame["rfi_sensitive_baseline"] = frame["baseline_code"].isin(
        PUBLISHED_RFI_BASELINE_CODES
    )
    frame["log_sigma_robust_z"] = 0.0

    for baseline_code in PUBLISHED_RFI_BASELINE_CODES:
        selected = frame["baseline_code"].to_numpy(dtype=int) == baseline_code
        valid = selected & np.isfinite(sigma) & (sigma > 0)
        if not np.any(valid):
            continue
        log_sigma = np.log(sigma[valid])
        median = float(np.median(log_sigma))
        mad = float(np.median(np.abs(log_sigma - median)))
        scale = 1.4826 * mad
        robust_z = np.zeros(log_sigma.shape, dtype=float)
        if scale > 0:
            robust_z = (log_sigma - median) / scale
        frame.loc[valid, "log_sigma_robust_z"] = robust_z

    frame["data_driven_rfi_flag"] = (
        frame["rfi_sensitive_baseline"]
        & (frame["log_sigma_robust_z"] > robust_z_threshold)
    )
    frame["likelihood_excluded"] = (
        frame["published_triangle_flag"]
        | frame["insufficient_pairs_flag"]
        | frame["data_driven_rfi_flag"]
    )
    reason = np.full(len(frame), "retained", dtype=object)
    reason[frame["insufficient_pairs_flag"].to_numpy()] = "insufficient_pairs"
    reason[frame["data_driven_rfi_flag"].to_numpy()] = "data_driven_rfi"
    reason[frame["published_triangle_flag"].to_numpy()] = "published_triangle"
    frame["exclusion_reason"] = reason
    return frame


@dataclass(frozen=True)
class VisibilityNoiseModel:
    """Thirty-minute per-baseline noise lookup for complex visibilities."""

    start_jd: float
    bin_minutes: float
    sigma_by_key: dict[tuple[int, int], float]
    excluded_keys: frozenset[tuple[int, int]]
    baseline_fallback: dict[int, float]
    global_fallback: float

    @classmethod
    def from_frame(
        cls,
        noise_bins: pd.DataFrame,
        *,
        start_jd: float,
        bin_minutes: float = 30.0,
        apply_classification: bool = True,
    ) -> "VisibilityNoiseModel":
        if not np.isfinite(start_jd):
            raise ValueError("start_jd must be finite")
        if bin_minutes <= 0:
            raise ValueError("bin_minutes must be positive")
        frame = (
            classify_noise_bins(noise_bins)
            if apply_classification and "likelihood_excluded" not in noise_bins
            else noise_bins.copy()
        )
        if "likelihood_excluded" not in frame:
            frame["likelihood_excluded"] = False
        valid_sigma = (
            np.isfinite(frame["adjacent_sigma_jy"])
            & (frame["adjacent_sigma_jy"] > 0)
        )
        retained = frame[valid_sigma & ~frame["likelihood_excluded"]]
        if retained.empty:
            raise ValueError("noise-bin table contains no retained measurements")
        sigma_by_key = {
            (int(row.baseline_code), int(row.time_bin)): float(row.adjacent_sigma_jy)
            for row in retained.itertuples()
        }
        excluded_keys = frozenset(
            (int(row.baseline_code), int(row.time_bin))
            for row in frame[frame["likelihood_excluded"]].itertuples()
        )
        baseline_fallback = {
            int(code): float(values.median())
            for code, values in retained.groupby("baseline_code")["adjacent_sigma_jy"]
        }
        return cls(
            start_jd=float(start_jd),
            bin_minutes=float(bin_minutes),
            sigma_by_key=sigma_by_key,
            excluded_keys=excluded_keys,
            baseline_fallback=baseline_fallback,
            global_fallback=float(retained["adjacent_sigma_jy"].median()),
        )

    @classmethod
    def from_csv(
        cls,
        path: str | Path,
        *,
        start_jd: float,
        bin_minutes: float = 30.0,
        apply_classification: bool = True,
    ) -> "VisibilityNoiseModel":
        return cls.from_frame(
            pd.read_csv(path), start_jd=start_jd, bin_minutes=bin_minutes,
            apply_classification=apply_classification,
        )

    def lookup(
        self,
        baseline_codes: np.ndarray | int,
        times_jd: np.ndarray | float,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return component sigma and retained mask for each visibility."""
        baselines, times = np.broadcast_arrays(
            np.asarray(baseline_codes, dtype=int),
            np.asarray(times_jd, dtype=float),
        )
        if not np.isfinite(times).all():
            raise ValueError("times_jd must be finite")
        bins = np.floor(
            (times - self.start_jd) * 1440.0 / self.bin_minutes
        ).astype(int)
        sigma = np.empty(baselines.shape, dtype=float)
        retained = np.ones(baselines.shape, dtype=bool)
        for index in np.ndindex(baselines.shape):
            baseline = int(baselines[index])
            key = (baseline, int(bins[index]))
            sigma[index] = self.sigma_by_key.get(
                key, self.baseline_fallback.get(baseline, self.global_fallback)
            )
            retained[index] = key not in self.excluded_keys
            if baseline in PUBLISHED_TRIANGLE_CODES:
                retained[index] = False
        return sigma, retained
