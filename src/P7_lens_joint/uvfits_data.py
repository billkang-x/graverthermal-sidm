"""Read and audit the public GM068 pipeline-calibrated UVFITS data."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits


@dataclass(frozen=True)
class UVFitsMetadata:
    path: str
    object_name: str
    telescope: str
    date_obs: str
    group_count: int
    if_count: int
    polarization_count: int
    frequencies_hz: tuple[float, ...]
    stokes_codes: tuple[int, ...]
    antenna_names: dict[int, str]

    def to_dict(self) -> dict:
        return asdict(self)


def decode_aips_baseline(codes: np.ndarray | float) -> tuple[np.ndarray, np.ndarray]:
    """Decode the standard AIPS ``256 * antenna_1 + antenna_2`` value."""
    values = np.rint(np.asarray(codes, dtype=float)).astype(np.int64)
    if np.any(values <= 0):
        raise ValueError("baseline codes must be positive")
    antenna_1 = values // 256
    antenna_2 = values % 256
    if np.any(antenna_1 <= 0) or np.any(antenna_2 <= 0):
        raise ValueError("invalid standard AIPS baseline code")
    return antenna_1, antenna_2


def weighted_parallel_hand_average(
    real: np.ndarray,
    imag: np.ndarray,
    weight: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Combine RR and LL into an inverse-variance weighted Stokes-I estimate."""
    real = np.asarray(real, dtype=float)
    imag = np.asarray(imag, dtype=float)
    weight = np.asarray(weight, dtype=float)
    if real.shape != imag.shape or real.shape != weight.shape:
        raise ValueError("real, imag, and weight arrays must have identical shapes")
    if real.ndim < 1 or real.shape[-1] != 2:
        raise ValueError("the final axis must contain RR and LL")

    valid = np.isfinite(real) & np.isfinite(imag) & np.isfinite(weight) & (weight > 0)
    safe_weight = np.where(valid, weight, 0.0)
    total_weight = np.sum(safe_weight, axis=-1)
    numerator = np.sum(safe_weight * (real + 1j * imag), axis=-1)
    visibility = np.zeros(total_weight.shape, dtype=np.complex128)
    np.divide(numerator, total_weight, out=visibility, where=total_weight > 0)
    return visibility, total_weight


def read_uvfits_metadata(path: str | Path) -> UVFitsMetadata:
    path = Path(path)
    with fits.open(path, memmap=True) as hdus:
        header = hdus[0].header
        if_count = int(header["NAXIS5"])
        polarization_count = int(header["NAXIS3"])
        frequency_reference = float(header["CRVAL4"])
        frequency_offsets = np.asarray(hdus["AIPS FQ"].data[0]["IF FREQ"], dtype=float)
        if frequency_offsets.size != if_count:
            raise ValueError("AIPS FQ table does not match the IF axis")
        frequencies = frequency_reference + frequency_offsets
        stokes = (
            float(header["CRVAL3"])
            + (np.arange(polarization_count) + 1.0 - float(header["CRPIX3"]))
            * float(header["CDELT3"])
        )
        antennas = {
            int(row["NOSTA"]): str(row["ANNAME"]).strip()
            for row in hdus["AIPS AN"].data
        }
        return UVFitsMetadata(
            path=str(path.resolve()),
            object_name=str(header.get("OBJECT", "")).strip(),
            telescope=str(header.get("TELESCOP", "")).strip(),
            date_obs=str(header.get("DATE-OBS", "")).strip(),
            group_count=int(header["GCOUNT"]),
            if_count=if_count,
            polarization_count=polarization_count,
            frequencies_hz=tuple(float(value) for value in frequencies),
            stokes_codes=tuple(int(round(value)) for value in stokes),
            antenna_names=antennas,
        )


def _format_baseline(code: int, antennas: dict[int, str]) -> tuple[int, int, str]:
    antenna_1, antenna_2 = decode_aips_baseline(float(code))
    first = int(antenna_1)
    second = int(antenna_2)
    return first, second, f"{antennas.get(first, first)}-{antennas.get(second, second)}"


def audit_visibility_noise(
    path: str | Path,
    *,
    chunk_groups: int = 250_000,
    bin_minutes: float = 30.0,
    max_pair_separation_seconds: float = 10.0,
    uv_sample_stride: int = 500,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    """Audit weights and adjacent-visibility noise without loading the full file."""
    if chunk_groups <= 0 or bin_minutes <= 0 or max_pair_separation_seconds <= 0:
        raise ValueError("chunk and time controls must be positive")
    if uv_sample_stride <= 0:
        raise ValueError("uv_sample_stride must be positive")

    metadata = read_uvfits_metadata(path)
    frequencies = np.asarray(metadata.frequencies_hz, dtype=float)
    noise_accumulator: dict[tuple[int, int], list[float]] = {}
    baseline_accumulator: dict[int, list[float]] = {}
    previous: dict[int, tuple[float, np.ndarray, np.ndarray]] = {}
    uv_rows: list[pd.DataFrame] = []
    total_valid = 0
    total_samples = 0
    start_jd: float | None = None
    end_jd: float | None = None

    with fits.open(path, memmap=True) as hdus:
        groups = hdus[0].data
        group_count = metadata.group_count
        for start in range(0, group_count, chunk_groups):
            stop = min(start + chunk_groups, group_count)
            times = np.asarray(groups.par("DATE")[start:stop], dtype=float)
            baselines = np.rint(groups.par("BASELINE")[start:stop]).astype(np.int64)
            uu_seconds = np.asarray(groups.par("UU---SIN")[start:stop], dtype=float)
            vv_seconds = np.asarray(groups.par("VV---SIN")[start:stop], dtype=float)
            raw = np.asarray(groups.data[start:stop])[:, 0, 0, :, 0, :, :]
            real = raw[..., 0].astype(float, copy=False)
            imag = raw[..., 1].astype(float, copy=False)
            weights = raw[..., 2].astype(float, copy=False)
            valid = (
                np.isfinite(real) & np.isfinite(imag) & np.isfinite(weights)
                & (weights > 0)
            )
            total_valid += int(np.count_nonzero(valid))
            total_samples += int(valid.size)
            if start_jd is None:
                finite_times = times[np.isfinite(times)]
                if finite_times.size == 0:
                    raise ValueError("UVFITS contains no finite group dates")
                start_jd = float(np.min(finite_times))
            finite_times = times[np.isfinite(times)]
            if finite_times.size:
                chunk_end_jd = float(np.max(finite_times))
                end_jd = chunk_end_jd if end_jd is None else max(end_jd, chunk_end_jd)

            for baseline in np.unique(baselines):
                indices = np.flatnonzero(baselines == baseline)
                current_times = times[indices]
                current_vis = real[indices] + 1j * imag[indices]
                current_weights = weights[indices]
                current_valid = valid[indices]
                if baseline in previous:
                    prior_time, prior_vis, prior_weights = previous[baseline]
                    current_times = np.concatenate(([prior_time], current_times))
                    current_vis = np.concatenate((prior_vis[None, ...], current_vis), axis=0)
                    current_weights = np.concatenate(
                        (prior_weights[None, ...], current_weights), axis=0
                    )
                    prior_valid = (
                        np.isfinite(prior_vis.real) & np.isfinite(prior_vis.imag)
                        & np.isfinite(prior_weights) & (prior_weights > 0)
                    )
                    current_valid = np.concatenate(
                        (prior_valid[None, ...], current_valid), axis=0
                    )

                previous[int(baseline)] = (
                    float(current_times[-1]),
                    np.asarray(current_vis[-1]),
                    np.asarray(current_weights[-1]),
                )
                if current_times.size < 2:
                    continue

                separation_seconds = np.diff(current_times) * 86400.0
                pair_valid = (
                    current_valid[1:] & current_valid[:-1]
                    & (separation_seconds[:, None, None] > 0)
                    & (separation_seconds[:, None, None] <= max_pair_separation_seconds)
                )
                differences = (current_vis[1:] - current_vis[:-1]) / np.sqrt(2.0)
                time_bins = np.floor(
                    (current_times[1:] - start_jd) * 1440.0 / bin_minutes
                ).astype(int)
                expanded_bins = np.broadcast_to(
                    time_bins[:, None, None], pair_valid.shape
                )
                selected_bins = expanded_bins[pair_valid]
                selected_power = np.abs(differences[pair_valid]) ** 2
                if selected_bins.size:
                    counts = np.bincount(selected_bins)
                    powers = np.bincount(selected_bins, weights=selected_power)
                    for bin_index in np.flatnonzero(counts):
                        key = (int(baseline), int(bin_index))
                        target = noise_accumulator.setdefault(key, [0.0, 0.0])
                        target[0] += float(counts[bin_index])
                        target[1] += float(powers[bin_index])

                positive_weights = current_weights[1:][current_valid[1:]]
                target = baseline_accumulator.setdefault(int(baseline), [0.0, 0.0, 0.0])
                target[0] += float(positive_weights.size)
                target[1] += float(np.sum(positive_weights))
                target[2] += float(np.sum(1.0 / positive_weights))

            sample_indices = np.arange(start, stop, uv_sample_stride, dtype=int) - start
            if sample_indices.size:
                sampled_baselines = baselines[sample_indices]
                sampled_u = uu_seconds[sample_indices, None] * frequencies[None, :] / 1e6
                sampled_v = vv_seconds[sample_indices, None] * frequencies[None, :] / 1e6
                sampled_real = real[sample_indices]
                sampled_imag = imag[sample_indices]
                sampled_weight = weights[sample_indices]
                stokes_i, stokes_weight = weighted_parallel_hand_average(
                    sampled_real, sampled_imag, sampled_weight
                )
                uv_rows.append(pd.DataFrame({
                    "group_index": np.repeat(
                        np.arange(start, stop, uv_sample_stride, dtype=int), frequencies.size
                    ),
                    "baseline_code": np.repeat(sampled_baselines, frequencies.size),
                    "if_index": np.tile(np.arange(frequencies.size), sample_indices.size),
                    "u_mlambda": sampled_u.ravel(),
                    "v_mlambda": sampled_v.ravel(),
                    "stokes_i_real_jy": stokes_i.ravel().real,
                    "stokes_i_imag_jy": stokes_i.ravel().imag,
                    "stokes_i_weight": stokes_weight.ravel(),
                }))

    noise_rows = []
    for (baseline, bin_index), (count, power) in sorted(noise_accumulator.items()):
        first, second, label = _format_baseline(baseline, metadata.antenna_names)
        noise_rows.append({
            "baseline_code": baseline,
            "antenna_1": first,
            "antenna_2": second,
            "baseline": label,
            "time_bin": bin_index,
            "start_hour": bin_index * bin_minutes / 60.0,
            "end_hour": (bin_index + 1) * bin_minutes / 60.0,
            "complex_difference_count": int(count),
            "adjacent_sigma_jy": float(np.sqrt(power / (2.0 * count))),
        })
    noise_frame = pd.DataFrame(noise_rows)

    baseline_rows = []
    for baseline, (count, weight_sum, inverse_weight_sum) in sorted(
        baseline_accumulator.items()
    ):
        first, second, label = _format_baseline(baseline, metadata.antenna_names)
        subset = noise_frame[noise_frame["baseline_code"] == baseline]
        baseline_rows.append({
            "baseline_code": baseline,
            "antenna_1": first,
            "antenna_2": second,
            "baseline": label,
            "valid_visibility_count": int(count),
            "mean_pipeline_weight": weight_sum / count if count else np.nan,
            "rms_pipeline_sigma_jy": (
                np.sqrt(inverse_weight_sum / count) if count else np.nan
            ),
            "median_adjacent_sigma_jy": (
                float(subset["adjacent_sigma_jy"].median()) if not subset.empty else np.nan
            ),
            "noise_bin_count": int(len(subset)),
        })
    baseline_frame = pd.DataFrame(baseline_rows)
    uv_frame = pd.concat(uv_rows, ignore_index=True) if uv_rows else pd.DataFrame()
    summary = {
        **metadata.to_dict(),
        "start_jd": start_jd,
        "end_jd": end_jd,
        "duration_hours": (
            (end_jd - start_jd) * 24.0
            if start_jd is not None and end_jd is not None else None
        ),
        "total_complex_visibility_count": total_samples,
        "valid_complex_visibility_count": total_valid,
        "valid_fraction": total_valid / total_samples,
        "noise_bin_minutes": bin_minutes,
        "max_pair_separation_seconds": max_pair_separation_seconds,
        "uv_sample_stride_groups": uv_sample_stride,
    }
    return noise_frame, baseline_frame, uv_frame, summary
