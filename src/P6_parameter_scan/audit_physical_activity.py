"""Audit threshold reachability and cooling times in direct physical runs.

The source-fluid velocity is not the velocity of the directly resimulated
physical halo after the elastic rescaling.  This audit therefore restores the
physical direct state and evaluates the threshold ratio and local cooling time
in the same frame used for the B1938 mass fit.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(_ROOT, "src", "fluid_runner"))

from run_born_valid import build_halo  # noqa: E402
from dsidm_models import born_expansion_parameter  # noqa: E402


_POINT_RE = re.compile(
    r"_eta(?P<eta>[0-9p]+)_sigma(?P<sigma>[0-9p]+)_vstar(?P<vstar>[0-9p]+)$"
)


def _decode_number(value: str) -> float:
    return float(value.replace("p", "."))


def parse_point_id(point_id: str) -> tuple[float, float, float]:
    """Extract eta, target sigma/m, and threshold from a scan point id."""
    match = _POINT_RE.search(str(point_id))
    if match is None:
        raise ValueError(f"unrecognized final-candidate point id: {point_id}")
    return tuple(_decode_number(match.group(name)) for name in ("eta", "sigma", "vstar"))


def _finite_max(values: pd.Series) -> float:
    arr = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    arr = arr[np.isfinite(arr)]
    return float(np.max(arr)) if arr.size else float("nan")


def audit_scan(scan_summary: str, output_csv: str, n_shells: int = 192) -> pd.DataFrame:
    """Audit one unified direct-time scan and write one row per point."""
    frame = pd.read_csv(scan_summary)
    required = {
        "scan_point_id", "model", "requested_source_time_gyr", "status",
        "direct_output_dir", "direct_vmax_kms", "direct_max_cooling_code",
        "direct_t_evo_phys_gyr",
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"scan summary lacks columns: {', '.join(sorted(missing))}")

    rows: list[dict] = []
    for point_id, group in frame.groupby("scan_point_id", sort=True):
        if not (group["status"] == "complete").all():
            rows.append({
                "scan_point_id": point_id,
                "model": str(group.iloc[0]["model"]),
                "status": "incomplete_scan",
            })
            continue

        eta, sigma_m_100, threshold = parse_point_id(str(point_id))
        initial = group.loc[
            group["requested_source_time_gyr"].astype(float).abs() < 1e-15
        ]
        if initial.empty:
            initial = group.sort_values("requested_source_time_gyr").iloc[[0]]
        direct_dir = str(initial.iloc[0]["direct_output_dir"])
        if not Path(direct_dir).is_dir():
            raise FileNotFoundError(f"direct output directory not found: {direct_dir}")

        # build_halo restores the latest saved state from this direct run and
        # recomputes the exact microphysical cooling array in physical units.
        model, halo, sigma_reference = build_halo(
            str(group.iloc[0]["model"]), eta, direct_dir, n_shells,
            sigma_m_100, True, threshold, 0.0, 0.0,
        )
        v_max_phys = _finite_max(group["direct_vmax_kms"])
        max_cooling_code = _finite_max(group["direct_max_cooling_code"])
        max_t_evo_phys = _finite_max(group["direct_t_evo_phys_gyr"])

        u = np.asarray(halo.u, dtype=float)
        cooling = np.asarray(getattr(halo, "C_cool", np.zeros_like(u)), dtype=float)
        finite_positive = np.isfinite(u) & np.isfinite(cooling) & (u > 0) & (cooling > 0)
        if np.any(finite_positive):
            min_cooling_time_code = float(np.min(u[finite_positive] / cooling[finite_positive]))
            scale_t_gyr = float(halo.scale_t.to("Gyr").value)
            min_cooling_time_gyr = min_cooling_time_code * scale_t_gyr
            max_fractional_cooling_rate_gyr = 1.0 / min_cooling_time_gyr
        else:
            min_cooling_time_code = float("inf")
            min_cooling_time_gyr = float("inf")
            max_fractional_cooling_rate_gyr = 0.0

        threshold_ratio = v_max_phys / threshold if threshold > 0 else float("nan")
        if threshold_ratio >= 1.0 and min_cooling_time_gyr <= 1.0:
            activity_class = "threshold_reachable_and_fast_cooling"
        elif threshold_ratio >= 1.0:
            activity_class = "threshold_reachable_but_slow_cooling"
        elif min_cooling_time_gyr <= 1.0:
            activity_class = "cooling_fast_below_nominal_threshold"
        else:
            activity_class = "threshold_unreachable_and_cooling_slow"

        rows.append({
            "scan_point_id": str(point_id),
            "model": str(group.iloc[0]["model"]),
            "eta_B": eta,
            "sigma_m_100_kms": float(sigma_reference),
            "vstar_kms": threshold,
            "n_time_rows": int(len(group)),
            "max_direct_t_evo_phys_gyr": max_t_evo_phys,
            "vmax_phys_kms_max_over_scan": v_max_phys,
            "threshold_ratio_vmax_over_vstar": threshold_ratio,
            "max_cooling_code_over_scan": max_cooling_code,
            "min_cooling_time_code_initial_state": min_cooling_time_code,
            "min_cooling_time_gyr_initial_state": min_cooling_time_gyr,
            "max_fractional_cooling_rate_per_gyr_initial_state": max_fractional_cooling_rate_gyr,
            "born_eta": float(born_expansion_parameter(model)),
            "m_chi_GeV": float(model.m_chi),
            "m_mediator_keV": float(model.m_mediator * 1e6),
            "activity_class": activity_class,
            "status": "complete",
        })

    result = pd.DataFrame(rows).sort_values("scan_point_id")
    path = Path(output_csv)
    path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(path, index=False, quoting=csv.QUOTE_MINIMAL)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scan-summary", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--n-shells", type=int, default=192)
    args = parser.parse_args()
    if args.n_shells < 16:
        parser.error("n-shells must be at least 16")
    frame = audit_scan(args.scan_summary, args.output, args.n_shells)
    print(f"wrote {len(frame)} rows to {args.output}")
    print(frame["activity_class"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
