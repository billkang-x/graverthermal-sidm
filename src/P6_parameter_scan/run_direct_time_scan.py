"""Run direct physical B1938 fits at several evolved source times.

This is the second stage of the P6 scan.  Fluid source runs are read from a
parameter-scan summary, and each requested source time is mapped to the
nearest saved snapshot.  The direct resimulation then evolves the reconstructed
physical halo with no post-hoc rescaling.  Results are written one row per
parameter point and source time, so time selection is explicit in later
likelihood or plotting code.
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
sys.path.insert(0, os.path.join(_ROOT, "src", "P5_resim"))

from resim_born_point import run_direct_resimulation  # noqa: E402


def _slug(value: float) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "", f"{value:g}".replace(".", "p"))


def _parse_times(value: str, allow_initial: bool = False) -> list[float]:
    times = []
    for token in value.split(","):
        try:
            time_gyr = float(token.strip())
        except ValueError as exc:
            raise ValueError(f"invalid source time: {token!r}") from exc
        if not np.isfinite(time_gyr) or time_gyr < 0:
            raise ValueError("source times must be finite and non-negative")
        if time_gyr == 0 and not allow_initial:
            raise ValueError("time zero requires --include-initial")
        times.append(time_gyr)
    if not times:
        raise ValueError("source time list cannot be empty")
    return sorted(set(times))


def _source_path(value: str) -> str:
    path = Path(value)
    if path.is_absolute():
        return str(path)
    return str((Path(_ROOT) / path).resolve())


def _nearest_snapshot(diagnostic_path: Path, requested_time_gyr: float,
                      include_initial: bool = False) -> tuple[int, float]:
    frame = pd.read_csv(diagnostic_path)
    lower_bound = 0.0 if include_initial else np.nextafter(0.0, 1.0)
    frame = frame[frame["snapshot_time_gyr"] >= lower_bound].copy()
    if frame.empty:
        raise ValueError("source diagnostics contain no evolved snapshots")
    times = frame.groupby("snapshot_idx")["snapshot_time_gyr"].first()
    snapshot_idx = int((times - requested_time_gyr).abs().idxmin())
    return snapshot_idx, float(times.loc[snapshot_idx])


def _read_rows(path: Path) -> dict[tuple[str, float], dict]:
    if not path.is_file():
        return {}
    frame = pd.read_csv(path)
    required = {"scan_point_id", "requested_source_time_gyr"}
    if not required.issubset(frame.columns):
        return {}
    return {
        (str(row["scan_point_id"]), float(row["requested_source_time_gyr"])): row.to_dict()
        for _, row in frame.iterrows()
    }


def _write_rows(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def run_scan(scan_summary: str, output_root: str, source_times_gyr: list[float],
             n_shells: int, t_epsilon: float, resume: bool,
             use_posthoc_best: bool = False,
             include_initial: bool = False,
             flag_dissipation: bool = True,
             diagnostic_filename: str = "b1938_snapshot_diagnostics.csv") -> pd.DataFrame:
    scan_frame = pd.read_csv(scan_summary)
    required = {"point_id", "model", "fluid_output_dir"}
    missing = required.difference(scan_frame.columns)
    if missing:
        raise ValueError(f"scan summary lacks columns: {', '.join(sorted(missing))}")
    output = Path(output_root)
    output.mkdir(parents=True, exist_ok=True)
    summary_path = output / "direct_time_scan_summary.csv"
    existing = _read_rows(summary_path)

    for _, scan_row in scan_frame.iterrows():
        point_id = str(scan_row["point_id"])
        source_run = _source_path(str(scan_row["fluid_output_dir"]))
        diagnostic_path = Path(source_run) / diagnostic_filename
        times = source_times_gyr
        if use_posthoc_best:
            best_time = scan_row.get("best_posthoc_time_gyr")
            if pd.isna(best_time):
                continue
            times = [float(best_time)]
        for requested_time in times:
            key = (point_id, float(requested_time))
            if resume and key in existing and existing[key].get("status") == "complete":
                continue
            row = {
                "scan_point_id": point_id,
                "model": str(scan_row["model"]),
                "requested_source_time_gyr": float(requested_time),
                "source_run_dir": source_run,
                "status": "pending",
            }
            try:
                snapshot_idx, selected_time = _nearest_snapshot(
                    diagnostic_path, requested_time, include_initial
                )
                direct_dir = output / point_id / (
                    f"t{_slug(requested_time)}_n{n_shells}"
                )
                summary_file = direct_dir / "summary.csv"
                if summary_file.is_file():
                    direct = pd.read_csv(summary_file).iloc[0].to_dict()
                else:
                    direct = run_direct_resimulation(
                        source_run, str(direct_dir), n_shells=n_shells,
                        snapshot_idx=snapshot_idx, t_epsilon=t_epsilon,
                        include_initial=include_initial,
                        flag_dissipation=flag_dissipation,
                        diagnostic_filename=diagnostic_filename,
                    )
                row.update({
                    "status": "complete",
                    "selected_snapshot_idx": snapshot_idx,
                    "selected_source_time_gyr": selected_time,
                    "direct_output_dir": direct.get("output_dir"),
                    "direct_mass_chi2": direct.get("mass_chi2_direct"),
                    "direct_mass_ratio": direct.get("mass_ratio_direct"),
                    "direct_ratio_pull": direct.get("mass_ratio_offset_sigma"),
                    "direct_inner_pull": direct.get("mass_residual_inner_sigma"),
                    "direct_outer_pull": direct.get("mass_residual_outer_sigma"),
                    "direct_t_evo_phys_gyr": direct.get("t_evo_phys_gyr"),
                    "direct_vmax_kms": direct.get("v_max_final_kms"),
                    "direct_max_cooling_code": direct.get("max_cooling_seen_code"),
                    "direct_steps": direct.get("steps"),
                    "flag_dissipation": direct.get("flag_dissipation", flag_dissipation),
                })
            except Exception as exc:  # noqa: BLE001 - retain failed points
                row["status"] = "failed"
                row["status_reason"] = f"{type(exc).__name__}: {exc}"
            existing[key] = row
            _write_rows(summary_path, list(existing.values()))
            print(point_id, f"t={requested_time:g}", row["status"])

    frame = pd.DataFrame(existing.values()).sort_values(
        ["scan_point_id", "requested_source_time_gyr"]
    )
    frame.to_csv(summary_path, index=False)
    return frame


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scan-summary", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--source-times-gyr", default="0.05,0.1,0.2,0.3,0.4,0.5")
    parser.add_argument("--n-shells", type=int, default=96)
    parser.add_argument("--t-epsilon", type=float, default=0.005)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--include-initial", action="store_true",
        help="allow t=0 and use the initial saved snapshot",
    )
    parser.add_argument(
        "--use-posthoc-best", action="store_true",
        help="use each source row's best_posthoc_time_gyr instead of a common time grid",
    )
    parser.add_argument(
        "--disable-cooling", action="store_true",
        help="run direct physical controls with the dissipative source disabled",
    )
    parser.add_argument(
        "--diagnostic-filename", default="b1938_snapshot_diagnostics.csv",
        help="candidate diagnostic CSV inside each source run directory",
    )
    args = parser.parse_args()
    if args.n_shells < 16 or args.t_epsilon <= 0:
        parser.error("n-shells must be >=16 and t-epsilon must be positive")
    times = _parse_times(args.source_times_gyr, args.include_initial)
    frame = run_scan(
        args.scan_summary, args.output_root, times,
        args.n_shells, args.t_epsilon, args.resume, args.use_posthoc_best,
        args.include_initial, not args.disable_cooling, args.diagnostic_filename,
    )
    print(f"wrote {len(frame)} rows to {args.output_root}\\direct_time_scan_summary.csv")
    print(frame["status"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
