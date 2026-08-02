"""Compare long fluid runs at common diagnostic levels."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

import pandas as pd


def _model_directories(root: str) -> list[Path]:
    return sorted(
        path for path in Path(root).iterdir()
        if path.is_dir() and (path / "fluid_trajectory.csv").is_file()
    )


def _summary_by_model(root: str) -> dict[str, dict]:
    path = Path(root) / "summary.csv"
    if not path.is_file():
        return {}
    frame = pd.read_csv(path)
    result = {}
    for _, row in frame.iterrows():
        value = row.to_dict()
        result[str(row["model"])] = value
        output_dir = str(row.get("output_dir", ""))
        if output_dir:
            result[Path(output_dir).name] = value
    return result


def _config_by_directory(model_dir: Path) -> dict:
    """Use the per-run config when an older root summary is incomplete."""
    path = model_dir / "run_config.json"
    if not path.is_file():
        return {}
    try:
        with path.open(encoding="utf-8") as handle:
            return json.load(handle)
    except (OSError, json.JSONDecodeError):
        return {}


def compare_run(label: str, root: str) -> list[dict]:
    """Extract final trajectory and best B1938 diagnostics for one run root."""
    summary = _summary_by_model(root)
    rows = []
    for model_dir in _model_directories(root):
        config = summary.get(model_dir.name, {}) or _config_by_directory(model_dir)
        model = str(config.get("model", model_dir.name))
        trajectory = pd.read_csv(model_dir / "fluid_trajectory.csv")
        initial = trajectory.iloc[0]
        final = trajectory.iloc[-1]
        diag_path = model_dir / "b1938_snapshot_diagnostics.csv"
        diagnostics = pd.read_csv(diag_path) if diag_path.is_file() else pd.DataFrame()
        best = diagnostics.loc[diagnostics["mass_chi2"].idxmin()] if len(diagnostics) else {}
        row = {
            "run_label": label,
            "run_root": os.path.normpath(root),
            "model": model,
            "n_shells": config.get("n_shells"),
            "t_epsilon": config.get("t_epsilon"),
            "target_time_gyr": config.get("target_time_gyr"),
            "final_time_gyr": final["snapshot_time_gyr"],
            "rho_r0p1rs_initial_code": initial["rho_r0p1rs_code"],
            "rho_r0p1rs_final_code": final["rho_r0p1rs_code"],
            "rho_r0p1rs_fractional_change": (
                final["rho_r0p1rs_code"] / initial["rho_r0p1rs_code"] - 1.0
            ),
            "u_r0p1rs_fractional_change": (
                final["u_r0p1rs_code"] / initial["u_r0p1rs_code"] - 1.0
            ),
            "v_max_final_kms": final["v_max_kms"],
            "v_threshold_kms": config.get("v_threshold_kms"),
            "max_cooling_seen_code": config.get("max_cooling_seen_code"),
            "best_mass_chi2": best.get("mass_chi2"),
            "best_mass_fit_time_gyr": best.get("snapshot_time_gyr"),
            "best_mass_ratio": best.get("mass_ratio"),
        }
        rows.append(row)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run", action="append", required=True, metavar="LABEL=ROOT",
        help="run root with a display label, e.g. baseline=data/P2_born_long...",
    )
    parser.add_argument("--output", default="data/P2_fluid_convergence_summary.csv")
    args = parser.parse_args()

    rows = []
    for spec in args.run:
        if "=" not in spec:
            parser.error("--run must use LABEL=ROOT")
        label, root = spec.split("=", 1)
        if not Path(root).is_dir():
            parser.error(f"run root does not exist: {root}")
        rows.extend(compare_run(label, root))
    if not rows:
        parser.error("no completed model directories found")

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    frame = pd.DataFrame(rows)
    print(f"wrote {len(frame)} rows to {output}")
    print(frame.to_string(index=False))


if __name__ == "__main__":
    main()
