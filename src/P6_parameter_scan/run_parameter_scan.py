"""Run a resumable Born-valid fluid and direct B1938 parameter scan.

The scan is deliberately two-stage.  Each point is constructed analytically
at fixed ``eta_B``, threshold velocity, and ``sigma/m`` at 100 km/s, then
filtered by the requested Born mask before any fluid evolution is attempted.
Passing points use the existing long fluid runner.  With ``--direct`` the best
evolved diagnostic candidate is then resimulated in physical units, without
post-hoc rescaling, using :mod:`P5_resim.resim_born_point`.

The output CSV is an audit table: post-hoc diagnostic columns are explicitly
named as such, while direct columns come only from physical resimulation.
Incomplete points are retained with ``status=failed`` so an interrupted scan
can be inspected and resumed safely.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "..", ".."))
sys.path[:0] = [
    os.path.join(_ROOT, "src", "cross_sections"),
    os.path.join(_ROOT, "src", "fluid_runner"),
    os.path.join(_ROOT, "src", "P3_rescaling"),
    os.path.join(_ROOT, "src", "P5_resim"),
]

from born_valid_scan import (  # noqa: E402
    parameters_for_target_sigma,
    parameters_for_target_sigma_at_threshold,
    threshold_velocity_km_s,
)
from dsidm_models import (  # noqa: E402
    benchmark_models,
    born_expansion_parameter,
    sigma_T_born,
)
from run_born_long import run_one_long  # noqa: E402
from resim_born_point import run_direct_resimulation  # noqa: E402


DEFAULT_MODELS = ("M1_dark_photon_massive", "M2_scalar_phi_massive")


def _parse_float_list(value: str, *, allow_none: bool = False) -> list[float | None]:
    """Parse a comma-separated positive float list."""
    result: list[float | None] = []
    for token in value.split(","):
        token = token.strip().lower()
        if allow_none and token in {"none", "benchmark", "default"}:
            result.append(None)
            continue
        try:
            number = float(token)
        except ValueError as exc:
            raise ValueError(f"invalid numeric scan value: {token!r}") from exc
        if not np.isfinite(number) or number <= 0:
            raise ValueError(f"scan values must be finite and positive: {token!r}")
        result.append(number)
    if not result:
        raise ValueError("scan value list cannot be empty")
    return result


def _slug(value: float | None) -> str:
    if value is None:
        return "benchmark"
    text = f"{value:g}"
    text = text.replace(".", "p").replace("-", "m")
    return re.sub(r"[^A-Za-z0-9_]+", "", text)


def _model_point(model_key: str, eta: float, sigma_m_100: float,
                 threshold_requested: float | None, born_max: float) -> dict:
    """Construct one physical point and evaluate its analytic diagnostics."""
    models = benchmark_models()
    if model_key not in DEFAULT_MODELS:
        raise ValueError(f"unsupported fluid-scan model: {model_key}")
    base = models[model_key]
    if threshold_requested is None:
        model = parameters_for_target_sigma(base, eta, sigma_m_100, 100.0)
    else:
        model = parameters_for_target_sigma_at_threshold(
            base, eta, sigma_m_100, threshold_requested, 100.0
        )
    eta_actual = float(born_expansion_parameter(model))
    sigma_actual = float(sigma_T_born(np.array([100.0]), model)[0])
    return {
        "model": model_key,
        "eta_target": float(eta),
        "eta_B": eta_actual,
        "born_valid": bool(eta_actual < born_max),
        "born_controlled": bool(eta_actual <= 0.1 + 1e-12),
        "target_sigma_m_100_kms": float(sigma_m_100),
        "sigma_m_100_kms": sigma_actual,
        "threshold_requested_kms": (
            None if threshold_requested is None else float(threshold_requested)
        ),
        "threshold_actual_kms": float(threshold_velocity_km_s(model)),
        "m_chi_GeV": float(model.m_chi),
        "m_mediator_keV": float(model.m_mediator * 1e6),
        "alpha_D": float(model.alpha_D),
    }


def _read_existing(path: Path) -> dict[str, dict]:
    if not path.is_file():
        return {}
    frame = pd.read_csv(path)
    if "point_id" not in frame:
        return {}
    return {str(row["point_id"]): row.to_dict() for _, row in frame.iterrows()}


def _write_rows(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _posthoc_best(diagnostic_path: str) -> dict:
    frame = pd.read_csv(diagnostic_path)
    frame = frame[frame["snapshot_time_gyr"] > 0]
    if frame.empty:
        return {}
    best = frame.loc[frame["mass_chi2"].idxmin()]
    return {
        "best_posthoc_mass_chi2": float(best["mass_chi2"]),
        "best_posthoc_time_gyr": float(best["snapshot_time_gyr"]),
        "best_posthoc_mass_ratio": float(best["mass_ratio"]),
        "best_posthoc_r2d_rs": float(best["r2D_rs"]),
    }


def run_point(point: dict, output_root: Path, *, n_shells: int,
              t_end_gyr: float, checkpoint_gyr: float, t_epsilon: float,
              direct: bool, direct_n_shells: int, max_steps: int | None) -> dict:
    """Run one point and return a single audit row."""
    model_key = str(point["model"])
    threshold = point["threshold_requested_kms"]
    point_id = (
        f"{model_key.split('_')[0]}_eta{_slug(point['eta_target'])}"
        f"_sigma{_slug(point['target_sigma_m_100_kms'])}"
        f"_vstar{_slug(threshold)}"
    )
    run_root = output_root / "runs" / point_id
    row = {"point_id": point_id, **point, "output_root": str(run_root)}
    row["status"] = "filtered" if not point["born_valid"] else "pending"
    if not point["born_valid"]:
        row["status_reason"] = "eta_B >= born_max"
        return row

    try:
        fluid = run_one_long(
            model_key, float(point["eta_target"]), str(run_root), n_shells,
            t_end_gyr, checkpoint_gyr, float(point["target_sigma_m_100_kms"]),
            True, max_steps, threshold, t_epsilon,
        )
        row.update({
            "status": "fluid_complete",
            "fluid_output_dir": fluid["output_dir"],
            "fluid_current_time_gyr": fluid["current_time_gyr"],
            "fluid_v_max_final_kms": fluid["v_max_final_kms"],
            "fluid_v_max_to_threshold": fluid["v_max_to_threshold_final"],
            "fluid_max_cooling_seen_code": fluid["max_cooling_seen_code"],
            "diagnostic_path": fluid["diagnostic_path"],
        })
        row.update(_posthoc_best(fluid["diagnostic_path"]))
        if direct:
            direct_dir = run_root / "direct" / (
                f"n{direct_n_shells}_t{_slug(t_end_gyr)}"
            )
            direct_summary = direct_dir / "summary.csv"
            if direct_summary.is_file():
                direct_row = pd.read_csv(direct_summary).iloc[0].to_dict()
            else:
                direct_row = run_direct_resimulation(
                    fluid["output_dir"], str(direct_dir),
                    n_shells=direct_n_shells, t_epsilon=t_epsilon,
                )
            for key, value in direct_row.items():
                if key == "output_dir":
                    row["direct_output_dir"] = value
                elif key != "model":
                    row[f"direct_{key}"] = value
            row["status"] = "direct_complete"
    except Exception as exc:  # noqa: BLE001 - retain failed points for audit
        row["status"] = "failed"
        row["status_reason"] = f"{type(exc).__name__}: {exc}"
    return row


def run_scan(args: argparse.Namespace) -> pd.DataFrame:
    output_root = Path(args.output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    eta_values = _parse_float_list(args.eta_values)
    sigma_values = _parse_float_list(args.sigma_values)
    threshold_values = _parse_float_list(args.threshold_values, allow_none=True)
    models = [value.strip() for value in args.models.split(",") if value.strip()]
    with (output_root / "scan_config.json").open("w", encoding="utf-8") as handle:
        json.dump(vars(args), handle, indent=2, sort_keys=True)
        handle.write("\n")
    existing = _read_existing(output_root / "parameter_scan_summary.csv")
    rows = dict(existing)
    for model_key in models:
        for eta in eta_values:
            for sigma_m_100 in sigma_values:
                for threshold in threshold_values:
                    point = _model_point(
                        model_key, float(eta), float(sigma_m_100), threshold,
                        args.born_max,
                    )
                    point_id = (
                        f"{model_key.split('_')[0]}_eta{_slug(eta)}"
                        f"_sigma{_slug(sigma_m_100)}_vstar{_slug(threshold)}"
                    )
                    if args.resume and point_id in rows:
                        old = rows[point_id]
                        old_time = float(old.get("fluid_current_time_gyr", 0.0) or 0.0)
                        complete = old.get("status") in {
                            "filtered", "fluid_complete", "direct_complete"
                        }
                        if complete and old_time >= float(args.t_end_gyr) * (1.0 - 1e-8):
                            continue
                    if not point["born_valid"]:
                        point["point_id"] = point_id
                        point["status"] = "filtered"
                        point["status_reason"] = "eta_B >= born_max"
                        rows[point_id] = point
                    else:
                        result = run_point(
                            point, output_root, n_shells=args.n_shells,
                            t_end_gyr=args.t_end_gyr,
                            checkpoint_gyr=args.checkpoint_gyr,
                            t_epsilon=args.t_epsilon,
                            direct=args.direct,
                            direct_n_shells=args.direct_n_shells,
                            max_steps=args.max_steps,
                        )
                        rows[point_id] = result
                    frame = pd.DataFrame(rows.values()).sort_values("point_id")
                    frame.to_csv(output_root / "parameter_scan_summary.csv", index=False)
                    print(point_id, rows[point_id]["status"])
    frame = pd.DataFrame(rows.values()).sort_values("point_id")
    frame.to_csv(output_root / "parameter_scan_summary.csv", index=False)
    return frame


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", default=os.path.join(_ROOT, "data", "P6_parameter_scan"))
    parser.add_argument("--models", default=",".join(DEFAULT_MODELS))
    parser.add_argument("--eta-values", default="0.1")
    parser.add_argument("--sigma-values", default="0.1,1,10")
    parser.add_argument("--threshold-values", default="20,100,500",
                        help="comma-separated threshold velocities; include 'none' for each model benchmark")
    parser.add_argument("--born-max", type=float, default=1.0,
                        help="strict eta_B upper bound for fluid evolution")
    parser.add_argument("--t-end-gyr", type=float, default=0.02)
    parser.add_argument("--checkpoint-gyr", type=float, default=0.005)
    parser.add_argument("--n-shells", type=int, default=48)
    parser.add_argument("--direct-n-shells", type=int, default=96)
    parser.add_argument("--t-epsilon", type=float, default=1e-2)
    parser.add_argument("--max-steps", type=int, default=None)
    parser.add_argument("--direct", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if args.t_end_gyr <= 0 or args.checkpoint_gyr <= 0:
        parser.error("time settings must be positive")
    if args.n_shells < 16 or args.direct_n_shells < 16:
        parser.error("shell counts must be at least 16")
    if not np.isfinite(args.born_max) or args.born_max <= 0:
        parser.error("--born-max must be finite and positive")
    frame = run_scan(args)
    print(f"wrote {len(frame)} rows to {args.output_root}\\parameter_scan_summary.csv")
    print(frame[["point_id", "status"]].to_string(index=False))


if __name__ == "__main__":
    main()
