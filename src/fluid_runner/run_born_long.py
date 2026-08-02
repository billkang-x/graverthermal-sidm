"""Run resumable Born-valid long-time fluid experiments.

The short regression runner is intentionally limited to a few steps.  This
driver keeps the same microphysical setup but evolves to a requested physical
time, saving explicit checkpoints and refusing to silently change the model
configuration of an existing run.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from pathlib import Path

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "..", ".."))
sys.path[:0] = [
    os.path.join(_ROOT, "src", "cross_sections"),
    _HERE,
    os.path.join(_ROOT, "src", "P3_rescaling"),
    os.path.join(_ROOT, "external", "gravothermalsidm"),
]

from born_valid_scan import (
    parameters_at_eta,
    parameters_at_mass_threshold,
    parameters_for_target_sigma,
    parameters_for_target_sigma_at_threshold,
    threshold_velocity_km_s,
)
from dsidm_models import benchmark_models, born_expansion_parameter, sigma_T_born
from run_born_valid import MODEL_KEYS, build_halo
from snapshot_diagnostics import (
    diagnose_directory,
    fluid_trajectory,
    write_diagnostic_csv,
    write_trajectory_csv,
)


def _snapshot_files(output_dir: str) -> list[str]:
    from SourcePy.record import HaloRecord

    record = HaloRecord(output_dir)
    files, _ = record.glob_pickle_files()
    return list(files)


def _run_config(model_key: str, model, eta: float, target_sigma_100,
                n_shells: int, flag_dissipation: bool,
                threshold_velocity_requested=None) -> dict:
    return {
        "model": model_key,
        "eta_B": float(born_expansion_parameter(model)),
        "alpha_D": float(model.alpha_D),
        "m_chi_GeV": float(model.m_chi),
        "m_mediator_keV": float(model.m_mediator * 1e6),
        "eta_requested": float(eta),
        "target_sigma_m_100_kms": (
            None if target_sigma_100 is None else float(target_sigma_100)
        ),
        "threshold_velocity_requested_kms": (
            None if threshold_velocity_requested is None
            else float(threshold_velocity_requested)
        ),
        "n_shells": int(n_shells),
        "flag_dissipation": bool(flag_dissipation),
    }


def _load_or_write_config(output_dir: str, config: dict) -> None:
    path = os.path.join(output_dir, "run_config.json")
    if os.path.isfile(path):
        with open(path, encoding="utf-8") as handle:
            old = json.load(handle)
        keys = tuple(config)
        mismatches = [
            key for key in keys
            if key in old and old.get(key) != config.get(key)
        ]
        if mismatches:
            raise ValueError(
                "existing run configuration differs for: " + ", ".join(mismatches)
            )
        return
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(config, handle, indent=2, sort_keys=True)
        handle.write("\n")


def _previous_summary(output_root: str, output_dir: str) -> dict:
    """Load the prior model summary when a run is being resumed.

    Fluid snapshots do not serialize the cooling diagnostic itself, so a
    resumed run must carry historical summary values forward explicitly.
    """
    path = os.path.join(output_root, "summary.csv")
    if not os.path.isfile(path):
        return {}
    try:
        with open(path, newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                old_dir = row.get("output_dir", "")
                if os.path.normcase(os.path.abspath(old_dir)) == os.path.normcase(os.path.abspath(output_dir)):
                    return row
    except (OSError, csv.Error):
        return {}
    return {}


def _write_summary(output_root: str, rows: list[dict]) -> str:
    """Merge current model rows with other configurations in the same root."""
    path = os.path.join(output_root, "summary.csv")
    existing = []
    if os.path.isfile(path):
        try:
            with open(path, newline="", encoding="utf-8") as handle:
                existing = list(csv.DictReader(handle))
        except (OSError, csv.Error):
            existing = []
    merged = {}
    for row in existing:
        key = os.path.normcase(os.path.abspath(row.get("output_dir", "")))
        if key:
            merged[key] = row
    for row in rows:
        key = os.path.normcase(os.path.abspath(row.get("output_dir", "")))
        merged[key] = row
    all_rows = list(merged.values())
    fields = list(rows[0])
    for row in all_rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_rows)
    return path


def _assert_finite_state(halo, model_key: str) -> None:
    arrays = (halo.r, halo.m, halo.rho, halo.p, halo.u)
    if not all(np.all(np.isfinite(array)) for array in arrays):
        raise FloatingPointError(f"non-finite fluid state for {model_key}")
    if not (
        np.all(halo.r > 0)
        and np.all(halo.m > 0)
        and np.all(halo.rho > 0)
        and np.all(halo.p > 0)
        and np.all(halo.u > 0)
    ):
        raise FloatingPointError(f"non-positive fluid state for {model_key}")


def run_one_long(model_key: str, eta: float, output_root: str,
                 n_shells: int, t_end_gyr: float,
                 checkpoint_gyr: float, target_sigma_100=None,
                 flag_dissipation: bool = True, max_steps: int | None = None,
                 threshold_velocity_requested=None,
                 t_epsilon: float = 1e-2, r_epsilon: float = 1e-8) -> dict:
    """Evolve one model and return an execution summary.

    Existing runs are resumed from their latest ``time*.pickle`` snapshot.
    The initialization directory and run configuration are never overwritten.
    """
    if t_end_gyr <= 0 or checkpoint_gyr <= 0:
        raise ValueError("t_end_gyr and checkpoint_gyr must be positive")
    if n_shells < 16:
        raise ValueError("n_shells must be at least 16")
    if t_epsilon <= 0 or r_epsilon <= 0:
        raise ValueError("t_epsilon and r_epsilon must be positive")

    tag = model_key.split("_")[0] + f"_eta{eta:g}".replace(".", "p")
    if target_sigma_100 is not None:
        tag += f"_sigma{target_sigma_100:g}".replace(".", "p")
    if threshold_velocity_requested is not None:
        tag += f"_vstar{threshold_velocity_requested:g}".replace(".", "p")
    if not flag_dissipation:
        tag += "_elastic"
    output_dir = os.path.join(output_root, tag)
    os.makedirs(output_root, exist_ok=True)

    had_snapshots = bool(_snapshot_files(output_dir)) if os.path.isdir(output_dir) else False
    base = benchmark_models()[model_key]
    if target_sigma_100 is None:
        model = (
            parameters_at_eta(base, eta)
            if threshold_velocity_requested is None
            else parameters_at_mass_threshold(
                base, base.m_chi, threshold_velocity_requested, eta
            )
        )
    else:
        model = (
            parameters_for_target_sigma(base, eta, target_sigma_100, 100.0)
            if threshold_velocity_requested is None
            else parameters_for_target_sigma_at_threshold(
                base, eta, target_sigma_100, threshold_velocity_requested, 100.0
            )
        )
    config = _run_config(
        model_key, model, eta, target_sigma_100, n_shells, flag_dissipation,
        threshold_velocity_requested,
    )
    config["t_epsilon"] = float(t_epsilon)
    config["r_epsilon"] = float(r_epsilon)
    os.makedirs(output_dir, exist_ok=True)
    _load_or_write_config(output_dir, config)
    previous = _previous_summary(output_root, output_dir) if had_snapshots else {}

    model, halo, sigma_reference = build_halo(
        model_key, eta, output_dir, n_shells, target_sigma_100,
        flag_dissipation, threshold_velocity_requested,
    )
    halo.t_epsilon = t_epsilon
    halo.r_epsilon = r_epsilon

    if not had_snapshots:
        if halo.flag_hydrostatic_initial:
            halo.hydrostatic_adjustment()
        halo.save_halo()

    _assert_finite_state(halo, model_key)
    initial_state = halo.record.get_halo_initialization()[1]
    rho_reference = float(initial_state["rho"][3])
    u_reference = float(1.5 * initial_state["p"][3] / initial_state["rho"][3])
    rho_start = float(halo.get_central_quantity(halo.rho))
    u_start = float(halo.get_central_quantity(halo.u))
    cooling_now = getattr(halo, "C_cool", np.zeros_like(halo.rho))
    max_cooling_initial = float(
        previous.get("max_cooling_initial_code", np.max(cooling_now))
    )
    max_cooling_seen = float(
        previous.get("max_cooling_seen_code", max_cooling_initial)
    )
    scale_v_kms = float(halo.scale_v.to("km/s").value)
    v_threshold_kms = float(threshold_velocity_km_s(model))
    initial_rho = np.asarray(initial_state["rho"], dtype=float)
    initial_p = np.asarray(initial_state["p"], dtype=float)
    initial_v = np.sqrt(np.clip(initial_p / np.clip(initial_rho, 1e-300, None), 0.0, None))
    v_max_initial = float(np.max(initial_v) * scale_v_kms)
    v_central_start = float(halo.get_central_quantity(halo.v) * scale_v_kms)
    scale_t_gyr = float(halo.scale_t.to("Gyr").value)
    target_t = t_end_gyr / scale_t_gyr
    checkpoint_t = checkpoint_gyr / scale_t_gyr
    next_checkpoint = (np.floor(halo.t / checkpoint_t) + 1.0) * checkpoint_t
    steps_this_call = 0

    while halo.t < target_t:
        if max_steps is not None and steps_this_call >= max_steps:
            break
        previous_t = float(halo.t)
        halo.conduct_heat()
        halo.hydrostatic_adjustment()
        steps_this_call += 1
        _assert_finite_state(halo, model_key)
        cooling_now = getattr(halo, "C_cool", np.zeros_like(halo.rho))
        max_cooling_seen = max(max_cooling_seen, float(np.max(cooling_now)))
        if not halo.t > previous_t:
            raise FloatingPointError("non-positive time step encountered")
        if halo.t >= next_checkpoint or halo.t >= target_t:
            halo.save_halo()
            while next_checkpoint <= halo.t:
                next_checkpoint += checkpoint_t

    # Always save the current state, including a max-step partial run.
    halo.save_halo()
    rho_final = float(halo.get_central_quantity(halo.rho))
    u_final = float(halo.get_central_quantity(halo.u))
    cooling_final = getattr(halo, "C_cool", np.zeros_like(halo.rho))
    v_max_final = float(np.max(halo.v) * scale_v_kms)
    v_central_final = float(halo.get_central_quantity(halo.v) * scale_v_kms)
    diagnostic_rows = diagnose_directory(output_dir, n_scan=80)
    diagnostic_path = os.path.join(output_dir, "b1938_snapshot_diagnostics.csv")
    write_diagnostic_csv(diagnostic_rows, diagnostic_path)
    trajectory_path = os.path.join(output_dir, "fluid_trajectory.csv")
    write_trajectory_csv(fluid_trajectory(output_dir), trajectory_path)

    summary = {
        **config,
        "sigma_m_100_kms": float(sigma_reference),
        "scale_t_gyr": scale_t_gyr,
        "target_time_gyr": float(t_end_gyr),
        "current_time_gyr": float(halo.t * scale_t_gyr),
        "current_time_dimensionless": float(halo.t),
        "steps_this_call": int(steps_this_call),
        "total_conduction_steps": int(halo.n_conduction),
        "n_snapshots": len(_snapshot_files(output_dir)),
        "rho_center_initial_code": rho_reference,
        "rho_center_start_code": rho_start,
        "rho_center_final_code": rho_final,
        "rho_center_fractional_change_from_initial": rho_final / rho_reference - 1.0,
        "u_center_initial_code": u_reference,
        "u_center_start_code": u_start,
        "u_center_final_code": u_final,
        "u_center_fractional_change_from_initial": u_final / u_reference - 1.0,
        "max_cooling_initial_code": max_cooling_initial,
        "max_cooling_final_code": float(np.max(cooling_final)),
        "max_cooling_seen_code": max_cooling_seen,
        "v_threshold_kms": v_threshold_kms,
        "v_max_initial_kms": v_max_initial,
        "v_max_final_kms": v_max_final,
        "v_central_start_kms": v_central_start,
        "v_central_final_kms": v_central_final,
        "v_max_to_threshold_final": v_max_final / v_threshold_kms,
        "diagnostic_path": diagnostic_path,
        "trajectory_path": trajectory_path,
        "output_dir": output_dir,
    }
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--eta", type=float, default=0.1)
    parser.add_argument("--t-end-gyr", type=float, default=0.05)
    parser.add_argument("--checkpoint-gyr", type=float, default=0.005)
    parser.add_argument("--steps", type=int, default=None,
                        help="optional per-call step cap; runs resume on rerun")
    parser.add_argument("--n-shells", type=int, default=48)
    parser.add_argument("--target-sigma-100", type=float, default=1.0,
                        help="target sigma/m at 100 km/s; use omitted value for eta-only points")
    parser.add_argument("--threshold-velocity-kms", type=float, default=None,
                        help="optional requested massive-emission threshold velocity")
    parser.add_argument("--disable-cooling", action="store_true")
    parser.add_argument("--t-epsilon", type=float, default=1e-2)
    parser.add_argument("--r-epsilon", type=float, default=1e-8)
    parser.add_argument(
        "--output-root",
        default=os.path.join(_ROOT, "data", "P2_born_long"),
    )
    args = parser.parse_args()
    if not (0 < args.eta < 1):
        parser.error("--eta must satisfy 0 < eta < 1")
    if args.steps is not None and args.steps <= 0:
        parser.error("--steps must be positive when provided")

    rows = []
    for model_key in MODEL_KEYS:
        row = run_one_long(
            model_key, args.eta, args.output_root, args.n_shells,
            args.t_end_gyr, args.checkpoint_gyr, args.target_sigma_100,
            not args.disable_cooling, args.steps, args.threshold_velocity_kms,
            args.t_epsilon, args.r_epsilon,
        )
        rows.append(row)
        print(
            model_key,
            f"t={row['current_time_gyr']:.6g}/{args.t_end_gyr:.6g} Gyr",
            f"steps_call={row['steps_this_call']}",
            f"snapshots={row['n_snapshots']}",
        )

    summary_path = _write_summary(args.output_root, rows)
    print(f"wrote {summary_path}")


if __name__ == "__main__":
    main()
