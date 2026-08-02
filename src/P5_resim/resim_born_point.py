"""Directly resimulate a Born-valid B1938 candidate in physical units.

The candidate is selected from the snapshot diagnostics of a long fluid run.
Its elastic rescaling parameters provide the physical NFW initial profile and
target evolution time.  The new halo is then evolved with the unchanged
velocity-dependent microscopic model, so the result does not assume that the
elastic rescaling symmetry remains valid.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
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
    os.path.join(_ROOT, "external", "gravothermalsidm"),
]

from born_valid_scan import threshold_velocity_km_s
from rescale import (
    MASS_RATIO_ERR,
    MASS_RATIO_OBS,
    M_20PC,
    M_20PC_ERR,
    M_90PC,
    M_90PC_ERR,
    RHO0_SIM,
    R_INNER_PC,
    absolute_mass_fit,
    projected_enclosed_mass,
)
from run_born_valid import build_halo


def select_candidate(diagnostic_path: str, snapshot_idx: int | None = None,
                     include_initial: bool = False) -> dict:
    """Select the minimum-chi-square evolved row, optionally by snapshot."""
    frame = pd.read_csv(diagnostic_path)
    frame = frame[frame["snapshot_time_gyr"] >= (0 if include_initial else np.nextafter(0.0, 1.0))].copy()
    if snapshot_idx is not None:
        frame = frame[frame["snapshot_idx"] == snapshot_idx]
    if frame.empty:
        raise ValueError("no evolved B1938 diagnostic candidate is available")
    return frame.loc[frame["mass_chi2"].idxmin()].to_dict()


def physical_candidate(candidate: dict, run_dir: str,
                       sigma_m_sim: float) -> dict:
    """Convert a ratio-selected snapshot row to its elastic physical guess."""
    from SourcePy.record import HaloRecord

    record = HaloRecord(run_dir)
    halo_ini, _ = record.get_halo_initialization()
    r_s_sim_kpc = float(halo_ini["r_s"])
    rho_s_sim = float(halo_ini["rho_s"])
    r2d_rs = float(candidate["r2D_rs"])
    mu = float(candidate["mass_mu"])
    lam = R_INNER_PC / (r2d_rs * r_s_sim_kpc * 1000.0)
    t_sim_gyr = float(candidate["snapshot_time_gyr"])
    return {
        "lambda": lam,
        "mu": mu,
        "r_s_sim_kpc": r_s_sim_kpc,
        "rho_s_sim_msun_pc3": rho_s_sim,
        "r_s_phys_kpc": lam * r_s_sim_kpc,
        "rho_s_phys_msun_pc3": mu * rho_s_sim / lam**3,
        "sigma_m_elastic_guess_cm2_g": lam**2 / mu * sigma_m_sim,
        "t_evo_phys_gyr": float(np.sqrt(lam**3 / mu) * t_sim_gyr),
    }


def _assert_state(halo) -> None:
    arrays = (halo.r, halo.m, halo.rho, halo.p, halo.u)
    if not all(np.all(np.isfinite(value)) for value in arrays):
        raise FloatingPointError("non-finite direct-resimulation state")
    if not all(np.all(value > 0) for value in arrays):
        raise FloatingPointError("non-positive direct-resimulation state")


def run_direct_resimulation(
    source_run_dir: str, output_dir: str, n_shells: int = 96,
    snapshot_idx: int | None = None, max_steps: int = 100_000,
    t_epsilon: float = 5e-3, r_epsilon: float = 1e-8,
    include_initial: bool = False, flag_dissipation: bool = True,
    diagnostic_filename: str = "b1938_snapshot_diagnostics.csv",
) -> dict:
    """Run one direct physical-parameter resimulation."""
    source_run_dir = os.path.abspath(source_run_dir)
    diagnostic_path = os.path.join(source_run_dir, diagnostic_filename)
    config_path = os.path.join(source_run_dir, "run_config.json")
    if not os.path.isfile(diagnostic_path) or not os.path.isfile(config_path):
        raise FileNotFoundError("source run lacks diagnostics or run_config.json")
    if os.path.exists(output_dir):
        raise FileExistsError(f"refusing to overwrite direct run: {output_dir}")

    with open(config_path, encoding="utf-8") as handle:
        config = json.load(handle)
    candidate = select_candidate(diagnostic_path, snapshot_idx, include_initial)
    sigma_m_sim = float(config.get("target_sigma_m_100_kms") or 0.0)
    if sigma_m_sim <= 0:
        raise ValueError("direct resimulation requires a positive target sigma/m")
    physical = physical_candidate(candidate, source_run_dir, sigma_m_sim)

    model_key = str(config["model"])
    eta = float(config["eta_requested"])
    threshold = config.get("threshold_velocity_requested_kms")
    os.makedirs(os.path.dirname(os.path.abspath(output_dir)), exist_ok=True)
    model, halo, sigma_reference = build_halo(
        model_key, eta, output_dir, n_shells, sigma_m_sim, flag_dissipation,
        threshold,
        physical["r_s_phys_kpc"], physical["rho_s_phys_msun_pc3"],
    )
    halo.t_epsilon = t_epsilon
    halo.r_epsilon = r_epsilon
    if halo.flag_hydrostatic_initial:
        halo.hydrostatic_adjustment()
    halo.save_halo()
    _assert_state(halo)

    scale_t_gyr = float(halo.scale_t.to("Gyr").value)
    target_t = physical["t_evo_phys_gyr"] / scale_t_gyr
    max_cooling_seen = float(np.max(getattr(halo, "C_cool", np.zeros_like(halo.rho))))
    completed = False
    steps = 0
    while halo.t < target_t and steps < max_steps:
        previous_t = float(halo.t)
        halo.conduct_heat()
        halo.hydrostatic_adjustment()
        steps += 1
        _assert_state(halo)
        if not halo.t > previous_t:
            raise FloatingPointError("non-positive direct-resimulation time step")
        max_cooling_seen = max(
            max_cooling_seen,
            float(np.max(getattr(halo, "C_cool", np.zeros_like(halo.rho)))),
        )
    completed = halo.t >= target_t
    halo.save_halo()

    r_kpc = halo.r * halo.scale_r.to("kpc").value
    rho_phys = halo.rho * halo.scale_rho.to("Msun/pc**3").value
    m_inner = projected_enclosed_mass(r_kpc, rho_phys, 0.02)
    m_outer = projected_enclosed_mass(r_kpc, rho_phys, 0.09)
    ratio = float(m_inner / m_outer) if m_outer > 0 else np.nan
    residual_inner = float((m_inner - M_20PC) / M_20PC_ERR)
    residual_outer = float((m_outer - M_90PC) / M_90PC_ERR)
    direct_chi2 = residual_inner**2 + residual_outer**2
    common_scale_fit = absolute_mass_fit(m_inner, m_outer)

    result = {
        "model": model_key,
        "source_run_dir": source_run_dir,
        "source_snapshot_idx": int(candidate["snapshot_idx"]),
        "source_snapshot_time_gyr": float(candidate["snapshot_time_gyr"]),
        "source_mass_chi2": float(candidate["mass_chi2"]),
        **physical,
        "eta_B": float(config["eta_B"]),
        "m_chi_GeV": float(model.m_chi),
        "m_mediator_keV": float(model.m_mediator * 1e6),
        "v_threshold_kms": float(threshold_velocity_km_s(model)),
        "sigma_m_100_kms": float(sigma_reference),
        "flag_dissipation": bool(flag_dissipation),
        "n_shells": int(n_shells),
        "t_epsilon": float(t_epsilon),
        "steps": int(steps),
        "completed": bool(completed),
        "t_final_gyr": float(halo.t * scale_t_gyr),
        "v_max_final_kms": float(np.max(halo.v) * halo.scale_v.to("km/s").value),
        "max_cooling_seen_code": max_cooling_seen,
        "M_inner_direct_msun": float(m_inner),
        "M_outer_direct_msun": float(m_outer),
        "mass_ratio_direct": ratio,
        "mass_ratio_offset_sigma": float((ratio - MASS_RATIO_OBS) / MASS_RATIO_ERR),
        "mass_residual_inner_sigma": residual_inner,
        "mass_residual_outer_sigma": residual_outer,
        "mass_chi2_direct": float(direct_chi2),
        "posthoc_common_mass_scale": common_scale_fit["mu"],
        "posthoc_common_scale_chi2": common_scale_fit["chi2"],
        "output_dir": os.path.abspath(output_dir),
    }
    summary_path = os.path.join(output_dir, "summary.csv")
    with open(summary_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(result))
        writer.writeheader()
        writer.writerow(result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_run_dir")
    parser.add_argument("output_dir")
    parser.add_argument("--snapshot-idx", type=int, default=None)
    parser.add_argument("--n-shells", type=int, default=96)
    parser.add_argument("--max-steps", type=int, default=100_000)
    parser.add_argument("--t-epsilon", type=float, default=5e-3)
    parser.add_argument("--r-epsilon", type=float, default=1e-8)
    parser.add_argument(
        "--include-initial", action="store_true",
        help="allow an explicit snapshot_idx=0 candidate for t=0 diagnostics",
    )
    parser.add_argument(
        "--disable-cooling", action="store_true",
        help="run the same physical halo with the dissipative source disabled",
    )
    parser.add_argument(
        "--diagnostic-filename", default="b1938_snapshot_diagnostics.csv",
        help="candidate diagnostic CSV inside the source run directory",
    )
    args = parser.parse_args()
    result = run_direct_resimulation(
        args.source_run_dir, args.output_dir, args.n_shells,
        args.snapshot_idx, args.max_steps, args.t_epsilon, args.r_epsilon,
        args.include_initial, not args.disable_cooling, args.diagnostic_filename,
    )
    print(
        f"{result['model']}: completed={result['completed']} "
        f"t={result['t_final_gyr']:.6g}/{result['t_evo_phys_gyr']:.6g} Gyr "
        f"chi2_direct={result['mass_chi2_direct']:.6g} "
        f"ratio={result['mass_ratio_direct']:.6g}"
    )


if __name__ == "__main__":
    main()
