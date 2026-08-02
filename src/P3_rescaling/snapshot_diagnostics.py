"""Snapshot-level projected-mass diagnostics for Born-valid fluid runs."""

from __future__ import annotations

import argparse
import csv
import os
import sys

import numpy as np
from astropy import constants as ct
from astropy import units as ut
from scipy.optimize import minimize_scalar

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(_ROOT, "external", "gravothermalsidm"))

from SourcePy.record import HaloRecord
from rescale import (
    MASS_RATIO_ERR,
    MASS_RATIO_OBS,
    RS_SIM_KPC,
    absolute_mass_fit,
    find_matching_radii,
    projected_enclosed_mass,
)


def scales_from_initialization(dir_data: str) -> tuple[float, float, float]:
    """Return simulation radius, density, and time scales in useful units."""
    record = HaloRecord(dir_data)
    halo_ini, _ = record.get_halo_initialization()
    r_s = float(halo_ini["r_s"]) * ut.kpc
    rho_s = float(halo_ini["rho_s"]) * ut.M_sun / ut.pc**3
    scale_t = (1.0 / np.sqrt(4.0 * np.pi * rho_s * ct.G)).to("Gyr").value
    return r_s.to_value("kpc"), rho_s.to_value(ut.M_sun / ut.pc**3), float(scale_t)


def diagnose_snapshot(data: dict, scale_r_kpc: float, scale_rho: float,
                      scale_t_gyr: float, r_s_kpc: float = RS_SIM_KPC,
                      snapshot_idx: int = 0, n_scan: int = 80) -> list[dict]:
    """Return ratio-preselected rows plus a continuous best-fit refinement."""
    r_kpc = np.asarray(data["r"], dtype=float) * scale_r_kpc
    rho = np.asarray(data["rho"], dtype=float) * scale_rho
    t_gyr = float(data["t"]) * scale_t_gyr
    scan = find_matching_radii(r_kpc, rho, r_s_kpc, n_scan=n_scan)

    def make_row(r2d_rs: float, m_inner: float | None = None,
                 m_outer: float | None = None, ratio: float | None = None) -> dict | None:
        if m_inner is None or m_outer is None:
            inner_kpc = r2d_rs * r_s_kpc
            m_inner = projected_enclosed_mass(r_kpc, rho, inner_kpc)
            m_outer = projected_enclosed_mass(r_kpc, rho, 4.5 * inner_kpc)
            ratio = m_inner / m_outer if m_outer > 0 else np.nan
        if not (np.isfinite(ratio) and m_inner > 0 and m_outer > 0):
            return None
        n_sigma = abs(ratio - MASS_RATIO_OBS) / MASS_RATIO_ERR
        if n_sigma > 3.0:
            return None
        fit = absolute_mass_fit(m_inner, m_outer)
        return {
            "snapshot_idx": int(snapshot_idx),
            "snapshot_time_gyr": t_gyr,
            "r2D_rs": float(r2d_rs),
            "M_inner_sim_msun": float(m_inner),
            "M_outer_sim_msun": float(m_outer),
            "mass_ratio": float(ratio),
            "ratio_nsigma": float(n_sigma),
            "mass_mu": fit["mu"],
            "M_inner_phys_msun": fit["M_inner_phys"],
            "M_outer_phys_msun": fit["M_outer_phys"],
            "mass_chi2": fit["chi2"],
            "mass_nsigma_dof1": fit["nsigma_dof1"],
            "mass_residual_inner_sigma": fit["residual_inner_sigma"],
            "mass_residual_outer_sigma": fit["residual_outer_sigma"],
        }

    rows = []
    grid_chi2 = []
    for r2d_rs, m_inner, m_outer, ratio in scan:
        row = make_row(float(r2d_rs), float(m_inner), float(m_outer), float(ratio))
        if row is not None:
            rows.append(row)
            grid_chi2.append((float(row["mass_chi2"]), float(r2d_rs)))

    # The grid locates the correct branch, while the bounded refinement avoids
    # candidate jumps when two adjacent r2D/r_s samples exchange rank.
    if grid_chi2:
        _, seed = min(grid_chi2)
        grid_radius = np.asarray(scan[:, 0], dtype=float)
        seed_index = int(np.argmin(np.abs(np.log(grid_radius / seed))))
        lower_index = max(0, seed_index - 1)
        upper_index = min(len(grid_radius) - 1, seed_index + 1)

        def objective(log_radius: float) -> float:
            row = make_row(float(np.exp(log_radius)))
            return float(row["mass_chi2"]) if row is not None else 1e12

        refined = minimize_scalar(
            objective,
            bounds=(np.log(grid_radius[lower_index]), np.log(grid_radius[upper_index])),
            method="bounded",
            options={"xatol": 1e-12},
        )
        if refined.success:
            refined_row = make_row(float(np.exp(refined.x)))
            if refined_row is not None:
                rows.append(refined_row)
    return rows


def diagnose_directory(dir_data: str, n_scan: int = 80,
                       snapshot_stride: int = 1) -> list[dict]:
    """Process saved snapshots in time order."""
    if snapshot_stride <= 0:
        raise ValueError("snapshot_stride must be positive")
    record = HaloRecord(dir_data)
    files, _ = record.glob_pickle_files()
    if not files:
        return []
    scale_r_kpc, scale_rho, scale_t_gyr = scales_from_initialization(dir_data)
    rows = []
    for snapshot_idx, file_name in enumerate(files[::snapshot_stride]):
        data = record.get_halo_state_pickled(file_halo=file_name)
        if not data or np.any(~np.isfinite(data.get("rho", np.array([np.nan])))):
            continue
        rows.extend(diagnose_snapshot(
            data, scale_r_kpc, scale_rho, scale_t_gyr,
            r_s_kpc=scale_r_kpc, snapshot_idx=snapshot_idx, n_scan=n_scan,
        ))
    return rows


def _log_interp_profile(r: np.ndarray, value: np.ndarray,
                        target: float) -> float:
    """Interpolate a positive shell profile at a dimensionless radius."""
    if target <= r[0]:
        return float(value[0])
    if target >= r[-1]:
        return float(value[-1])
    return float(np.exp(np.interp(
        np.log(target), np.log(r), np.log(np.clip(value, 1e-300, None))
    )))


def fluid_trajectory(dir_data: str, snapshot_stride: int = 1,
                     radii_over_rs: tuple[float, ...] = (0.1, 0.2, 0.5)) -> list[dict]:
    """Return fixed-radius fluid diagnostics for every saved snapshot."""
    if snapshot_stride <= 0:
        raise ValueError("snapshot_stride must be positive")
    record = HaloRecord(dir_data)
    files, _ = record.glob_pickle_files()
    if not files:
        return []
    r_s_kpc, rho_scale, scale_t_gyr = scales_from_initialization(dir_data)
    halo_ini, _ = record.get_halo_initialization()
    rho_s = float(halo_ini["rho_s"]) * ut.M_sun / ut.pc**3
    r_s = float(halo_ini["r_s"]) * ut.kpc
    scale_v_kms = np.sqrt(ct.G * 4.0 * np.pi * rho_s * r_s**2).to("km/s").value
    rows = []
    for snapshot_idx, file_name in enumerate(files[::snapshot_stride]):
        data = record.get_halo_state_pickled(file_halo=file_name)
        if not data:
            continue
        r = np.asarray(data["r"], dtype=float)
        rho = np.asarray(data["rho"], dtype=float)
        pressure = np.asarray(data["p"], dtype=float)
        u_profile = 1.5 * pressure / np.clip(rho, 1e-300, None)
        v_profile = np.sqrt(np.clip(pressure / np.clip(rho, 1e-300, None), 0.0, None))
        row = {
            "snapshot_idx": int(snapshot_idx),
            "snapshot_time_gyr": float(data["t"] * scale_t_gyr),
            "rho_center_shell3_code": float(rho[3]),
            "u_center_shell3_code": float(u_profile[3]),
            "v_center_shell3_kms": float(v_profile[3] * scale_v_kms),
            "v_max_kms": float(np.max(v_profile) * scale_v_kms),
        }
        for radius_over_rs in radii_over_rs:
            label = str(radius_over_rs).replace(".", "p")
            rho_value = _log_interp_profile(r, rho, radius_over_rs)
            u_value = _log_interp_profile(r, u_profile, radius_over_rs)
            row[f"rho_r{label}rs_code"] = rho_value
            row[f"u_r{label}rs_code"] = u_value
        rows.append(row)
    return rows


def write_trajectory_csv(rows: list[dict], path: str) -> None:
    """Write fixed-radius trajectory diagnostics."""
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    if not fields:
        fields = ["snapshot_idx", "snapshot_time_gyr"]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def best_mass_fit(rows: list[dict]) -> dict | None:
    """Return the minimum-chi-square row, or None when no ratio match exists."""
    if not rows:
        return None
    return min(rows, key=lambda row: float(row["mass_chi2"]))


def write_diagnostic_csv(rows: list[dict], path: str) -> None:
    """Write diagnostics with a stable header, including empty outputs."""
    fields = [
        "snapshot_idx", "snapshot_time_gyr", "r2D_rs",
        "M_inner_sim_msun", "M_outer_sim_msun", "mass_ratio", "ratio_nsigma",
        "mass_mu", "M_inner_phys_msun", "M_outer_phys_msun", "mass_chi2",
        "mass_nsigma_dof1", "mass_residual_inner_sigma",
        "mass_residual_outer_sigma",
    ]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory")
    parser.add_argument("--output", default=None)
    parser.add_argument("--n-scan", type=int, default=80)
    args = parser.parse_args()
    rows = diagnose_directory(args.directory, n_scan=args.n_scan)
    output = args.output or os.path.join(args.directory, "b1938_snapshot_diagnostics.csv")
    write_diagnostic_csv(rows, output)
    best = best_mass_fit(rows)
    print(f"wrote {output} ({len(rows)} ratio-preselected rows)")
    if best is not None:
        print(
            f"best snapshot={best['snapshot_idx']} "
            f"t={best['snapshot_time_gyr']:.6g} Gyr "
            f"chi2={best['mass_chi2']:.6g} "
            f"r2D/rs={best['r2D_rs']:.6g}"
        )


if __name__ == "__main__":
    main()
