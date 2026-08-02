"""Run short gravothermal regressions at controlled Born-valid points.

These runs exercise the complete velocity-dependent elastic transport and
direct emission-kernel cooling path. They are intentionally short numerical
regressions, not production B1938 simulations.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "..", ".."))
sys.path[:0] = [
    os.path.join(_ROOT, "src", "cross_sections"),
    _HERE,
    os.path.join(_ROOT, "external", "gravothermalsidm"),
]

from born_valid_scan import (
    parameters_at_eta,
    parameters_at_mass_threshold,
    parameters_for_target_sigma,
    parameters_for_target_sigma_at_threshold,
)
from dsidm_models import benchmark_models, born_expansion_parameter, sigma_T_born
from emission_kernel import microscopic_cooling_sigma_m
from thermal_avg import effective_cooling_sigma_m_from_kernel, effective_sigma_m_and_rdiss
from dissipative_halo import DissipativeHalo
from SourcePy.record import HaloRecord


MODEL_KEYS = ("M1_dark_photon_massive", "M2_scalar_phi_massive")


def build_halo(model_key: str, eta: float, output_dir: str, n_shells: int,
               target_sigma_100: float | None = None,
               flag_dissipation: bool = True,
               threshold_velocity_km_s: float | None = None,
               r_s_kpc: float = 3.6,
               rho_s_msun_pc3: float = 7.09e-3):
    base = benchmark_models()[model_key]
    if target_sigma_100 is None:
        model = (
            parameters_at_eta(base, eta)
            if threshold_velocity_km_s is None
            else parameters_at_mass_threshold(
                base, base.m_chi, threshold_velocity_km_s, eta
            )
        )
    else:
        model = (
            parameters_for_target_sigma(base, eta, target_sigma_100, 100.0)
            if threshold_velocity_km_s is None
            else parameters_for_target_sigma_at_threshold(
                base, eta, target_sigma_100, threshold_velocity_km_s, 100.0
            )
        )
    temperature_grid = np.logspace(0.5, 6.5, 28)
    sigma_eff, _, _, _ = effective_sigma_m_and_rdiss(
        lambda velocity: sigma_T_born(np.atleast_1d(velocity), model),
        lambda velocity: np.ones_like(np.atleast_1d(velocity), dtype=float),
        temperature_grid,
    )
    cooling_eff, cooling_grid = effective_cooling_sigma_m_from_kernel(
        lambda velocity: microscopic_cooling_sigma_m(velocity, model),
        temperature_grid,
    )
    if not np.all(np.isfinite(cooling_grid)) or np.any(cooling_grid < 0):
        raise FloatingPointError("invalid direct-cooling thermal grid")

    sigma_reference = float(sigma_T_born(np.array([100.0]), model)[0])
    record = HaloRecord(output_dir)
    halo = DissipativeHalo(
        record,
        sigma_m_eff_callable=sigma_eff,
        cooling_sigma_m_eff_callable=cooling_eff,
        flag_dissipation=flag_dissipation,
        profile="NFW",
        r_s=r_s_kpc,
        rho_s=rho_s_msun_pc3,
        sigma_m_with_units=max(sigma_reference, 1e-12),
        w_units=100.0,
        n_shells=n_shells,
        r_max=50.0,
        r_min=0.02,
        n_adjustment_max=50,
        flag_hydrostatic_initial=True,
        flag_timestep_use_relaxation=True,
        flag_timestep_use_energy=True,
    )
    halo.t_epsilon = 1e-2
    halo.r_epsilon = 1e-8
    return model, halo, sigma_reference


def run_one(model_key: str, eta: float, output_root: str,
            n_shells: int, steps: int, target_sigma_100: float | None = None,
            flag_dissipation: bool = True):
    tag = model_key.split("_")[0] + f"_eta{eta:g}".replace(".", "p")
    if target_sigma_100 is not None:
        tag += f"_sigma{target_sigma_100:g}".replace(".", "p")
    if not flag_dissipation:
        tag += "_elastic"
    output_dir = os.path.join(output_root, tag)
    if os.path.exists(output_dir):
        raise FileExistsError(
            f"refusing to overwrite existing regression directory: {output_dir}"
        )

    model, halo, sigma_reference = build_halo(
        model_key, eta, output_dir, n_shells, target_sigma_100,
        flag_dissipation,
    )
    if halo.flag_hydrostatic_initial:
        halo.hydrostatic_adjustment()
    rho_initial = float(halo.get_central_quantity(halo.rho))
    u_initial = float(halo.get_central_quantity(halo.u))
    max_cooling_initial = float(
        np.max(halo.C_cool) if hasattr(halo, "C_cool") else 0.0
    )
    halo.save_halo()

    for _ in range(steps):
        halo.conduct_heat()
        halo.hydrostatic_adjustment()
        if not (
            np.all(np.isfinite(halo.rho))
            and np.all(np.isfinite(halo.p))
            and np.all(np.isfinite(halo.u))
            and np.all(halo.rho > 0)
            and np.all(halo.p > 0)
            and np.all(halo.u > 0)
        ):
            raise FloatingPointError(f"non-finite fluid state for {model_key}")

    halo.save_halo()
    rho_final = float(halo.get_central_quantity(halo.rho))
    u_final = float(halo.get_central_quantity(halo.u))
    return {
        "model": model_key,
        "eta_B": born_expansion_parameter(model),
        "alpha_D": model.alpha_D,
        "m_chi_GeV": model.m_chi,
        "m_mediator_keV": model.m_mediator * 1e6,
        "sigma_m_100_kms": sigma_reference,
        "n_shells": n_shells,
        "steps": steps,
        "flag_dissipation": flag_dissipation,
        "time_dimensionless": float(halo.t),
        "time_gyr": float((halo.t * halo.scale_t).to("Gyr").value),
        "rho_center_initial_code": rho_initial,
        "rho_center_final_code": rho_final,
        "rho_center_fractional_change": rho_final / rho_initial - 1.0,
        "u_center_initial_code": u_initial,
        "u_center_final_code": u_final,
        "u_center_fractional_change": u_final / u_initial - 1.0,
        "max_cooling_initial_code": max_cooling_initial,
        "output_dir": output_dir,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--eta", type=float, default=0.1)
    parser.add_argument("--steps", type=int, default=25)
    parser.add_argument("--n-shells", type=int, default=48)
    parser.add_argument(
        "--target-sigma-100", type=float, default=None,
        help="adjust m_chi at fixed threshold and eta_B to hit sigma/m at 100 km/s",
    )
    parser.add_argument("--disable-cooling", action="store_true")
    parser.add_argument(
        "--output-root",
        default=os.path.join(_ROOT, "data", "P2_born_valid_regression"),
    )
    args = parser.parse_args()
    if not (0 < args.eta < 1):
        parser.error("--eta must satisfy 0 < eta < 1")
    if args.steps <= 0 or args.n_shells < 16:
        parser.error("--steps must be positive and --n-shells must be at least 16")

    os.makedirs(args.output_root, exist_ok=True)
    rows = [
        run_one(
            key, args.eta, args.output_root, args.n_shells, args.steps,
            args.target_sigma_100, not args.disable_cooling,
        )
        for key in MODEL_KEYS
    ]
    summary_path = os.path.join(args.output_root, "summary.csv")
    with open(summary_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {summary_path}")
    for row in rows:
        print(
            row["model"],
            f"eta={row['eta_B']:.3g}",
            f"sigma100={row['sigma_m_100_kms']:.6g}",
            f"t={row['time_gyr']:.6g} Gyr",
            f"delta_rho={row['rho_center_fractional_change']:.3e}",
        )


if __name__ == "__main__":
    main()
