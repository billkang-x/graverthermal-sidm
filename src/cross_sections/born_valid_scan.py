"""Scan the controlled Born region for benchmark dark-sector models.

The nominal P2 calibration fixes ``sigma_T/m`` at 50 cm^2/g and drives M1/M2
far outside the Born expansion. This utility fixes the expansion parameter
``eta_B = alpha_D * mu / m_med`` and records the resulting elastic and direct
cooling combinations at representative halo velocities.

The output is a diagnostic parameter map, not an exclusion result. Points with
``eta_B <= 0.1`` are marked ``born_controlled``; ``eta_B < 1`` is the lenient
validity mask used by the fluid runner.
"""

from __future__ import annotations

import argparse
import dataclasses
import os
import sys

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "..", ".."))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from dsidm_models import (
    KM_S_TO_NAT,
    benchmark_models,
    born_expansion_parameter,
    sigma_T_born,
)
from emission_kernel import microscopic_cooling_sigma_m


DEFAULT_ETA = (0.01, 0.03, 0.1, 0.3, 0.5, 0.9, 0.99)
DEFAULT_VELOCITIES = (50.0, 100.0, 500.0, 1000.0)
DEFAULT_MASSES = np.logspace(-1.0, 2.0, 61)


def parameters_at_eta(model, eta: float):
    """Return a copy of ``model`` with the requested Born parameter."""
    if eta <= 0 or not np.isfinite(eta):
        raise ValueError("eta must be finite and positive")
    mediator_mass = model.m_mediator
    if mediator_mass is None or mediator_mass <= 0:
        raise ValueError("the scan requires a positive elastic mediator mass")
    alpha_D = eta * mediator_mass / model.mu
    return dataclasses.replace(model, alpha_D=alpha_D)


def scan_model(model_key: str, eta_values=DEFAULT_ETA,
               velocities=DEFAULT_VELOCITIES) -> pd.DataFrame:
    """Scan one benchmark and return a tidy DataFrame."""
    models = benchmark_models()
    if model_key not in models:
        raise KeyError(f"unknown benchmark model: {model_key}")

    base = models[model_key]
    rows = []
    for eta in eta_values:
        model = parameters_at_eta(base, float(eta))
        eta_actual = born_expansion_parameter(model)
        row = {
            "model": model_key,
            "eta_B": eta_actual,
            "alpha_D": model.alpha_D,
            "born_valid": eta_actual < 1.0,
            "born_controlled": eta_actual <= 0.1 + 1e-12,
        }
        for velocity in velocities:
            tag = str(int(velocity)) if float(velocity).is_integer() else str(velocity)
            row[f"sigma_m_{tag}_kms"] = float(
                sigma_T_born(np.array([velocity]), model)[0]
            )
            row[f"cooling_sigma_m_{tag}_kms"] = float(
                microscopic_cooling_sigma_m(velocity, model)[0]
            )
        rows.append(row)
    return pd.DataFrame(rows)


def scan_all_models(eta_values=DEFAULT_ETA,
                    velocities=DEFAULT_VELOCITIES) -> pd.DataFrame:
    frames = [scan_model(key, eta_values, velocities)
              for key in benchmark_models()]
    return pd.concat(frames, ignore_index=True)


def threshold_velocity_km_s(model) -> float:
    """Return the massive-emission threshold velocity."""
    if model.m_mediator <= 0:
        raise ValueError("threshold velocity requires a positive emitted mass")
    return float(np.sqrt(2.0 * model.m_mediator / model.mu) / KM_S_TO_NAT)


def parameters_at_mass_threshold(model, m_chi: float, v_star_km_s: float,
                                 eta: float):
    """Change the DM mass while preserving threshold velocity and ``eta_B``."""
    if m_chi <= 0 or v_star_km_s <= 0:
        raise ValueError("m_chi and v_star_km_s must be positive")
    mu = 0.5 * float(m_chi)
    mediator_mass = 0.5 * mu * (v_star_km_s * KM_S_TO_NAT) ** 2
    alpha_D = eta * mediator_mass / mu
    return dataclasses.replace(
        model,
        m_chi=float(m_chi),
        m_mediator=float(mediator_mass),
        alpha_D=float(alpha_D),
    )


def parameters_for_target_sigma(model, eta: float, target_sigma_m: float,
                                velocity_km_s: float = 100.0):
    """Solve for ``m_chi`` at fixed threshold and Born parameter.

    At fixed ``eta_B`` and threshold velocity the dimensionless Yukawa shape
    is unchanged and ``sigma_T/m_chi`` scales exactly as ``m_chi**-3``.
    """
    if target_sigma_m <= 0 or velocity_km_s <= 0:
        raise ValueError("target_sigma_m and velocity_km_s must be positive")
    v_star = threshold_velocity_km_s(model)
    reference = parameters_at_mass_threshold(
        model, model.m_chi, v_star, eta
    )
    sigma_reference = float(
        sigma_T_born(np.array([velocity_km_s]), reference)[0]
    )
    solved_mass = model.m_chi * (sigma_reference / target_sigma_m) ** (1.0 / 3.0)
    solved = parameters_at_mass_threshold(model, solved_mass, v_star, eta)
    return solved


def parameters_for_target_sigma_at_threshold(
    model, eta: float, target_sigma_m: float,
    v_star_km_s: float, velocity_km_s: float = 100.0,
):
    """Solve for ``m_chi`` at a requested threshold and Born parameter."""
    if v_star_km_s <= 0:
        raise ValueError("v_star_km_s must be positive")
    reference = parameters_at_mass_threshold(
        model, model.m_chi, v_star_km_s, eta
    )
    sigma_reference = float(
        sigma_T_born(np.array([velocity_km_s]), reference)[0]
    )
    solved_mass = model.m_chi * (sigma_reference / target_sigma_m) ** (1.0 / 3.0)
    return parameters_at_mass_threshold(model, solved_mass, v_star_km_s, eta)


def scan_mass_hierarchy(masses=DEFAULT_MASSES, eta_values=(0.1, 0.99),
                        velocity=100.0) -> pd.DataFrame:
    """Scan DM mass at fixed benchmark threshold velocity."""
    rows = []
    models = benchmark_models()
    for model_key in ("M1_dark_photon_massive", "M2_scalar_phi_massive"):
        base = models[model_key]
        v_star = threshold_velocity_km_s(base)
        for eta in eta_values:
            for m_chi in masses:
                model = parameters_at_mass_threshold(base, m_chi, v_star, eta)
                sigma_m = float(sigma_T_born(np.array([velocity]), model)[0])
                rows.append({
                    "model": model_key,
                    "eta_target": float(eta),
                    "eta_B": born_expansion_parameter(model),
                    "m_chi_GeV": model.m_chi,
                    "m_mediator_keV": model.m_mediator * 1e6,
                    "v_star_km_s": v_star,
                    "alpha_D": model.alpha_D,
                    "velocity_km_s": velocity,
                    "sigma_m_cm2_g": sigma_m,
                    "astrophysical_overlap": 0.1 <= sigma_m <= 10.0,
                })
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default=os.path.join(_ROOT, "data", "born_valid_parameter_scan.csv"),
        help="CSV output path",
    )
    parser.add_argument(
        "--mass-output",
        default=os.path.join(_ROOT, "data", "born_valid_mass_hierarchy_scan.csv"),
        help="fixed-threshold mass hierarchy CSV output path",
    )
    args = parser.parse_args()
    frame = scan_all_models()
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    frame.to_csv(args.output, index=False)
    print(f"wrote {len(frame)} rows to {args.output}")
    for key, group in frame.groupby("model", sort=False):
        max_sigma = group.loc[group["born_valid"], "sigma_m_100_kms"].max()
        print(f"{key}: sampled max sigma_m(100 km/s), eta_B<1 = {max_sigma:.6g} cm^2/g")

    mass_frame = scan_mass_hierarchy()
    mass_frame.to_csv(args.mass_output, index=False)
    print(f"wrote {len(mass_frame)} rows to {args.mass_output}")
    overlap = mass_frame[mass_frame["astrophysical_overlap"]]
    for (key, eta), group in overlap.groupby(["model", "eta_target"], sort=False):
        print(
            f"{key}, eta_B={eta:g}: astrophysical overlap at "
            f"m_chi={group['m_chi_GeV'].min():.3g}-"
            f"{group['m_chi_GeV'].max():.3g} GeV"
        )


if __name__ == "__main__":
    main()
