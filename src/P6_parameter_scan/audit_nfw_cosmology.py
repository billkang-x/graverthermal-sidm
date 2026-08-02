"""Quantify the isolated-NFW extrapolation of a direct B1938 initial fit.

The inferred concentration is an extrapolation diagnostic, not a subhalo
concentration measurement: the direct fluid domain is truncated and does not
model a cosmological infall or tidal history.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from scipy.optimize import brentq


def nfw_mass_factor(concentration: float) -> float:
    return float(np.log1p(concentration) - concentration / (1.0 + concentration))


def isolated_nfw_extrapolation(
    r_s_pc: float,
    rho_s_msun_pc3: float,
    redshift: float,
    *,
    h0_km_s_mpc: float = 67.4,
    omega_m: float = 0.315,
    overdensity: float = 200.0,
    modeled_rmax_over_rs: float = 50.0,
) -> dict[str, float]:
    """Map ``(r_s, rho_s)`` to an isolated spherical-overdensity halo."""
    values = (r_s_pc, rho_s_msun_pc3, redshift, h0_km_s_mpc, omega_m, overdensity)
    if not all(np.isfinite(values)) or r_s_pc <= 0 or rho_s_msun_pc3 <= 0:
        raise ValueError("NFW parameters and cosmology must be finite and physical")
    if redshift < 0 or h0_km_s_mpc <= 0 or not 0 < omega_m < 1 or overdensity <= 0:
        raise ValueError("NFW parameters and cosmology must be finite and physical")

    cosmology = FlatLambdaCDM(H0=h0_km_s_mpc, Om0=omega_m, Tcmb0=2.7255)
    rho_critical = float(cosmology.critical_density(redshift).to_value(u.Msun / u.pc**3))
    characteristic_overdensity = rho_s_msun_pc3 / rho_critical

    def residual(concentration: float) -> float:
        predicted = (
            overdensity / 3.0 * concentration**3 / nfw_mass_factor(concentration)
        )
        return predicted - characteristic_overdensity

    concentration = float(brentq(residual, 1e-4, 1e6))
    r_delta_pc = concentration * r_s_pc
    m_delta_msun = (
        4.0 * np.pi / 3.0 * overdensity * rho_critical * r_delta_pc**3
    )
    return {
        "redshift": redshift,
        "H0_km_s_Mpc": h0_km_s_mpc,
        "Omega_m": omega_m,
        "overdensity_critical": overdensity,
        "rho_critical_msun_pc3": rho_critical,
        "r_s_pc": r_s_pc,
        "rho_s_msun_pc3": rho_s_msun_pc3,
        "c_delta_isolated_extrapolation": concentration,
        "r_delta_pc_isolated_extrapolation": r_delta_pc,
        "M_delta_msun_isolated_extrapolation": m_delta_msun,
        "modeled_rmax_over_rs": modeled_rmax_over_rs,
        "r_delta_over_modeled_rmax": concentration / modeled_rmax_over_rs,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--direct-summary", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--redshift", type=float, default=0.881)
    parser.add_argument("--h0", type=float, default=67.4)
    parser.add_argument("--omega-m", type=float, default=0.315)
    parser.add_argument("--modeled-rmax-over-rs", type=float, default=50.0)
    args = parser.parse_args()

    source = pd.read_csv(args.direct_summary).iloc[0]
    result = isolated_nfw_extrapolation(
        1000.0 * float(source["r_s_phys_kpc"]),
        float(source["rho_s_phys_msun_pc3"]),
        args.redshift,
        h0_km_s_mpc=args.h0,
        omega_m=args.omega_m,
        modeled_rmax_over_rs=args.modeled_rmax_over_rs,
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([result]).to_csv(output, index=False)
    print(f"wrote isolated-NFW extrapolation to {output}")
    print(pd.Series(result).to_string())


if __name__ == "__main__":
    main()
