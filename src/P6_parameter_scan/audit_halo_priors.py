"""Audit Lambda-CDM concentration tension with a broad tidal envelope.

The initial two-mass NFW fit is not assigned an isolated-halo concentration
prior directly.  Instead, this script forward models an infall-normalized NFW
halo, applies the smooth Baltz-Marshall-Oguri (BMO) truncation

    rho(r) = rho_NFW(r) [r_t^2 / (r^2 + r_t^2)]^2,

and profiles the two B1938 projected masses over halo mass and
``tau = r_t/r_s``.  The Dutton-Maccio Planck concentration relation supplies
the median and its 0.11-dex intrinsic scatter.  Because no host orbit is used,
the truncation interval is a deliberately generous sensitivity envelope, not
a cosmological measurement of the tidal radius.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from scipy.optimize import differential_evolution


MASS_OBS = np.array([4.25e5, 1.167e6], dtype=float)
MASS_ERR = np.array([0.21e5, 0.039e6], dtype=float)
PROJECTED_RADII_PC = np.array([20.0, 90.0], dtype=float)
DUTTON_SCATTER_DEX = 0.11


def dutton_maccio_c200(
    mass_msun: np.ndarray | float,
    redshift: float,
    *,
    h: float = 0.674,
) -> np.ndarray | float:
    """Return the Planck NFW median concentration from Dutton & Maccio.

    This is their redshift-dependent fit for ``c_200``:
    ``log10(c_200) = a + b log10[M_200/(1e12 h^-1 Msun)]`` with
    ``b=-0.101+0.026z`` and
    ``a=0.520+(0.905-0.520) exp(-0.617 z^1.21)``.
    """
    mass = np.asarray(mass_msun, dtype=float)
    if np.any(~np.isfinite(mass)) or np.any(mass <= 0):
        raise ValueError("mass_msun must be finite and positive")
    if not np.isfinite(redshift) or redshift < 0 or not 0 < h < 2:
        raise ValueError("redshift and h must be physical")
    a = 0.520 + (0.905 - 0.520) * np.exp(-0.617 * redshift**1.21)
    b = -0.101 + 0.026 * redshift
    concentration = 10.0 ** (a + b * np.log10(mass * h / 1.0e12))
    if np.ndim(mass_msun) == 0:
        return float(concentration)
    return concentration


def nfw_scales_from_m200_c200(
    mass_msun: float,
    concentration: float,
    rho_critical_msun_pc3: float,
) -> tuple[float, float, float]:
    """Return ``(r200, r_s, rho_s)`` in pc and Msun/pc^3."""
    values = (mass_msun, concentration, rho_critical_msun_pc3)
    if not all(np.isfinite(values)) or any(value <= 0 for value in values):
        raise ValueError("mass, concentration, and critical density must be positive")
    r200_pc = (3.0 * mass_msun / (4.0 * np.pi * 200.0 * rho_critical_msun_pc3)) ** (1.0 / 3.0)
    r_s_pc = r200_pc / concentration
    mass_factor = np.log1p(concentration) - concentration / (1.0 + concentration)
    rho_s = mass_msun / (4.0 * np.pi * r_s_pc**3 * mass_factor)
    return float(r200_pc), float(r_s_pc), float(rho_s)


def _cylindrical_solid_angle_fraction(x: np.ndarray, projected_x: float) -> np.ndarray:
    """Fraction of a spherical shell lying inside a projected cylinder."""
    fraction = np.ones_like(x)
    outside = x > projected_x
    q = (projected_x / x[outside]) ** 2
    fraction[outside] = q / (1.0 + np.sqrt(np.clip(1.0 - q, 0.0, None)))
    return fraction


def bmo_projected_masses(
    mass_msun: float,
    concentration: float,
    tau: float,
    rho_critical_msun_pc3: float,
    *,
    radial_count: int = 1800,
) -> tuple[np.ndarray, dict[str, float]]:
    """Project the smoothly truncated NFW profile at 20 and 90 pc."""
    if not np.isfinite(tau) or tau <= 0:
        raise ValueError("tau must be finite and positive")
    if radial_count < 400:
        raise ValueError("radial_count must be at least 400")
    r200_pc, r_s_pc, rho_s = nfw_scales_from_m200_c200(
        mass_msun, concentration, rho_critical_msun_pc3
    )
    x = np.logspace(-6.0, 4.0, radial_count)
    bmo_factor = (tau**2 / (x**2 + tau**2)) ** 2
    radial_mass_kernel = x * bmo_factor / (1.0 + x) ** 2
    normalization = 4.0 * np.pi * rho_s * r_s_pc**3
    masses = []
    for radius_pc in PROJECTED_RADII_PC:
        fraction = _cylindrical_solid_angle_fraction(x, radius_pc / r_s_pc)
        integrate = getattr(np, "trapezoid", None)
        if integrate is None:
            from scipy.integrate import trapezoid as integrate
        masses.append(normalization * integrate(radial_mass_kernel * fraction, x))
    return np.asarray(masses), {
        "r200_pc": r200_pc,
        "r_s_pc": r_s_pc,
        "rho_s_msun_pc3": rho_s,
        "r_t_pc": tau * r_s_pc,
    }


def projected_mass_chi2(masses: np.ndarray) -> float:
    masses = np.asarray(masses, dtype=float)
    if masses.shape != MASS_OBS.shape or np.any(~np.isfinite(masses)):
        raise ValueError("masses must contain the two finite projected masses")
    return float(np.sum(((masses - MASS_OBS) / MASS_ERR) ** 2))


def profile_at_concentration_offset(
    z_concentration: float,
    redshift: float,
    rho_critical_msun_pc3: float,
    *,
    h: float = 0.674,
    scatter_dex: float = DUTTON_SCATTER_DEX,
    log_mass_bounds: tuple[float, float] = (5.5, 10.0),
    log_tau_bounds: tuple[float, float] = (0.0, 2.0),
    seed: int = 1938666,
) -> dict[str, float]:
    """Profile the two-mass likelihood over mass and BMO truncation."""
    if scatter_dex <= 0:
        raise ValueError("scatter_dex must be positive")

    def objective(parameters: np.ndarray) -> float:
        log_mass, log_tau = parameters
        mass = 10.0**log_mass
        concentration = float(
            dutton_maccio_c200(mass, redshift, h=h)
            * 10.0 ** (scatter_dex * z_concentration)
        )
        masses, _ = bmo_projected_masses(
            mass, concentration, 10.0**log_tau, rho_critical_msun_pc3
        )
        return projected_mass_chi2(masses)

    result = differential_evolution(
        objective,
        bounds=(log_mass_bounds, log_tau_bounds),
        seed=seed,
        popsize=10,
        maxiter=90,
        tol=2e-7,
        polish=True,
        workers=1,
        updating="immediate",
    )
    log_mass, log_tau = result.x
    mass = float(10.0**log_mass)
    tau = float(10.0**log_tau)
    median_concentration = float(dutton_maccio_c200(mass, redshift, h=h))
    concentration = float(median_concentration * 10.0 ** (scatter_dex * z_concentration))
    masses, scales = bmo_projected_masses(
        mass, concentration, tau, rho_critical_msun_pc3
    )
    return {
        "concentration_offset_sigma": float(z_concentration),
        "profile_mass_chi2": projected_mass_chi2(masses),
        "M200_msun": mass,
        "c200_median": median_concentration,
        "c200_candidate": concentration,
        "tau_rt_over_rs": tau,
        "M2D_20pc_msun": float(masses[0]),
        "M2D_90pc_msun": float(masses[1]),
        **scales,
    }


def profile_grid(
    offsets: np.ndarray,
    redshift: float,
    rho_critical_msun_pc3: float,
    **kwargs,
) -> pd.DataFrame:
    rows = [
        profile_at_concentration_offset(
            float(offset), redshift, rho_critical_msun_pc3,
            seed=1938666 + index, **kwargs,
        )
        for index, offset in enumerate(np.asarray(offsets, dtype=float))
    ]
    return pd.DataFrame(rows)


def _best_penalized_candidate(profile: pd.DataFrame) -> pd.Series:
    score = (
        profile["profile_mass_chi2"].to_numpy(dtype=float)
        + profile["concentration_offset_sigma"].to_numpy(dtype=float) ** 2
    )
    return profile.iloc[int(np.argmin(score))]


def summarize_prior_audit(
    profile: pd.DataFrame,
    isolated: pd.Series,
    *,
    redshift: float,
    h: float,
    scatter_dex: float,
    tau_bounds: tuple[float, float],
    chi2_threshold: float = 2.30,
) -> pd.DataFrame:
    isolated_mass = float(isolated["M_delta_msun_isolated_extrapolation"])
    isolated_concentration = float(isolated["c_delta_isolated_extrapolation"])
    isolated_median = float(dutton_maccio_c200(isolated_mass, redshift, h=h))
    isolated_offset = float(
        np.log10(isolated_concentration / isolated_median) / scatter_dex
    )
    accepted = profile[profile["profile_mass_chi2"] <= chi2_threshold]
    if accepted.empty:
        threshold = profile.loc[profile["profile_mass_chi2"].idxmin()]
        threshold_reached = False
        interpolated_offset = float(threshold["concentration_offset_sigma"])
    else:
        ordered = profile.sort_values("concentration_offset_sigma").reset_index(drop=True)
        accepted_indices = ordered.index[
            ordered["profile_mass_chi2"] <= chi2_threshold
        ].to_numpy()
        first_index = int(accepted_indices[0])
        threshold = ordered.iloc[first_index]
        interpolated_offset = float(threshold["concentration_offset_sigma"])
        if first_index > 0:
            previous = ordered.iloc[first_index - 1]
            x0 = float(previous["concentration_offset_sigma"])
            x1 = float(threshold["concentration_offset_sigma"])
            y0 = float(previous["profile_mass_chi2"])
            y1 = float(threshold["profile_mass_chi2"])
            if (y0 - chi2_threshold) * (y1 - chi2_threshold) <= 0 and y1 != y0:
                interpolated_offset = x0 + (
                    (chi2_threshold - y0) * (x1 - x0) / (y1 - y0)
                )
        threshold_reached = True
    penalized = _best_penalized_candidate(profile)
    return pd.DataFrame([{
        "redshift": redshift,
        "h": h,
        "concentration_scatter_dex": scatter_dex,
        "tau_min": tau_bounds[0],
        "tau_max": tau_bounds[1],
        "profile_chi2_threshold": chi2_threshold,
        "isolated_M200_msun": isolated_mass,
        "isolated_c200": isolated_concentration,
        "isolated_c200_median": isolated_median,
        "isolated_concentration_offset_sigma": isolated_offset,
        "threshold_reached": threshold_reached,
        "minimum_offset_sigma_at_threshold_grid": float(
            threshold["concentration_offset_sigma"]
        ),
        "minimum_offset_sigma_at_threshold_interpolated": interpolated_offset,
        "threshold_profile_chi2": float(threshold["profile_mass_chi2"]),
        "threshold_M200_msun": float(threshold["M200_msun"]),
        "threshold_c200": float(threshold["c200_candidate"]),
        "threshold_tau_rt_over_rs": float(threshold["tau_rt_over_rs"]),
        "threshold_r_t_pc": float(threshold["r_t_pc"]),
        "penalized_offset_sigma": float(penalized["concentration_offset_sigma"]),
        "penalized_profile_chi2": float(penalized["profile_mass_chi2"]),
        "penalized_score_chi2_plus_z2": float(
            penalized["profile_mass_chi2"]
            + penalized["concentration_offset_sigma"] ** 2
        ),
        "penalized_M200_msun": float(penalized["M200_msun"]),
        "penalized_c200": float(penalized["c200_candidate"]),
        "penalized_tau_rt_over_rs": float(penalized["tau_rt_over_rs"]),
    }])


def plot_prior_audit(
    profile: pd.DataFrame,
    summary: pd.DataFrame,
    output: str,
    sensitivity_profile: pd.DataFrame | None = None,
    sensitivity_summary: pd.DataFrame | None = None,
) -> None:
    row = summary.iloc[0]
    redshift = float(row["redshift"])
    h = float(row["h"])
    scatter = float(row["concentration_scatter_dex"])
    masses = np.logspace(5.5, 10.0, 240)
    median = np.asarray(dutton_maccio_c200(masses, redshift, h=h))

    fig, axes = plt.subplots(1, 2, figsize=(9.8, 3.9), constrained_layout=True)
    ax = axes[0]
    ax.fill_between(
        masses, median * 10 ** (-3 * scatter), median * 10 ** (3 * scatter),
        color="#D9E4F0", label=r"$\pm3\sigma_{\log c}$",
    )
    ax.fill_between(
        masses, median * 10 ** (-scatter), median * 10 ** scatter,
        color="#76A5CF", label=r"$\pm1\sigma_{\log c}$",
    )
    ax.plot(masses, median, color="#173F5F", linewidth=1.8, label="median")
    ax.scatter(
        row["isolated_M200_msun"], row["isolated_c200"], marker="*", s=95,
        color="#C43C39", edgecolor="black", linewidth=0.5,
        label="isolated extrapolation", zorder=5,
    )
    ax.scatter(
        row["threshold_M200_msun"], row["threshold_c200"], marker="D", s=48,
        color="#E09F3E", edgecolor="black", linewidth=0.5,
        label="tidal-envelope fit", zorder=5,
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$M_{200}$ [$M_\odot$]")
    ax.set_ylabel(r"$c_{200}$")
    ax.grid(which="both", alpha=0.2, linewidth=0.6)
    ax.legend(frameon=False, fontsize=7.7, loc="upper right")
    ax.text(0.03, 0.96, "(a)", transform=ax.transAxes, va="top", fontweight="bold")

    ax = axes[1]
    x = profile["concentration_offset_sigma"].to_numpy(dtype=float)
    y = profile["profile_mass_chi2"].to_numpy(dtype=float)
    y_plot = np.clip(y, 1e-5, None)
    ax.plot(
        x, y_plot, color="#2A9D8F", marker="o", markersize=3.2,
        linewidth=1.5, label=r"baseline $r_t/r_s\geq1$",
    )
    if sensitivity_profile is not None:
        sensitivity_x = sensitivity_profile[
            "concentration_offset_sigma"
        ].to_numpy(dtype=float)
        sensitivity_y = np.clip(
            sensitivity_profile["profile_mass_chi2"].to_numpy(dtype=float),
            1e-5,
            None,
        )
        ax.plot(
            sensitivity_x, sensitivity_y, color="#E09F3E", linestyle="--",
            linewidth=1.4, label=r"aggressive $r_t/r_s\geq0.5$",
        )
    ax.axhline(
        float(row["profile_chi2_threshold"]), color="black", linestyle="--",
        linewidth=1.1, label=r"$\chi^2_{2M}=2.30$",
    )
    ax.axvline(
        float(row["minimum_offset_sigma_at_threshold_interpolated"]), color="#C43C39",
        linestyle=":", linewidth=1.4,
        label=rf"baseline minimum $\simeq{row['minimum_offset_sigma_at_threshold_interpolated']:.2g}\sigma$",
    )
    if sensitivity_summary is not None:
        sensitivity_row = sensitivity_summary.iloc[0]
        ax.axvline(
            float(sensitivity_row["minimum_offset_sigma_at_threshold_interpolated"]),
            color="#E09F3E", linestyle="-.", linewidth=1.25,
            label=rf"aggressive minimum $\simeq{sensitivity_row['minimum_offset_sigma_at_threshold_interpolated']:.2g}\sigma$",
        )
    ax.set_yscale("log")
    ax.set_ylim(1e-5, 3e2)
    ax.set_xlabel(r"Concentration offset $\Delta\log_{10}c/(0.11\,{\rm dex})$")
    ax.set_ylabel(r"Profiled two-mass $\chi^2$")
    ax.grid(which="both", alpha=0.2, linewidth=0.6)
    ax.legend(
        frameon=True, framealpha=0.92, facecolor="white", edgecolor="none",
        fontsize=7.7, loc="upper right",
    )
    ax.text(0.03, 0.96, "(b)", transform=ax.transAxes, va="top", fontweight="bold")

    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    if path.suffix.lower() != ".pdf":
        fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--isolated-csv", required=True)
    parser.add_argument("--profile-output", required=True)
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--figure", required=True)
    parser.add_argument("--sensitivity-profile", default=None)
    parser.add_argument("--sensitivity-summary", default=None)
    parser.add_argument("--redshift", type=float, default=0.881)
    parser.add_argument("--h0", type=float, default=67.4)
    parser.add_argument("--omega-m", type=float, default=0.315)
    parser.add_argument("--scatter-dex", type=float, default=DUTTON_SCATTER_DEX)
    parser.add_argument("--tau-min", type=float, default=1.0)
    parser.add_argument("--tau-max", type=float, default=100.0)
    parser.add_argument("--offset-min", type=float, default=0.0)
    parser.add_argument("--offset-max", type=float, default=11.0)
    parser.add_argument("--offset-step", type=float, default=0.25)
    args = parser.parse_args()
    if args.h0 <= 0 or not 0 < args.omega_m < 1:
        parser.error("cosmology must be physical")
    if args.tau_min <= 0 or args.tau_max <= args.tau_min:
        parser.error("tau bounds must be positive and increasing")
    if args.offset_step <= 0 or args.offset_max <= args.offset_min:
        parser.error("offset grid must be increasing")

    h = args.h0 / 100.0
    cosmology = FlatLambdaCDM(H0=args.h0, Om0=args.omega_m, Tcmb0=2.7255)
    rho_critical = float(
        cosmology.critical_density(args.redshift).to_value(u.Msun / u.pc**3)
    )
    offsets = np.arange(
        args.offset_min, args.offset_max + 0.5 * args.offset_step, args.offset_step
    )
    log_tau_bounds = (np.log10(args.tau_min), np.log10(args.tau_max))
    profile = profile_grid(
        offsets,
        args.redshift,
        rho_critical,
        h=h,
        scatter_dex=args.scatter_dex,
        log_tau_bounds=log_tau_bounds,
    )
    isolated = pd.read_csv(args.isolated_csv).iloc[0]
    summary = summarize_prior_audit(
        profile,
        isolated,
        redshift=args.redshift,
        h=h,
        scatter_dex=args.scatter_dex,
        tau_bounds=(args.tau_min, args.tau_max),
    )

    profile_path = Path(args.profile_output)
    summary_path = Path(args.summary_output)
    profile_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    profile.to_csv(profile_path, index=False, quoting=csv.QUOTE_MINIMAL)
    summary.to_csv(summary_path, index=False, quoting=csv.QUOTE_MINIMAL)
    sensitivity_profile = (
        pd.read_csv(args.sensitivity_profile) if args.sensitivity_profile else None
    )
    sensitivity_summary = (
        pd.read_csv(args.sensitivity_summary) if args.sensitivity_summary else None
    )
    if (sensitivity_profile is None) != (sensitivity_summary is None):
        parser.error("sensitivity profile and summary must be supplied together")
    plot_prior_audit(
        profile, summary, args.figure, sensitivity_profile, sensitivity_summary
    )
    print(f"wrote {len(profile)} concentration-profile rows to {profile_path}")
    print(summary.iloc[0].to_string())
    print(f"wrote {args.figure}")


if __name__ == "__main__":
    main()
