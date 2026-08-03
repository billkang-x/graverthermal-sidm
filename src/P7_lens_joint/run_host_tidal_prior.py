"""Run the B1938 host-orbit/tidal sensitivity analysis."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from host_tidal_prior import (
    ORBIT_PRIORS,
    FreePseudoJaffeSummary,
    PublishedTidalReference,
    importance_sample_orbit_posterior,
    minimum_current_radius_pc,
    weighted_quantile,
)


COLORS = {
    "circular_upper_envelope": "#1B4F72",
    "phase_mixed": "#2A9D8F",
    "radial_sensitivity": "#C43C39",
}
LABELS = {
    "circular_upper_envelope": "circular upper envelope",
    "phase_mixed": r"phase mixed: $f_p\sim\mathrm{Beta}(4,2)$",
    "radial_sensitivity": r"radial sensitivity: $f_p\sim\mathrm{Beta}(2,4)$",
}


def posterior_summary(result: dict, reference, imaging) -> dict[str, float | str]:
    weight = np.asarray(result["weight"])
    current = np.asarray(result["current_radius_pc"])
    line_of_sight = np.abs(np.asarray(result["line_of_sight_pc"]))
    pericentre = np.asarray(result["pericentre_pc"])
    predicted_tidal = np.asarray(result["predicted_tidal_radius_pc"])
    quantiles = (0.16, 0.50, 0.84)
    current_q = weighted_quantile(current, quantiles, weight)
    los_q = weighted_quantile(line_of_sight, quantiles, weight)
    pericentre_q = weighted_quantile(pericentre, quantiles, weight)
    tidal_q = weighted_quantile(predicted_tidal, quantiles, weight)
    return {
        "orbit_prior": result["orbit_prior"],
        "sample_count": len(weight),
        "effective_sample_size": result["effective_sample_size"],
        "host_scale_radius_pc": result["host_scale_radius_pc"],
        "host_max_radius_pc": result["host_max_radius_pc"],
        "projected_radius_pc": reference.projected_radius_pc,
        "current_radius_p16_pc": current_q[0],
        "current_radius_median_pc": current_q[1],
        "current_radius_p84_pc": current_q[2],
        "abs_los_p16_pc": los_q[0],
        "abs_los_median_pc": los_q[1],
        "abs_los_p84_pc": los_q[2],
        "pericentre_p16_pc": pericentre_q[0],
        "pericentre_median_pc": pericentre_q[1],
        "pericentre_p84_pc": pericentre_q[2],
        "predicted_rt_p16_pc": tidal_q[0],
        "predicted_rt_median_pc": tidal_q[1],
        "predicted_rt_p84_pc": tidal_q[2],
        "posterior_probability_current_r_lt_10kpc": float(
            np.sum(weight[current < 10_000.0])
        ),
        "free_pj_rt_mean_pc": imaging.tidal_radius_pc,
        "free_pj_rt_sigma_pc": imaging.tidal_radius_sigma_pc,
    }


def plot_results(results: list[dict], reference, imaging, output: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(11.1, 3.6), constrained_layout=True)
    bins = np.linspace(reference.projected_radius_pc / 1000.0, 60.0, 100)
    for result in results:
        name = str(result["orbit_prior"])
        axes[0].hist(
            np.asarray(result["current_radius_pc"]) / 1000.0,
            bins=bins, weights=np.asarray(result["weight"]), histtype="step",
            linewidth=1.7, color=COLORS[name], label=LABELS[name], density=True,
        )
        axes[1].hist(
            np.asarray(result["pericentre_pc"]) / 1000.0,
            bins=bins, weights=np.asarray(result["weight"]), histtype="step",
            linewidth=1.7, color=COLORS[name], density=True,
        )
        axes[2].hist(
            np.asarray(result["predicted_tidal_radius_pc"]),
            bins=np.linspace(40.0, 210.0, 100),
            weights=np.asarray(result["weight"]), histtype="step",
            linewidth=1.7, color=COLORS[name], density=True,
        )
    lower_bound = float(minimum_current_radius_pc(
        imaging.tidal_radius_pc, imaging.total_mass_msun, reference=reference
    ))
    axes[0].axvline(
        lower_bound / 1000.0, color="black", linestyle="--", linewidth=1.1,
        label="mean-parameter lower bound",
    )
    axes[0].set_xlabel(r"current 3D radius [kpc]")
    axes[0].set_ylabel("posterior density")
    axes[0].set_xlim(0.0, 60.0)
    axes[0].legend(frameon=False, fontsize=7.2)
    axes[0].set_title("(a) deprojected location", fontsize=9.5)
    axes[1].set_xlabel(r"pericentre [kpc]")
    axes[1].set_ylabel("posterior density")
    axes[1].set_xlim(0.0, 60.0)
    axes[1].set_title("(b) tidal-setting orbit", fontsize=9.5)
    axes[2].axvspan(
        imaging.tidal_radius_pc - imaging.tidal_radius_sigma_pc,
        imaging.tidal_radius_pc + imaging.tidal_radius_sigma_pc,
        color="#D9D9D9", alpha=0.8, label=r"free-PJ $1\sigma$",
    )
    axes[2].axvline(
        reference.tidal_radius_pc, color="#E09F3E", linestyle=":",
        linewidth=1.4, label="projected-radius tide",
    )
    axes[2].set_xlabel(r"tidal radius [pc]")
    axes[2].set_ylabel("posterior density")
    axes[2].set_xlim(40.0, 210.0)
    axes[2].set_title("(c) imaging-conditioned tide", fontsize=9.5)
    axes[2].legend(frameon=False, fontsize=7.2)
    for axis in axes:
        axis.grid(alpha=0.18, linewidth=0.6)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--figure", required=True)
    parser.add_argument("--samples", type=int, default=500_000)
    parser.add_argument("--seed", type=int, default=1938666)
    parser.add_argument("--posterior-draws", type=int, default=20_000)
    parser.add_argument(
        "--host-scale-radii-kpc", default="10,30,100",
        help="NFW-like tracer scale radii for the phase-mixed sensitivity table",
    )
    args = parser.parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    reference = PublishedTidalReference()
    imaging = FreePseudoJaffeSummary()
    results = []
    summaries = []
    rng = np.random.default_rng(args.seed + 101)
    posterior_rows = []
    for index, orbit_prior in enumerate(ORBIT_PRIORS):
        result = importance_sample_orbit_posterior(
            args.samples, orbit_prior, seed=args.seed + index
        )
        results.append(result)
        summaries.append(posterior_summary(result, reference, imaging))
        chosen = rng.choice(
            args.samples, size=args.posterior_draws, replace=True,
            p=np.asarray(result["weight"]),
        )
        posterior_rows.append(pd.DataFrame({
            "orbit_prior": orbit_prior,
            "line_of_sight_pc": np.asarray(result["line_of_sight_pc"])[chosen],
            "current_radius_pc": np.asarray(result["current_radius_pc"])[chosen],
            "pericentre_fraction": np.asarray(result["pericentre_fraction"])[chosen],
            "pericentre_pc": np.asarray(result["pericentre_pc"])[chosen],
            "total_mass_msun": np.asarray(result["total_mass_msun"])[chosen],
            "predicted_tidal_radius_pc": np.asarray(
                result["predicted_tidal_radius_pc"]
            )[chosen],
        }))
    summary = pd.DataFrame(summaries)
    summary.to_csv(output_dir / "host_tidal_summary.csv", index=False)
    pd.concat(posterior_rows, ignore_index=True).to_csv(
        output_dir / "host_tidal_posterior_draws.csv", index=False
    )
    scale_radii_kpc = [
        float(value) for value in args.host_scale_radii_kpc.split(",") if value.strip()
    ]
    if not scale_radii_kpc or any(value <= 0 for value in scale_radii_kpc):
        parser.error("host scale radii must be a comma-separated positive list")
    scale_summaries = []
    for index, scale_kpc in enumerate(scale_radii_kpc):
        result = importance_sample_orbit_posterior(
            max(args.samples // 2, 100_000), "phase_mixed",
            seed=args.seed + 1000 + index,
            host_scale_radius_pc=scale_kpc * 1000.0,
        )
        scale_summaries.append(posterior_summary(result, reference, imaging))
    pd.DataFrame(scale_summaries).to_csv(
        output_dir / "host_scale_sensitivity.csv", index=False
    )
    lower_bound = float(minimum_current_radius_pc(
        imaging.tidal_radius_pc, imaging.total_mass_msun, reference=reference
    ))
    provenance = {
        "reference": reference.to_dict(),
        "free_pseudo_jaffe_summary": asdict(imaging),
        "orbit_independent_current_radius_lower_bound_pc": lower_bound,
        "orbit_independent_abs_los_lower_bound_pc": float(np.sqrt(
            max(lower_bound**2 - reference.projected_radius_pc**2, 0.0)
        )),
        "important_limitation": (
            "The public free-PJ mass and radius summaries are treated as "
            "independent Gaussians because the PRONTO posterior chain is private."
        ),
        "host_scale_sensitivity_kpc": scale_radii_kpc,
    }
    with (output_dir / "host_tidal_metadata.json").open("w", encoding="ascii") as f:
        json.dump(provenance, f, indent=2, sort_keys=True)
    plot_results(results, reference, imaging, Path(args.figure))
    print(summary.to_string(index=False))
    print(f"orbit-independent current-radius lower bound: {lower_bound:.1f} pc")


if __name__ == "__main__":
    main()
