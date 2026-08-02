"""Plot the final 192-shell B1938 time scan and physical activity audit."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


COLORS = {20: "#0072B2", 100: "#009E73", 500: "#D55E00"}
MARKERS = {20: "o", 100: "s", 500: "^"}


def _threshold_from_id(value: str) -> int:
    return int(str(value).rsplit("vstar", 1)[1])


def make_figure(time_csv: str, audit_csv: str, output: str) -> None:
    time_frame = pd.read_csv(time_csv)
    audit = pd.read_csv(audit_csv)
    time_frame["vstar_kms"] = time_frame["scan_point_id"].map(_threshold_from_id)

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.8), constrained_layout=True)

    ax = axes[0]
    for threshold in (20, 100, 500):
        group = time_frame[time_frame["vstar_kms"] == threshold]
        curve = group.groupby("requested_source_time_gyr", as_index=False).agg(
            chi2=("direct_mass_chi2", "mean"),
            spread=("direct_mass_chi2", lambda values: float(np.ptp(values))),
        )
        ax.plot(
            curve["requested_source_time_gyr"], curve["chi2"],
            color=COLORS[threshold], marker=MARKERS[threshold], linewidth=1.8,
            markersize=5.2, label=rf"$v_\star={threshold}\,\mathrm{{km\,s^{{-1}}}}$",
        )
    ax.set_xlabel(r"Source snapshot time $t_{\rm source}$ [Gyr]")
    ax.set_ylabel(r"Direct two-mass $\chi^2$")
    ax.set_xlim(-0.003, 0.103)
    ax.grid(alpha=0.22, linewidth=0.7)
    ax.legend(frameon=False, fontsize=8.5)
    ax.text(0.02, 0.95, "(a)", transform=ax.transAxes, va="top", fontweight="bold")

    ax = axes[1]
    ratios = audit.groupby("vstar_kms", as_index=False)[
        "threshold_ratio_vmax_over_vstar"
    ].mean().sort_values("vstar_kms")
    x = np.arange(len(ratios))
    values = ratios["threshold_ratio_vmax_over_vstar"].to_numpy()
    bars = ax.bar(
        x, values, width=0.62,
        color=[COLORS[int(value)] for value in ratios["vstar_kms"]],
        edgecolor="black", linewidth=0.6,
    )
    ax.axhline(1.0, color="black", linestyle="--", linewidth=1.1,
               label="Emission threshold")
    ax.set_yscale("log")
    ax.set_ylim(7e-3, 1.8)
    ax.set_xticks(x, [rf"${int(value)}$" for value in ratios["vstar_kms"]])
    ax.set_xlabel(r"Threshold velocity $v_\star$ [km s$^{-1}$]")
    ax.set_ylabel(r"$v_{\max,\rm phys}/v_\star$")
    ax.grid(axis="y", which="both", alpha=0.22, linewidth=0.7)
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")
    for bar, value in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2, value * 1.13, f"{value:.3f}",
            ha="center", va="bottom", fontsize=8.5,
        )
    ax.text(0.02, 0.95, "(b)", transform=ax.transAxes, va="top", fontweight="bold")

    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    if path.suffix.lower() != ".pdf":
        fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--time-csv", required=True)
    parser.add_argument("--audit-csv", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    make_figure(args.time_csv, args.audit_csv, args.output)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
