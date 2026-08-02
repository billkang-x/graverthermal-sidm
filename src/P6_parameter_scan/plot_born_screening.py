"""Plot the Born-valid mass-hierarchy and reference-cross-section screening."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize


MODEL_LABELS = {
    "M1_dark_photon_massive": "M1 vector emission",
    "M2_scalar_phi_massive": "M2 scalar emission",
}
MODEL_MARKERS = {
    "M1_dark_photon_massive": "o",
    "M2_scalar_phi_massive": "s",
}
MODEL_COLORS = {
    "M1_dark_photon_massive": "#0072B2",
    "M2_scalar_phi_massive": "#D55E00",
}


def make_figure(input_csv: str, output: str) -> None:
    frame = pd.read_csv(input_csv)
    frame = frame[frame["model"].isin(MODEL_LABELS)].copy()
    frame["model_label"] = frame["model"].map(MODEL_LABELS)

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.8), constrained_layout=True)
    norm = Normalize(vmin=0.1, vmax=1.0)
    cmap = plt.get_cmap("viridis")

    ax = axes[0]
    for model, group in frame.groupby("model", sort=False):
        color = MODEL_COLORS[model]
        marker = MODEL_MARKERS[model]
        for overlap, points in group.groupby("astrophysical_overlap"):
            if overlap:
                ax.scatter(
                    points["m_chi_GeV"], points["m_mediator_keV"],
                    c=points["eta_target"], cmap=cmap, norm=norm,
                    marker=marker, s=34, linewidths=0.65,
                    edgecolors="black", alpha=0.9, label=None,
                )
            else:
                ax.scatter(
                    points["m_chi_GeV"], points["m_mediator_keV"],
                    marker=marker, s=34, linewidths=0.65,
                    edgecolors=color, facecolors="none", alpha=0.9,
                    label=None,
                )
        ax.scatter([], [], marker=marker, s=34, color=color,
                   label=MODEL_LABELS[model])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Dark-matter mass $m_\chi$ [GeV]")
    ax.set_ylabel(r"Mediator mass $m_{\rm med}$ [keV]")
    ax.grid(which="both", alpha=0.2, linewidth=0.65)
    ax.legend(frameon=False, fontsize=8.2, loc="upper left")
    ax.text(0.03, 0.95, "(a)", transform=ax.transAxes,
            va="top", fontweight="bold")
    ax.text(0.97, 0.04, "black edge: overlap", transform=ax.transAxes,
            ha="right", fontsize=7.5)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array(np.array([]))
    cbar = fig.colorbar(sm, ax=ax, pad=0.02, fraction=0.046)
    cbar.set_label(r"Hierarchy parameter $\eta$")

    ax = axes[1]
    for model, group in frame.groupby("model", sort=False):
        marker = MODEL_MARKERS[model]
        color = MODEL_COLORS[model]
        for overlap, points in group.groupby("astrophysical_overlap"):
            if overlap:
                ax.scatter(
                    points["v_star_km_s"], points["sigma_m_cm2_g"],
                    c=points["eta_target"], cmap=cmap, norm=norm,
                    marker=marker, s=34, linewidths=0.65,
                    edgecolors="black", alpha=0.9,
                )
            else:
                ax.scatter(
                    points["v_star_km_s"], points["sigma_m_cm2_g"],
                    marker=marker, s=34, linewidths=0.65,
                    edgecolors=color, facecolors="none", alpha=0.9,
                )
    ax.axhspan(0.1, 10.0, color="#009E73", alpha=0.10,
               label=r"reference scan: $0.1$--$10$")
    ax.set_yscale("log")
    ax.set_xlabel(r"Threshold velocity $v_\star$ [km s$^{-1}$]")
    ax.set_ylabel(r"$\sigma_T/m_\chi$ at $v_\star$ [cm$^2$ g$^{-1}$]")
    ax.set_xlim(170, 530)
    ax.grid(which="both", alpha=0.2, linewidth=0.65)
    ax.legend(frameon=False, fontsize=8.2, loc="lower right")
    ax.text(0.03, 0.95, "(b)", transform=ax.transAxes,
            va="top", fontweight="bold")

    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    if path.suffix.lower() != ".pdf":
        fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-csv", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    make_figure(args.input_csv, args.output)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
