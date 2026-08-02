"""Plot the extended 192-shell threshold and cross-section audits."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from audit_physical_activity import parse_point_id


THRESHOLD_COLORS = {
    1.0: "#CC79A7",
    5.0: "#0072B2",
    20.0: "#009E73",
    100.0: "#E69F00",
    500.0: "#D55E00",
}
SIGMA_COLORS = {0.1: "#009E73", 1.0: "#0072B2", 10.0: "#D55E00"}
MARKERS = {1.0: "D", 5.0: "o", 20.0: "s", 100.0: "^", 500.0: "v"}


def _read_many(paths: list[str]) -> pd.DataFrame:
    if not paths:
        raise ValueError("at least one input CSV is required")
    return pd.concat([pd.read_csv(path) for path in paths], ignore_index=True)


def _with_coordinates(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    coordinates = result["scan_point_id"].map(parse_point_id)
    result[["eta_B_scan", "sigma_scan", "vstar_scan"]] = pd.DataFrame(
        coordinates.tolist(), index=result.index
    )
    return result


def _save(fig: plt.Figure, output: str) -> None:
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    if path.suffix.lower() != ".pdf":
        fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_time_scans(time_frame: pd.DataFrame, output: str) -> None:
    frame = _with_coordinates(time_frame)
    frame = frame.loc[frame["status"] == "complete"].copy()
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.8), constrained_layout=True)

    ax = axes[0]
    selected = frame[np.isclose(frame["sigma_scan"], 0.1)]
    for threshold in sorted(selected["vstar_scan"].unique()):
        curve = selected[np.isclose(selected["vstar_scan"], threshold)].groupby(
            "requested_source_time_gyr", as_index=False
        )["direct_mass_chi2"].mean()
        ax.plot(
            curve["requested_source_time_gyr"], curve["direct_mass_chi2"],
            color=THRESHOLD_COLORS[float(threshold)], marker=MARKERS[float(threshold)],
            linewidth=1.7, markersize=4.8,
            label=rf"$v_\star={threshold:g}$",
        )
    ax.set_xlabel(r"Source snapshot time $t_{\rm source}$ [Gyr]")
    ax.set_ylabel(r"Direct two-mass $\chi^2$")
    ax.set_xlim(-0.003, 0.103)
    ax.grid(alpha=0.22, linewidth=0.7)
    ax.legend(frameon=False, fontsize=8.0, ncol=2, title=r"km s$^{-1}$")
    ax.text(
        0.98, 0.95, "(a)", transform=ax.transAxes, ha="right", va="top",
        fontweight="bold",
    )

    ax = axes[1]
    selected = frame[np.isclose(frame["vstar_scan"], 5.0)]
    for sigma in sorted(selected["sigma_scan"].unique()):
        curve = selected[np.isclose(selected["sigma_scan"], sigma)].groupby(
            "requested_source_time_gyr", as_index=False
        )["direct_mass_chi2"].mean()
        ax.plot(
            curve["requested_source_time_gyr"], curve["direct_mass_chi2"],
            color=SIGMA_COLORS[float(sigma)], marker="o", linewidth=1.7,
            markersize=4.8, label=rf"${sigma:g}$",
        )
    ax.set_xlabel(r"Source snapshot time $t_{\rm source}$ [Gyr]")
    ax.set_ylabel(r"Direct two-mass $\chi^2$")
    ax.set_xlim(-0.003, 0.103)
    ax.grid(alpha=0.22, linewidth=0.7)
    ax.legend(
        frameon=False, fontsize=8.0,
        title=r"$\sigma_T/m_\chi(100)$ [cm$^2$ g$^{-1}$]",
    )
    ax.text(0.02, 0.95, "(b)", transform=ax.transAxes, va="top", fontweight="bold")

    _save(fig, output)


def plot_activity(audit_frame: pd.DataFrame, output: str) -> None:
    frame = _with_coordinates(audit_frame)
    frame = frame.loc[frame["status"] == "complete"].copy()
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.8), constrained_layout=True)

    ax = axes[0]
    ratios = frame.groupby("vstar_scan", as_index=False)[
        "threshold_ratio_vmax_over_vstar"
    ].mean().sort_values("vstar_scan")
    x = np.arange(len(ratios))
    values = ratios["threshold_ratio_vmax_over_vstar"].to_numpy()
    bars = ax.bar(
        x, values, width=0.64,
        color=[THRESHOLD_COLORS[float(value)] for value in ratios["vstar_scan"]],
        edgecolor="black", linewidth=0.6,
    )
    ax.axhline(1.0, color="black", linestyle="--", linewidth=1.1)
    ax.set_yscale("log")
    ax.set_ylim(7e-3, 8.0)
    ax.set_xticks(x, [rf"${value:g}$" for value in ratios["vstar_scan"]])
    ax.set_xlabel(r"Threshold velocity $v_\star$ [km s$^{-1}$]")
    ax.set_ylabel(r"$v_{\max,\rm phys}/v_\star$")
    ax.grid(axis="y", which="both", alpha=0.22, linewidth=0.7)
    for bar, value in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2, value * 1.14, f"{value:.3g}",
            ha="center", va="bottom", fontsize=8.0,
        )
    ax.text(0.02, 0.95, "(a)", transform=ax.transAxes, va="top", fontweight="bold")

    ax = axes[1]
    finite = frame[np.isfinite(frame["min_cooling_time_gyr_initial_state"])].copy()
    finite = finite[finite["vstar_scan"].isin([1.0, 5.0, 20.0])]
    for threshold in sorted(finite["vstar_scan"].unique()):
        for model, linestyle, marker in (("M1", "-", "o"), ("M2", "--", "s")):
            group = finite[
                np.isclose(finite["vstar_scan"], threshold)
                & finite["model"].astype(str).str.startswith(model)
            ].sort_values("sigma_scan")
            ax.plot(
                group["sigma_scan"], group["min_cooling_time_gyr_initial_state"],
                color=THRESHOLD_COLORS[float(threshold)], linestyle=linestyle,
                marker=marker, linewidth=1.6, markersize=4.8,
            )
    ax.axhline(6.37, color="black", linestyle=":", linewidth=1.2)
    ax.text(0.12, 8.2, "6.37 Gyr", fontsize=8.0, ha="left", va="bottom")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.075, 14.0)
    ax.set_xlabel(r"$\sigma_T/m_\chi(100)$ [cm$^2$ g$^{-1}$]")
    ax.set_ylabel(r"Minimum local cooling time [Gyr]")
    ax.grid(which="both", alpha=0.22, linewidth=0.7)
    threshold_handles = [
        Line2D([0], [0], color=THRESHOLD_COLORS[value], lw=1.8,
               label=rf"$v_\star={value:g}$")
        for value in (1.0, 5.0, 20.0)
    ]
    model_handles = [
        Line2D([0], [0], color="black", linestyle="-", marker="o", label="M1"),
        Line2D([0], [0], color="black", linestyle="--", marker="s", label="M2"),
    ]
    first = ax.legend(handles=threshold_handles, frameon=False, fontsize=7.8,
                      loc="lower left")
    ax.add_artist(first)
    ax.legend(handles=model_handles, frameon=False, fontsize=7.8, loc="upper right")
    ax.text(0.02, 0.95, "(b)", transform=ax.transAxes, va="top", fontweight="bold")

    _save(fig, output)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--time-csv", nargs="+", required=True)
    parser.add_argument("--audit-csv", nargs="+", required=True)
    parser.add_argument("--time-output", required=True)
    parser.add_argument("--activity-output", required=True)
    args = parser.parse_args()
    plot_time_scans(_read_many(args.time_csv), args.time_output)
    plot_activity(_read_many(args.audit_csv), args.activity_output)
    print(f"wrote {args.time_output} and {args.activity_output}")


if __name__ == "__main__":
    main()
