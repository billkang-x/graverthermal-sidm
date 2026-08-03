"""Generate a reproducible data-quality audit for the public GM068 UVFITS."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from uvfits_data import audit_visibility_noise
from visibility_likelihood import classify_noise_bins


ARCHIVE_DOI = "https://doi.org/10.48717/wch4-m437"
ARCHIVE_URL = (
    "https://archive.jive.eu/exp/GM068_111106/pipe/"
    "gm068_B1938+666.UVDATA.FITS"
)
PUBLISHER_URL = "https://doi.org/10.1038/s41550-025-02651-2"
SENSITIVE_TRIANGLE = {"EF-JB", "EF-WB", "JB-WB"}


def sha256(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(chunk_size):
            digest.update(block)
    return digest.hexdigest()


def plot_audit(noise, baselines, uv_sample, output: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.7), constrained_layout=True)

    valid_uv = uv_sample[uv_sample["stokes_i_weight"] > 0]
    axes[0].scatter(
        valid_uv["u_mlambda"], valid_uv["v_mlambda"], s=0.25,
        color="#1B4F72", alpha=0.25, rasterized=True,
    )
    axes[0].scatter(
        -valid_uv["u_mlambda"], -valid_uv["v_mlambda"], s=0.25,
        color="#1B4F72", alpha=0.25, rasterized=True,
    )
    axes[0].set_xlabel(r"$u$ [M$\lambda$]")
    axes[0].set_ylabel(r"$v$ [M$\lambda$]")
    axes[0].set_aspect("equal", adjustable="box")
    axes[0].set_title("(a) sampled $uv$ coverage", fontsize=9.5)

    plotted = noise[np.isfinite(noise["adjacent_sigma_jy"])]
    axes[1].scatter(
        plotted["start_hour"], plotted["adjacent_sigma_jy"] * 1e3,
        s=3.0, color="#2A9D8F", alpha=0.35, rasterized=True,
    )
    median_by_time = plotted.groupby("time_bin", as_index=False).agg(
        start_hour=("start_hour", "first"),
        adjacent_sigma_jy=("adjacent_sigma_jy", "median"),
    )
    axes[1].plot(
        median_by_time["start_hour"], median_by_time["adjacent_sigma_jy"] * 1e3,
        color="#C43C39", linewidth=1.4, label="baseline median",
    )
    axes[1].set_yscale("log")
    axes[1].set_xlabel("hours from first visibility")
    axes[1].set_ylabel("adjacent-difference RMS [mJy]")
    axes[1].set_title("(b) 30-minute noise audit", fontsize=9.5)
    axes[1].legend(frameon=False, fontsize=7.5)

    ordered = baselines.sort_values("median_adjacent_sigma_jy")
    colors = [
        "#C43C39" if label in SENSITIVE_TRIANGLE else "#E09F3E"
        for label in ordered["baseline"]
    ]
    axes[2].scatter(
        ordered["rms_pipeline_sigma_jy"] * 1e3,
        ordered["median_adjacent_sigma_jy"] * 1e3,
        c=colors, s=15, edgecolor="black", linewidth=0.25,
    )
    finite = ordered[["rms_pipeline_sigma_jy", "median_adjacent_sigma_jy"]].to_numpy()
    finite = finite[np.isfinite(finite).all(axis=1) & (finite > 0).all(axis=1)]
    if finite.size:
        lower = float(np.min(finite) * 1e3 * 0.8)
        upper = float(np.max(finite) * 1e3 * 1.25)
        axes[2].plot([lower, upper], [lower, upper], "--", color="black", linewidth=0.8)
        axes[2].set_xlim(lower, upper)
        axes[2].set_ylim(lower, upper)
    axes[2].set_xscale("log")
    axes[2].set_yscale("log")
    axes[2].set_xlabel("pipeline-weight RMS [mJy]")
    axes[2].set_ylabel("adjacent-difference RMS [mJy]")
    axes[2].set_title("(c) pipeline-weight comparison", fontsize=9.5)

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight", dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("uvfits")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--figure", required=True)
    parser.add_argument("--chunk-groups", type=int, default=250_000)
    parser.add_argument("--bin-minutes", type=float, default=30.0)
    parser.add_argument("--max-pair-seconds", type=float, default=10.0)
    parser.add_argument("--uv-sample-stride", type=int, default=500)
    args = parser.parse_args()

    input_path = Path(args.uvfits)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    noise, baselines, uv_sample, summary = audit_visibility_noise(
        input_path,
        chunk_groups=args.chunk_groups,
        bin_minutes=args.bin_minutes,
        max_pair_separation_seconds=args.max_pair_seconds,
        uv_sample_stride=args.uv_sample_stride,
    )
    noise = classify_noise_bins(noise)
    noise.to_csv(output_dir / "noise_bins.csv", index=False)
    baselines.to_csv(output_dir / "baseline_summary.csv", index=False)
    uv_sample.to_csv(output_dir / "uv_sample.csv", index=False)
    summary.update({
        "sha256": sha256(input_path),
        "archive_doi": ARCHIVE_DOI,
        "archive_url": ARCHIVE_URL,
        "publisher_url": PUBLISHER_URL,
        "sensitive_triangle_removed_by_published_analysis": sorted(SENSITIVE_TRIANGLE),
    })
    with (output_dir / "metadata.json").open("w", encoding="ascii") as stream:
        json.dump(summary, stream, indent=2, sort_keys=True)
    plot_audit(noise, baselines, uv_sample, Path(args.figure))
    print(f"valid visibilities: {summary['valid_complex_visibility_count']:,}")
    print(f"valid fraction: {summary['valid_fraction']:.6f}")
    print(f"noise bins: {len(noise):,}; baselines: {len(baselines):,}")


if __name__ == "__main__":
    main()
