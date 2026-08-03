"""Make a full-data naturally weighted GM068 dirty image for coordinate QA."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import finufft
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits

from lens_forward_model import ARCSEC_TO_RADIAN, ImageGrid
from uvfits_data import read_uvfits_metadata
from visibility_likelihood import VisibilityNoiseModel


def grid_dirty_image(
    uvfits_path: str | Path,
    noise_csv: str | Path,
    metadata_json: str | Path,
    image_grid: ImageGrid,
    *,
    chunk_groups: int = 250_000,
    tolerance: float = 1.0e-6,
) -> tuple[np.ndarray, dict]:
    """Adjoint-grid all retained Stokes-I samples without deconvolution."""
    if chunk_groups <= 0 or tolerance <= 0:
        raise ValueError("chunk size and NUFFT tolerance must be positive")
    with Path(metadata_json).open(encoding="ascii") as stream:
        audit_metadata = json.load(stream)
    noise_model = VisibilityNoiseModel.from_csv(
        noise_csv,
        start_jd=float(audit_metadata["start_jd"]),
        bin_minutes=float(audit_metadata["noise_bin_minutes"]),
    )
    metadata = read_uvfits_metadata(uvfits_path)
    frequencies = np.asarray(metadata.frequencies_hz, dtype=float)
    pixel_radians = image_grid.pixel_scale_arcsec * ARCSEC_TO_RADIAN
    center_x_radians = image_grid.center_x_arcsec * ARCSEC_TO_RADIAN
    center_y_radians = image_grid.center_y_arcsec * ARCSEC_TO_RADIAN
    accumulated = np.zeros((image_grid.nx, image_grid.ny), dtype=np.complex128)
    total_weight = 0.0
    gridded_samples = 0
    retained_groups = 0

    with fits.open(uvfits_path, memmap=True) as hdus:
        groups = hdus[0].data
        for start in range(0, metadata.group_count, chunk_groups):
            stop = min(start + chunk_groups, metadata.group_count)
            times = np.asarray(groups.par("DATE")[start:stop], dtype=float)
            baselines = np.rint(groups.par("BASELINE")[start:stop]).astype(np.int64)
            uu_seconds = np.asarray(groups.par("UU---SIN")[start:stop], dtype=float)
            vv_seconds = np.asarray(groups.par("VV---SIN")[start:stop], dtype=float)
            sigma, retained = noise_model.lookup(baselines, times)
            retained_groups += int(np.count_nonzero(retained))
            raw = np.asarray(groups.data[start:stop])[:, 0, 0, :, 0, :, :]
            real = raw[..., 0].astype(float, copy=False)
            imag = raw[..., 1].astype(float, copy=False)
            pipeline_weight = raw[..., 2].astype(float, copy=False)
            valid = (
                np.isfinite(real) & np.isfinite(imag)
                & np.isfinite(pipeline_weight) & (pipeline_weight > 0)
                & retained[:, None, None]
            )
            hand_count = np.sum(valid, axis=-1)
            visibility_sum = np.sum(
                np.where(valid, real + 1j * imag, 0.0), axis=-1
            )
            stokes_i = np.zeros(hand_count.shape, dtype=np.complex128)
            np.divide(
                visibility_sum, hand_count, out=stokes_i, where=hand_count > 0
            )
            inverse_variance = hand_count / sigma[:, None] ** 2
            selected = hand_count > 0
            if not np.any(selected):
                continue
            frequency_grid = np.broadcast_to(
                frequencies[None, :], hand_count.shape
            )[selected]
            uu_grid = np.broadcast_to(
                uu_seconds[:, None], hand_count.shape
            )[selected]
            vv_grid = np.broadcast_to(
                vv_seconds[:, None], hand_count.shape
            )[selected]
            weight = inverse_variance[selected]
            values = stokes_i[selected]
            x = 2.0 * np.pi * uu_grid * frequency_grid * pixel_radians
            y = 2.0 * np.pi * vv_grid * frequency_grid * pixel_radians
            phase_conjugate = np.exp(
                2j * np.pi * frequency_grid
                * (uu_grid * center_x_radians + vv_grid * center_y_radians)
            )
            accumulated += finufft.nufft2d1(
                x, y, values * weight * phase_conjugate,
                (image_grid.nx, image_grid.ny), isign=1, eps=tolerance,
            )
            total_weight += float(np.sum(weight))
            gridded_samples += int(np.count_nonzero(selected))
            print(
                f"gridded groups {start:,}:{stop:,}; "
                f"Stokes-I samples={gridded_samples:,}",
                flush=True,
            )
    if total_weight <= 0:
        raise ValueError("no visibility samples survived the likelihood mask")
    dirty = np.real(accumulated.T) / total_weight
    summary = {
        "uvfits": str(Path(uvfits_path).resolve()),
        "noise_csv": str(Path(noise_csv).resolve()),
        "nx": image_grid.nx,
        "ny": image_grid.ny,
        "pixel_scale_mas": image_grid.pixel_scale_arcsec * 1000.0,
        "center_x_arcsec": image_grid.center_x_arcsec,
        "center_y_arcsec": image_grid.center_y_arcsec,
        "retained_group_count": retained_groups,
        "gridded_stokes_i_count": gridded_samples,
        "natural_weight_sum": total_weight,
        "dirty_peak_jy_per_beam": float(np.max(dirty)),
        "dirty_minimum_jy_per_beam": float(np.min(dirty)),
        "purpose": "coordinate and phase-centre QA only; not deconvolved",
        "model_coordinate_registration": (
            "not published; PRONTO x/y values must not be overplotted directly"
        ),
    }
    return dirty, summary


def write_fits(image: np.ndarray, grid: ImageGrid, output: Path) -> None:
    header = fits.Header()
    header["BUNIT"] = "Jy/beam"
    header["CTYPE1"] = "RA---SIN"
    header["CTYPE2"] = "DEC--SIN"
    header["CUNIT1"] = "deg"
    header["CUNIT2"] = "deg"
    header["CRPIX1"] = grid.nx // 2 + 1
    header["CRPIX2"] = grid.ny // 2 + 1
    header["CRVAL1"] = grid.center_x_arcsec / 3600.0
    header["CRVAL2"] = grid.center_y_arcsec / 3600.0
    header["CDELT1"] = -grid.pixel_scale_arcsec / 3600.0
    header["CDELT2"] = grid.pixel_scale_arcsec / 3600.0
    header["HISTORY"] = "Natural-weight dirty image for GM068 likelihood QA."
    output.parent.mkdir(parents=True, exist_ok=True)
    fits.PrimaryHDU(np.asarray(image, dtype=np.float32), header=header).writeto(
        output, overwrite=True
    )


def plot_dirty_image(image: np.ndarray, grid: ImageGrid, output: Path) -> None:
    xx, yy = grid.coordinates()
    extent = [xx.min(), xx.max(), yy.min(), yy.max()]
    scale = float(np.quantile(np.abs(image), 0.9995)) * 1e3
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.0), constrained_layout=True)
    for axis in axes:
        displayed = axis.imshow(
            image * 1e3, origin="lower", extent=extent, cmap="RdBu_r",
            vmin=-scale, vmax=scale, interpolation="nearest",
        )
        axis.set_xlabel(r"$x$ [arcsec]")
        axis.set_ylabel(r"$y$ [arcsec]")
        axis.set_aspect("equal")
    axes[0].set_title("(a) full naturally weighted dirty image", fontsize=9.5)
    axes[1].set_xlim(-0.08, 0.08)
    axes[1].set_ylim(-0.08, 0.08)
    axes[1].set_title("(b) UVFITS phase-centre region", fontsize=9.5)
    colorbar = fig.colorbar(displayed, ax=axes, shrink=0.88, pad=0.02)
    colorbar.set_label("dirty intensity [mJy beam$^{-1}$]")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("uvfits")
    parser.add_argument("--noise-csv", required=True)
    parser.add_argument("--metadata-json", required=True)
    parser.add_argument("--fits-output", required=True)
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--figure", required=True)
    parser.add_argument("--pixels", type=int, default=768)
    parser.add_argument("--pixel-scale-mas", type=float, default=1.5)
    parser.add_argument("--center-x", type=float, default=0.12)
    parser.add_argument("--center-y", type=float, default=0.12)
    parser.add_argument("--chunk-groups", type=int, default=250_000)
    args = parser.parse_args()
    grid = ImageGrid(
        args.pixels, args.pixels, args.pixel_scale_mas / 1000.0,
        args.center_x, args.center_y,
    )
    image, summary = grid_dirty_image(
        args.uvfits, args.noise_csv, args.metadata_json, grid,
        chunk_groups=args.chunk_groups,
    )
    write_fits(image, grid, Path(args.fits_output))
    summary_path = Path(args.summary_output)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w", encoding="ascii") as stream:
        json.dump(summary, stream, indent=2, sort_keys=True)
    plot_dirty_image(image, grid, Path(args.figure))
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
