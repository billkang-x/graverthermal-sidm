"""Build a Born-controlled cooling-failure map for the physical B1938 halo.

The dense map is a post-processing calculation on the common, unevolved
192-shell physical halo.  At fixed threshold and ``eta_B``, the microscopic
cooling source scales linearly with the requested reference cross section, so
one exact thermal-kernel evaluation per threshold is sufficient.  Existing
direct simulations are overplotted as validation points; the map itself is
not presented as an additional gravothermal evolution scan.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "..", ".."))
sys.path[:0] = [
    os.path.join(_ROOT, "src", "cross_sections"),
    os.path.join(_ROOT, "src", "fluid_runner"),
]

from audit_physical_activity import (  # noqa: E402
    COSMIC_AGE_AT_LENS_GYR,
    local_timescale_diagnostics,
)
from born_valid_scan import parameters_for_target_sigma_at_threshold  # noqa: E402
from dsidm_models import (  # noqa: E402
    benchmark_models,
    born_expansion_parameter,
    sigma_T_born,
)
from run_born_valid import build_halo  # noqa: E402


MODEL_KEYS = ("M1_dark_photon_massive", "M2_scalar_phi_massive")
MODEL_LABELS = {
    "M1_dark_photon_massive": "M1: massive-vector emission",
    "M2_scalar_phi_massive": "M2: massive-scalar emission",
}


def _read_many(paths: list[str]) -> pd.DataFrame:
    if not paths:
        raise ValueError("at least one input CSV is required")
    return pd.concat([pd.read_csv(path) for path in paths], ignore_index=True)


def _reference_directories(frame: pd.DataFrame) -> dict[str, str]:
    required = {
        "model", "requested_source_time_gyr", "direct_output_dir", "status",
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"time summaries lack columns: {', '.join(sorted(missing))}")
    initial = frame[
        (frame["status"] == "complete")
        & np.isclose(frame["requested_source_time_gyr"].astype(float), 0.0)
    ]
    result: dict[str, str] = {}
    for model_key in MODEL_KEYS:
        rows = initial[initial["model"].astype(str) == model_key]
        if rows.empty:
            raise ValueError(f"no complete t=0 direct state for {model_key}")
        direct_dir = str(rows.iloc[0]["direct_output_dir"])
        if not Path(direct_dir).is_dir():
            raise FileNotFoundError(f"reference direct directory not found: {direct_dir}")
        result[model_key] = direct_dir
    return result


def _threshold_reference(
    model_key: str,
    eta: float,
    threshold_kms: float,
    reference_dir: str,
    n_shells: int,
) -> tuple[dict[str, float | int], float]:
    model, halo, sigma_reference = build_halo(
        model_key,
        eta,
        reference_dir,
        n_shells,
        1.0,
        True,
        threshold_kms,
        0.0,
        0.0,
    )
    timescales = local_timescale_diagnostics(
        np.asarray(halo.u, dtype=float),
        np.asarray(getattr(halo, "C_cool", np.zeros_like(halo.u)), dtype=float),
        np.asarray(halo.rho, dtype=float),
        float(halo.scale_t.to("Gyr").value),
    )
    if not np.isclose(sigma_reference, 1.0, rtol=2e-10, atol=0.0):
        raise RuntimeError(f"reference cross section failed to close: {sigma_reference}")
    return timescales, float(born_expansion_parameter(model))


def build_map(
    time_frame: pd.DataFrame,
    *,
    eta: float = 0.1,
    n_shells: int = 192,
    ratio_grid: np.ndarray | None = None,
    sigma_grid: np.ndarray | None = None,
    vmax_kms: float | None = None,
) -> pd.DataFrame:
    """Return the dense failure map as one tidy machine-readable table."""
    if ratio_grid is None:
        ratio_grid = np.logspace(-1.0, 2.0, 49)
    if sigma_grid is None:
        sigma_grid = np.logspace(-1.0, 1.0, 41)
    ratio_grid = np.asarray(ratio_grid, dtype=float)
    sigma_grid = np.asarray(sigma_grid, dtype=float)
    if np.any(~np.isfinite(ratio_grid)) or np.any(ratio_grid <= 0):
        raise ValueError("ratio_grid must contain finite positive values")
    if np.any(~np.isfinite(sigma_grid)) or np.any(sigma_grid <= 0):
        raise ValueError("sigma_grid must contain finite positive values")

    if vmax_kms is None:
        velocities = pd.to_numeric(time_frame["direct_vmax_kms"], errors="coerce")
        velocities = velocities[np.isfinite(velocities)]
        if velocities.empty:
            raise ValueError("time summaries do not contain a finite direct_vmax_kms")
        vmax_kms = float(np.median(velocities))
    if not np.isfinite(vmax_kms) or vmax_kms <= 0:
        raise ValueError("vmax_kms must be finite and positive")

    references = _reference_directories(time_frame)
    base_models = benchmark_models()
    rows: list[dict] = []
    for model_key in MODEL_KEYS:
        for ratio in ratio_grid:
            threshold = float(ratio * vmax_kms)
            timescales, eta_actual = _threshold_reference(
                model_key, eta, threshold, references[model_key], n_shells,
            )
            for sigma_target in sigma_grid:
                model = parameters_for_target_sigma_at_threshold(
                    base_models[model_key], eta, float(sigma_target), threshold, 100.0
                )
                sigma_actual = float(sigma_T_born(np.array([100.0]), model)[0])
                cooling_time_gyr = float(timescales["cooling_time_gyr"])
                if np.isfinite(cooling_time_gyr):
                    cooling_time_gyr /= float(sigma_target)
                    cooling_to_dynamical = (
                        float(timescales["cooling_to_dynamical"])
                        / float(sigma_target)
                    )
                    cooling_to_age = cooling_time_gyr / COSMIC_AGE_AT_LENS_GYR
                else:
                    cooling_to_dynamical = float("inf")
                    cooling_to_age = float("inf")
                rows.append({
                    "model": model_key,
                    "eta_B": eta_actual,
                    "vmax_phys_kms": vmax_kms,
                    "vstar_over_vmax": float(ratio),
                    "vstar_kms": threshold,
                    "sigma_m_100_target": float(sigma_target),
                    "sigma_m_100_kms": sigma_actual,
                    "min_cooling_time_gyr_initial_state": cooling_time_gyr,
                    "local_dynamical_time_gyr_initial_state": float(
                        timescales["dynamical_time_gyr"]
                    ),
                    "cooling_to_dynamical_time_initial_state": cooling_to_dynamical,
                    "cosmic_age_at_lens_gyr": COSMIC_AGE_AT_LENS_GYR,
                    "cooling_to_cosmic_age_initial_state": cooling_to_age,
                    "cooling_min_shell_index_initial_state": int(
                        timescales["shell_index"]
                    ),
                    "m_chi_GeV": float(model.m_chi),
                    "m_mediator_keV": float(model.m_mediator * 1e6),
                    "born_controlled": bool(eta_actual <= 0.1 + 1e-12),
                })
    return pd.DataFrame(rows)


def plot_map(map_frame: pd.DataFrame, audit_frame: pd.DataFrame, output: str) -> None:
    """Plot M1/M2 failure margins and overlay the direct 192-shell points."""
    finite = map_frame[
        np.isfinite(map_frame["cooling_to_dynamical_time_initial_state"])
        & (map_frame["cooling_to_dynamical_time_initial_state"] > 0)
    ]
    if finite.empty:
        raise ValueError("failure map contains no finite cooling-time ratios")
    log_values = np.log10(finite["cooling_to_dynamical_time_initial_state"])
    vmin = float(np.floor(log_values.min()))
    vmax = float(np.ceil(log_values.max()))
    if vmax <= vmin:
        vmax = vmin + 1.0

    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.0), constrained_layout=True)
    image = None
    for ax, model_key, panel in zip(axes, MODEL_KEYS, ("(a)", "(b)")):
        group = map_frame[map_frame["model"] == model_key].copy()
        if "sigma_m_100_target" not in group:
            # Older generated tables can be replotted without repeating the
            # expensive kernel integrations.  Target values are logarithmic
            # grid coordinates; the actual values differ only at roundoff.
            group["sigma_m_100_target"] = 10.0 ** np.round(
                np.log10(group["sigma_m_100_kms"].to_numpy(dtype=float)), 12
            )
        x = np.sort(group["vstar_over_vmax"].unique())
        y = np.sort(group["sigma_m_100_target"].unique())
        values = group.pivot(
            index="sigma_m_100_target",
            columns="vstar_over_vmax",
            values="cooling_to_dynamical_time_initial_state",
        ).reindex(index=y, columns=x).to_numpy(dtype=float)
        log_ratio = np.log10(np.clip(values, 10**vmin, 10**vmax))
        log_ratio[~np.isfinite(log_ratio)] = vmax
        image = ax.pcolormesh(
            x, y, log_ratio, shading="auto", cmap="viridis", vmin=vmin, vmax=vmax,
        )
        ax.axvline(1.0, color="white", linestyle="--", linewidth=1.25)
        ax.text(
            0.93, 0.12, "radiatively open", color="white", fontsize=8.0,
            ha="right", va="bottom", rotation=90,
        )

        direct = audit_frame[
            (audit_frame["status"] == "complete")
            & (audit_frame["model"].astype(str) == model_key)
        ].copy()
        direct_x = 1.0 / direct["threshold_ratio_vmax_over_vstar"].to_numpy(float)
        ax.scatter(
            direct_x,
            direct["sigma_m_100_kms"],
            marker="o",
            s=29,
            facecolors="none",
            edgecolors="white",
            linewidths=0.9,
            label="192-shell candidates",
        )
        minimum = group["cooling_to_dynamical_time_initial_state"].replace(
            [np.inf, -np.inf], np.nan
        ).min()
        ax.text(
            0.03, 0.05,
            rf"minimum $t_{{\rm cool}}/t_{{\rm dyn}}={minimum:.1e}$",
            transform=ax.transAxes, color="white", fontsize=8.0,
        )
        ax.text(
            0.03, 0.96, panel, transform=ax.transAxes, color="white",
            va="top", fontweight="bold",
        )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$v_\star/v_{\max}$")
        ax.set_title(MODEL_LABELS[model_key], fontsize=10.0)
        ax.grid(which="both", color="white", alpha=0.12, linewidth=0.5)
    axes[0].set_ylabel(r"$\sigma_T/m_\chi(100\,{\rm km\,s^{-1}})$ "
                       r"[cm$^2$ g$^{-1}$]")
    axes[1].legend(frameon=False, fontsize=8.0, loc="upper right", labelcolor="white")
    assert image is not None
    colorbar = fig.colorbar(image, ax=axes, pad=0.02, aspect=28)
    colorbar.set_label(r"$\log_{10}(t_{\rm cool}/t_{\rm dyn})$")

    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    if path.suffix.lower() != ".pdf":
        fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--time-csv", nargs="+", required=True)
    parser.add_argument("--audit-csv", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--figure", required=True)
    parser.add_argument("--eta", type=float, default=0.1)
    parser.add_argument("--n-shells", type=int, default=192)
    parser.add_argument("--vmax-kms", type=float, default=None)
    parser.add_argument("--ratio-min", type=float, default=0.1)
    parser.add_argument("--ratio-max", type=float, default=100.0)
    parser.add_argument("--ratio-count", type=int, default=49)
    parser.add_argument("--sigma-min", type=float, default=0.1)
    parser.add_argument("--sigma-max", type=float, default=10.0)
    parser.add_argument("--sigma-count", type=int, default=41)
    args = parser.parse_args()
    if args.ratio_count < 2 or args.sigma_count < 2:
        parser.error("grid counts must be at least two")
    if args.ratio_min <= 0 or args.ratio_max <= args.ratio_min:
        parser.error("ratio bounds must be positive and increasing")
    if args.sigma_min <= 0 or args.sigma_max <= args.sigma_min:
        parser.error("sigma bounds must be positive and increasing")

    time_frame = _read_many(args.time_csv)
    audit_frame = _read_many(args.audit_csv)
    result = build_map(
        time_frame,
        eta=args.eta,
        n_shells=args.n_shells,
        ratio_grid=np.logspace(
            np.log10(args.ratio_min), np.log10(args.ratio_max), args.ratio_count
        ),
        sigma_grid=np.logspace(
            np.log10(args.sigma_min), np.log10(args.sigma_max), args.sigma_count
        ),
        vmax_kms=args.vmax_kms,
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, index=False, quoting=csv.QUOTE_MINIMAL)
    plot_map(result, audit_frame, args.figure)
    finite = result["cooling_to_dynamical_time_initial_state"].replace(
        [np.inf, -np.inf], np.nan
    )
    print(f"wrote {len(result)} map rows to {output}")
    print(f"minimum t_cool/t_dyn = {finite.min():.6g}")
    print(f"wrote {args.figure}")


if __name__ == "__main__":
    main()
