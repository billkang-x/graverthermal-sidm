"""
P4 exclusion plot generator.

Produces the main multi-model exclusion figure:
  - X-axis: σ_T/m at low velocity (dwarf scale, v ~ 50 km/s)
  - Y-axis: σ_T/m at high velocity (cluster scale, v ~ 1000-3000 km/s)
  - Each model is a CURVE (parameterized by alpha_D) on this plane
  - Velocity-independent models are DIAGONAL lines (σ/m_low = σ/m_high)
  - Observational constraints are SHADED REGIONS (allowed vs excluded)

Key insight: velocity-dependent models trace distinct curves in this plane,
while velocity-independent models lie on the diagonal. The slope of each
model curve is the "smoking gun" signature.

Output: figures/P4_exclusion.png
"""
from __future__ import annotations

import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJ_ROOT = os.path.normpath(os.path.join(_HERE, '..', '..'))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_PROJ_ROOT, 'src', 'cross_sections'))
sys.path.insert(0, os.path.join(_PROJ_ROOT, 'src', 'fluid_runner'))

from constraints import CONSTRAINTS, get_constraint, ObservationalConstraint
from dsidm_models import benchmark_models, sigma_T_born, r_diss as r_diss_func, DSIDMParameters
from model_eval import compute_effective_sigma_m, scan_alpha_for_model, CONSTRAINT_T

OUT_FIG = os.path.join(_PROJ_ROOT, 'figures')
OUT_DATA = os.path.join(_PROJ_ROOT, 'data')
os.makedirs(OUT_FIG, exist_ok=True)
os.makedirs(OUT_DATA, exist_ok=True)


# ----------------------------------------------------------------------
# Plot configuration: (v_low, v_high) axes
# ----------------------------------------------------------------------
V_LOW_LABEL = 'dwarf'        # characteristic v ~ 50 km/s
V_HIGH_LABEL = 'cluster'     # characteristic v ~ 1000-3000 km/s

# Indices into CONSTRAINT_T / model_eval arrays
CONSTRAINT_IDX = list(CONSTRAINT_T.keys())  # ['bullet', 'cluster_cores', 'dwarf_cores', 'b1938']
IDX_DWARF = CONSTRAINT_IDX.index('dwarf_cores')
IDX_CLUSTER = CONSTRAINT_IDX.index('cluster_cores')
IDX_BULLET = CONSTRAINT_IDX.index('bullet')
IDX_B1938 = CONSTRAINT_IDX.index('b1938')


# ----------------------------------------------------------------------
# Constraint regions in (σ_low, σ_high) plane
# ----------------------------------------------------------------------
def get_excluded_regions():
    """Return dict describing allowed vs excluded regions in the
    (σ/m at dwarf, σ/m at cluster) plane.

    For each constraint we get a vertical or horizontal band depending on
    whether the constraint applies to σ_low (dwarf) or σ_high (cluster).
    """
    dwarf = get_constraint('dwarf_cores')
    cluster = get_constraint('cluster_cores')
    bullet = get_constraint('bullet')
    b1938 = get_constraint('b1938')
    return {
        'dwarf_upper':  ('v', dwarf.sigma_upper),  # exclude σ_low > 50
        'dwarf_lower':  ('v', dwarf.sigma_lower),  # require σ_low > 1 (allowed region)
        'cluster_upper':('h', cluster.sigma_upper),
        'bullet_upper': ('h', bullet.sigma_upper),
        'b1938_range':  ('v_band', (b1938.sigma_lower, b1938.sigma_upper)),
    }


# ----------------------------------------------------------------------
# Plot main exclusion figure
# ----------------------------------------------------------------------
def plot_exclusion(out_path):
    """Build the multi-model exclusion plot."""
    fig, ax = plt.subplots(figsize=(10, 8))

    # ----- 1. Constraint shading -----
    # Allowed region: σ_low in [1, 50], σ_high in [0, 1] (for cluster)
    # We shade excluded regions lightly.
    sigma_low_max = 1e4    # plot bounds
    sigma_low_min = 1e-3
    sigma_high_max = 1e3
    sigma_high_min = 1e-4

    # Shaded: σ_low > 50 (dwarfs excluded) — light red band on right
    ax.axvspan(50, sigma_low_max, color='red', alpha=0.07)
    ax.axvline(50, color='teal', linestyle='--', lw=1.2, alpha=0.7)
    ax.text(70, 0.01, 'Dwarfs\nexcluded\n($\\sigma/m > 50$)', fontsize=9,
            color='teal', ha='left', va='bottom')

    # Shaded: σ_low < 1 (no cores in dwarfs) — left band
    ax.axvspan(sigma_low_min, 1.0, color='gray', alpha=0.07)
    ax.axvline(1.0, color='teal', linestyle='--', lw=1.2, alpha=0.7)
    ax.text(0.5, 0.01, 'Cuspy\ndwarfs\n(no cores)', fontsize=9,
            color='teal', ha='center', va='bottom')

    # σ_high > 1.0 (cluster cores excluded) — top band
    ax.axhspan(1.0, sigma_high_max, color='purple', alpha=0.06)
    ax.axhline(1.0, color='purple', linestyle='--', lw=1.2, alpha=0.7)
    ax.text(0.005, 1.5, 'Cluster cores excluded ($\\sigma/m_{\\rm cl} > 1$)',
            fontsize=9, color='purple', ha='left', va='bottom')

    # σ_high > 1.25 (Bullet) — darker top band
    ax.axhspan(1.25, sigma_high_max, color='darkred', alpha=0.08)
    ax.axhline(1.25, color='darkred', linestyle='--', lw=1.2, alpha=0.7)
    ax.text(0.005, 2.0, 'Bullet Cluster ($\\sigma/m > 1.25$)',
            fontsize=9, color='darkred', ha='left', va='bottom')

    # B1938 band: SELF-CONSISTENT viable region from P5 grid scan.
    # The grid scan directly simulated (sigma/m, r_s) with the DissipativeHalo
    # fluid code and found the viable sigma/m range where t_cross <= t_obs.
    # This replaces the previous elastic-approximation band [0.07, 7.8].
    #
    # The grid scan sigma/m is at v_ref=100 km/s; the P4 x-axis is at v=50 km/s.
    # Conversion: sigma(v=50) = sigma(v=100) * (sigma(50)/sigma(100)).
    # For M1/M2 (velocity-dependent), this ratio is ~14-15 (steep rise at low v).
    # For M3 (constant), the ratio is 1.0.
    #
    # Self-consistent viable sigma/m at v=50 km/s (dwarf scale):
    #   M1: [0.073, 14.5] cm²/g
    #   M2: [0.070, 14.0] cm²/g
    #   M3: [0.011, 1.0]  cm²/g (constant)
    #
    # We shade the union of viable regions and annotate per-model boundaries.
    b1938_viable_ranges = {
        'M1': (0.073, 14.5),
        'M2': (0.070, 14.0),
        'M3': (0.011, 1.0),
    }
    # Union for shading
    all_lo = min(r[0] for r in b1938_viable_ranges.values())
    all_hi = max(r[1] for r in b1938_viable_ranges.values())
    ax.axvspan(all_lo, all_hi, color='orange', alpha=0.06)
    # Per-model boundaries (dashed)
    model_colors = {'M1': 'tab:green', 'M2': 'tab:red', 'M3': 'tab:gray'}
    for ml, (lo, hi) in b1938_viable_ranges.items():
        ax.axvline(lo, color=model_colors[ml], linestyle=':', lw=1.0, alpha=0.7)
        ax.axvline(hi, color=model_colors[ml], linestyle=':', lw=1.0, alpha=0.7)
    ax.text(0.5, 0.0003, 'B1938+666\n(self-consistent\ngrid scan)',
            fontsize=9, color='darkorange', ha='center', va='bottom',
            style='italic')

    # ----- 2. Model curves -----
    bm = benchmark_models()
    model_styles = {
        'M1_dark_photon_massive':   {'color': 'tab:green', 'marker': 'o', 'label': r'M1: dark photon ($m_V=1.1$ keV)'},
        'M2_scalar_phi_massive':    {'color': 'tab:red',   'marker': 's', 'label': r'M2: scalar $\phi$ ($m_\phi=6.95$ keV)'},
        'M3_massless_control':      {'color': 'tab:gray',  'marker': '^', 'label': r'M3: massless control'},
    }

    T_arr = np.array([CONSTRAINT_T[k] for k in CONSTRAINT_IDX])

    # Scan alpha_D and trace each model's curve in (σ_low, σ_high) plane
    all_curves = []
    for mk, p in bm.items():
        if mk not in model_styles:
            continue
        result = scan_alpha_for_model(mk, p, alpha_log_range=(-3, 3), n_alpha=80)
        sigma_low = result['sigma_eff'][:, IDX_DWARF]
        sigma_high_cluster = result['sigma_eff'][:, IDX_CLUSTER]
        sigma_high_bullet = result['sigma_eff'][:, IDX_BULLET]
        # Mask NaN
        mask = np.isfinite(sigma_low) & np.isfinite(sigma_high_cluster)
        sigma_low = sigma_low[mask]
        sigma_high_cluster = sigma_high_cluster[mask]
        sigma_high_bullet = sigma_high_bullet[mask]
        # Sort by sigma_low for clean line plot
        order = np.argsort(sigma_low)
        sigma_low = sigma_low[order]
        sigma_high_cluster = sigma_high_cluster[order]
        sigma_high_bullet = sigma_high_bullet[order]

        st = model_styles[mk]
        ax.plot(sigma_low, sigma_high_cluster, color=st['color'], lw=2.0,
                label=st['label'] + r' ($v_{\rm high}=1000$ km/s)', alpha=0.85)
        ax.plot(sigma_low, sigma_high_bullet, color=st['color'], lw=1.5,
                linestyle=':', alpha=0.6,
                label=st['label'] + r' ($v_{\rm high}=3000$ km/s)')
        # Scatter points along the curve (selected alpha values)
        idx_pts = np.linspace(0, len(sigma_low)-1, 8).astype(int)
        ax.scatter(sigma_low[idx_pts], sigma_high_cluster[idx_pts],
                   color=st['color'], marker=st['marker'], s=30, zorder=5)

        all_curves.append({
            'model': mk, 'sigma_low': sigma_low,
            'sigma_high_cluster': sigma_high_cluster,
            'sigma_high_bullet': sigma_high_bullet,
        })

    # ----- 3. Velocity-independent reference (diagonal) -----
    # σ_low = σ_high = const; we draw a few reference iso-σ contours
    sig_iso_arr = np.logspace(-2, 2, 5)
    for s_iso in sig_iso_arr:
        ax.plot([s_iso, s_iso], [sigma_high_min, sigma_high_max],
                color='lightgray', lw=0.6, zorder=1)
        ax.plot([sigma_low_min, sigma_low_max], [s_iso, s_iso],
                color='lightgray', lw=0.6, zorder=1)

    # Velocity-independent model = diagonal line
    sig_diag = np.logspace(-3, 3, 100)
    ax.plot(sig_diag, sig_diag, color='black', lw=1.5, linestyle='--',
            label=r'Velocity-independent ($\sigma_{\rm low}=\sigma_{\rm high}$)',
            alpha=0.6)

    # ----- 4. Decorations -----
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(sigma_low_min, sigma_low_max)
    ax.set_ylim(sigma_high_min, sigma_high_max)
    ax.set_xlabel(r'$\sigma_T/m_\chi$ at dwarf scale ($v\sim 50$ km/s) '
                  '[cm$^2$/g]', fontsize=12)
    ax.set_ylabel(r'$\sigma_T/m_\chi$ at cluster scale ($v\sim 1000$ km/s) '
                  '[cm$^2$/g]', fontsize=12)
    ax.set_title('SIDM model exclusion: velocity-dependent models trace '
                 'distinct curves\n'
                 r'(each curve = family parameterized by $\alpha_D$)',
                 fontsize=12)

    # Custom legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc='upper left', fontsize=9, framealpha=0.9,
              ncol=1)

    ax.grid(True, which='both', alpha=0.25, zorder=0)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    print(f"Saved {out_path}")
    return all_curves


# ----------------------------------------------------------------------
# Side figure: σ/m vs velocity for each model (illustrates velocity dependence)
# ----------------------------------------------------------------------
def plot_sigma_vs_v(out_path, alpha_rescale=1.0):
    """Plot σ_T/m as a function of velocity for each model, with constraint
    bands overlaid."""
    fig, ax = plt.subplots(figsize=(9, 6.5))

    v_arr = np.logspace(0.5, 4.0, 200)  # ~3 to 10000 km/s
    bm = benchmark_models()
    model_styles = {
        'M1_dark_photon_massive':   {'color': 'tab:green', 'label': r'M1: dark photon ($m_V=1.1$ keV, $v^*\!\approx\!200$ km/s)'},
        'M2_scalar_phi_massive':    {'color': 'tab:red',   'label': r'M2: scalar $\phi$ ($m_\phi=6.95$ keV, $v^*\!\approx\!500$ km/s)'},
        'M3_massless_control':      {'color': 'tab:gray',  'label': r'M3: massless control'},
    }

    # Rescale alpha so σ/m at v=100 km/s matches ~5 cm²/g for visibility
    target_sigma_at_100 = 5.0
    for mk, p in bm.items():
        if mk not in model_styles:
            continue
        # Find alpha that gives target_sigma_at_100
        sig_at_100 = sigma_T_born(np.array([100.0]), p)[0]
        alpha_factor = np.sqrt(np.sqrt(target_sigma_at_100 / max(sig_at_100, 1e-30)))
        # alpha_D scales as g^2; σ ∝ alpha_D^2 ∝ g^4
        p_rescaled = DSIDMParameters(
            model=p.model, m_chi=p.m_chi, m_mediator=p.m_mediator,
            alpha_D=p.alpha_D * alpha_factor**2,
            m_mediator_heavy=p.m_mediator_heavy,
            mediation=p.mediation, emission_type=p.emission_type,
        )
        sig = sigma_T_born(v_arr, p_rescaled)
        st = model_styles[mk]
        ax.plot(v_arr, sig, color=st['color'], lw=2.0, label=st['label'])

    # Constraint bands
    for c in CONSTRAINTS:
        if c.sigma_upper is not None:
            ax.axhline(c.sigma_upper, color=c.color, linestyle=c.linestyle,
                       lw=1.3, alpha=0.8)
            ax.text(1.2, c.sigma_upper * 1.05, f'{c.label} (upper: {c.sigma_upper})',
                    color=c.color, fontsize=9, va='bottom', ha='left')
        if c.sigma_lower is not None:
            ax.axhline(c.sigma_lower, color=c.color, linestyle=':',
                       lw=1.3, alpha=0.6)
            ax.text(1.2, c.sigma_lower * 1.05, f'{c.label} (lower: {c.sigma_lower})',
                    color=c.color, fontsize=9, va='bottom', ha='left')

    # Characteristic velocity markers
    for c in CONSTRAINTS:
        ax.axvline(c.v_char_kms, color=c.color, linestyle=':', lw=1.0, alpha=0.4)
        ax.text(c.v_char_kms * 1.05, 0.02, c.label,
                color=c.color, fontsize=9, rotation=90, va='bottom', ha='left')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(3, 10000)
    ax.set_ylim(0.01, 1e4)
    ax.set_xlabel(r'relative velocity $v$ [km/s]', fontsize=12)
    ax.set_ylabel(r'$\sigma_T/m_\chi$ [cm$^2$/g]', fontsize=12)
    ax.set_title(r'Velocity-dependent cross sections: each model has '
                 r'distinct $\sigma(v)$ slope', fontsize=12)
    ax.legend(loc='upper right', fontsize=9, framealpha=0.9)
    ax.grid(True, which='both', alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    print(f"Saved {out_path}")


# ----------------------------------------------------------------------
# Supplementary figure: r_diss dimension (addresses review #5)
# ----------------------------------------------------------------------
# The main exclusion plot (P4_exclusion.png) lives in the (sigma_low, sigma_high)
# plane, which captures only the elastic-scattering dimension. Dissipative SIDM
# has a second, independent dimension: r_diss(v) (the fractional kinetic-energy
# loss per collision). Two models can share the same sigma_T(v) yet have very
# different r_diss(v), and therefore different halo thermodynamics.
#
# This supplementary figure has TWO panels to make the physics explicit:
#   LEFT  : r_diss at B1938/dwarf velocity scale (~80 km/s) — shows that
#           massive emission is KINEMATICALLY SUPPRESSED at low v (r_diss ≈ 1
#           for M1, M2; only M3 has finite r_diss). This is the symmetry-
#           breaking suppression documented in P3.
#   RIGHT : r_diss at cluster velocity scale (~1200 km/s) — shows MODEL
#           DISCRIMINATION: M1 (vector, v*=200) is now well above threshold and
#           has r_diss > 1; M2 (scalar, v*=500) is also above threshold but
#           suppressed by C_phi; M3 remains constant.
#
# Both panels share the X-axis (sigma/m at dwarf scale), so a model family
# (parameterized by alpha_D) traces a horizontal line in each panel.

V_RDISS_LOW_KEY = 'b1938'        # ~80 km/s  (low-velocity panel)
V_RDISS_HIGH_KEY = 'cluster_cores'  # ~1200 km/s (high-velocity panel)


def plot_rdiss_dimension(out_path):
    """Build the 2-panel supplementary (sigma_dwarf, r_diss_eff) figure.

    Layout:
    - Legend placed in the upper-right corner of the right panel, where
      the data region is empty (r_diss values cluster near 1.00-1.05 at
      the bottom of the plot).
    - Panel titles carry the velocity-scale info; no suptitle.
    """
    fig, (ax_lo, ax_hi) = plt.subplots(1, 2, figsize=(14, 6), sharey=False)

    sigma_low_min = 1e-3
    sigma_low_max = 1e4

    bm = benchmark_models()
    model_styles = {
        'M1_dark_photon_massive':   {'color': 'tab:green', 'marker': 'o',
                                     'label': r'Dark photon ($m_V=1.1$ keV, $v^*\!\approx\!200$ km/s)'},
        'M2_scalar_phi_massive':    {'color': 'tab:red',   'marker': 's',
                                     'label': r'Scalar $\phi$ ($m_\phi=6.95$ keV, $v^*\!\approx\!500$ km/s)'},
        'M3_massless_control':      {'color': 'tab:gray',  'marker': '^',
                                     'label': r'Massless control ($r_{\rm diss}=$const)'},
    }

    all_curves = []
    for mk, p in bm.items():
        if mk not in model_styles:
            continue
        result = scan_alpha_for_model(mk, p, alpha_log_range=(-3, 3), n_alpha=80)
        sigma_low = result['sigma_eff'][:, IDX_DWARF]
        rdiss_lo = result['rdiss_eff'][:, IDX_B1938]
        rdiss_hi = result['rdiss_eff'][:, IDX_CLUSTER]

        mask = np.isfinite(sigma_low) & np.isfinite(rdiss_lo) & np.isfinite(rdiss_hi)
        sigma_low = sigma_low[mask]
        rdiss_lo = rdiss_lo[mask]
        rdiss_hi = rdiss_hi[mask]
        order = np.argsort(sigma_low)
        sigma_low = sigma_low[order]
        rdiss_lo = rdiss_lo[order]
        rdiss_hi = rdiss_hi[order]

        st = model_styles[mk]
        for ax, rdiss_arr in [
            (ax_lo, rdiss_lo),
            (ax_hi, rdiss_hi),
        ]:
            ax.plot(sigma_low, rdiss_arr, color=st['color'], lw=2.0,
                    alpha=0.85)
            idx_pts = np.linspace(0, len(sigma_low) - 1, 6).astype(int)
            ax.scatter(sigma_low[idx_pts], rdiss_arr[idx_pts],
                       color=st['color'], marker=st['marker'], s=25, zorder=5)

        all_curves.append({
            'model': mk, 'sigma_low': sigma_low,
            'rdiss_eff_b1938': rdiss_lo, 'rdiss_eff_cluster': rdiss_hi,
        })

    # ---- Shared decorations on both panels ----
    for ax, v_scale_label, rdiss_arrs in [
        (ax_lo, r'$v \sim 80$ km/s (dwarf scale)',
         [c['rdiss_eff_b1938'] for c in all_curves]),
        (ax_hi, r'$v \sim 1200$ km/s (cluster scale)',
         [c['rdiss_eff_cluster'] for c in all_curves]),
    ]:
        # Constraint bands on x-axis
        ax.axvspan(50, sigma_low_max, color='red', alpha=0.04)
        ax.axvline(50, color='teal', linestyle='--', lw=0.9, alpha=0.5)
        ax.axvspan(sigma_low_min, 1.0, color='gray', alpha=0.04)
        ax.axvline(1.0, color='teal', linestyle='--', lw=0.9, alpha=0.5)

        # Elastic limit (r_diss = 1)
        ax.axhline(1.0, color='black', lw=1.0, linestyle='--', alpha=0.5)

        # Fiducial r_diss range
        ax.axhspan(1.01, 1.3, color='blue', alpha=0.03)

        ax.set_xscale('log')
        ax.set_xlim(sigma_low_min, sigma_low_max)

        # Y range: data-driven
        rmax = 1.05
        for arr in rdiss_arrs:
            if len(arr):
                rmax = max(rmax, float(np.max(arr)))
        ax.set_ylim(0.998, max(1.15, rmax * 1.08))

        ax.set_xlabel(r'$\sigma_T/m_\chi$ at dwarf scale ($v\sim 50$ km/s) [cm$^2$/g]',
                      fontsize=11)
        ax.set_title(v_scale_label, fontsize=11)
        ax.grid(True, which='both', alpha=0.2, zorder=0)

    ax_lo.set_ylabel(r'Effective $r_{\rm diss}$', fontsize=12)

    # Legend on the right panel only, placed in upper-right corner.
    # Data (r_diss ~ 1.00-1.05) sits near the bottom of the plot, so the
    # upper region is empty and the legend will not cover any curves.
    handles = []
    for mk, st in model_styles.items():
        handles.append(plt.Line2D([0], [0], color=st['color'], lw=2.0,
                                  marker=st['marker'], markersize=5,
                                  label=st['label']))
    ax_hi.legend(handles=handles, loc='upper right', fontsize=8,
                 framealpha=0.9, ncol=1)

    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {out_path}")
    return all_curves


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    print("[P4] Generating multi-model exclusion plot...")
    curves = plot_exclusion(os.path.join(OUT_FIG, 'P4_exclusion.png'))

    print("\n[P4] Generating σ/m vs velocity figure...")
    plot_sigma_vs_v(os.path.join(OUT_FIG, 'P4_sigma_vs_v.png'))

    print("\n[P4] Generating r_diss dimension figure...")
    rdiss_curves = plot_rdiss_dimension(
        os.path.join(OUT_FIG, 'P4_rdiss_dimension.png'))

    # Save curve data
    import pandas as pd
    rows = []
    for c in curves:
        for sl, sh_cl, sh_bul in zip(c['sigma_low'], c['sigma_high_cluster'],
                                      c['sigma_high_bullet']):
            rows.append({
                'model': c['model'],
                'sigma_low_dwarf': sl,
                'sigma_high_cluster': sh_cl,
                'sigma_high_bullet': sh_bul,
            })
    df = pd.DataFrame(rows)
    csv_path = os.path.join(OUT_DATA, 'P4_model_curves.csv')
    df.to_csv(csv_path, index=False)
    print(f"\nSaved curve data → {csv_path}")

    # Save r_diss dimension data
    rows_rd = []
    for c in rdiss_curves:
        for sl, rd_lo, rd_hi in zip(c['sigma_low'], c['rdiss_eff_b1938'],
                                    c['rdiss_eff_cluster']):
            rows_rd.append({
                'model': c['model'],
                'sigma_low_dwarf': sl,
                'rdiss_eff_b1938': rd_lo,
                'rdiss_eff_cluster': rd_hi,
            })
    df_rd = pd.DataFrame(rows_rd)
    csv_path_rd = os.path.join(OUT_DATA, 'P4_rdiss_curves.csv')
    df_rd.to_csv(csv_path_rd, index=False)
    print(f"Saved r_diss data → {csv_path_rd}")

    # ----- Quick assessment -----
    print("\n=== P4 Assessment ===")
    for c in curves:
        sl = c['sigma_low']
        sh = c['sigma_high_cluster']
        # Find points where σ_low in [1, 50] AND σ_high_cluster < 1
        allowed = (sl >= 1) & (sl <= 50) & (sh < 1)
        n_allowed = int(np.sum(allowed))
        print(f"  {c['model']}: {n_allowed} allowed points "
              f"(σ_low∈[1,50], σ_high_cl<1)")
        if n_allowed > 0:
            idx = np.where(allowed)[0]
            print(f"    σ_low range: {sl[idx].min():.3e} - {sl[idx].max():.3e}")
            print(f"    σ_high_cl range: {sh[idx].min():.3e} - {sh[idx].max():.3e}")


if __name__ == '__main__':
    main()
