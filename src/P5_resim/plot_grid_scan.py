"""
Plot the self-consistent B1938+666 viable region from the 2D grid scan.

Generates:
  1. (sigma/m, r_s) viable region per model, with t_cross contours
  2. Comparison of self-consistent vs elastic-approximation viable region
  3. Corrected exclusion curve in (sigma/m, t_evo) plane
"""
from __future__ import annotations

import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.interpolate import griddata

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJ_ROOT = os.path.normpath(os.path.join(_HERE, '..', '..'))
GRID_CSV = os.path.join(_PROJ_ROOT, 'data', 'P5_grid_scan', 'grid_summary.csv')
P3_CSV = os.path.join(_PROJ_ROOT, 'data', 'P3_rescaled_params.csv')
OUT_DIR = os.path.join(_PROJ_ROOT, 'figures')

T_OBS_GYR = 6.37
MASS_RATIO_OBS = 0.36418166238217653
MASS_RATIO_ERR = 0.022

MODEL_CFG = {
    'M1': {'color': 'tab:green', 'marker': 'o',
           'label': r'M1: dark photon ($v^*\!=\!200$ km/s)'},
    'M2': {'color': 'tab:red', 'marker': 's',
           'label': r'M2: scalar $\phi$ ($v^*\!=\!500$ km/s)'},
    'M3': {'color': 'tab:blue', 'marker': '^',
           'label': r'M3: const $r_{\rm diss}\!=\!1.05$'},
}


def load_grid():
    if not os.path.exists(GRID_CSV):
        print(f"ERROR: {GRID_CSV} not found")
        sys.exit(1)
    df = pd.read_csv(GRID_CSV)
    # Filter out errors
    df = df[df.t_cross_gyr.notna()].copy()
    return df


def load_p3():
    if not os.path.exists(P3_CSV):
        return None
    return pd.read_csv(P3_CSV)


def plot_viable_region(df, out_path):
    """Viable region in (sigma/m, r_s) plane with t_cross contours."""
    fig, ax = plt.subplots(figsize=(9, 7))

    for ml, cfg in MODEL_CFG.items():
        sub = df[df.model_label == ml]
        if len(sub) == 0:
            continue

        # Plot all points colored by t_cross
        sig = sub['sigma_m_cm2_g'].values
        rs = sub['r_s_kpc'].values
        tc = sub['t_cross_gyr'].values
        # Replace inf with large value for plotting
        tc_plot = np.where(np.isinf(tc), 100.0, tc)
        tc_plot = np.where(tc_plot == 0, 0.001, tc_plot)

        sc = ax.scatter(sig, rs, c=tc_plot, cmap='viridis_r', s=30,
                        norm=LogNorm(vmin=0.01, vmax=100),
                        marker=cfg['marker'], edgecolors='k',
                        linewidths=0.3, alpha=0.8, label=cfg['label'])

        # Highlight viable points (t_cross <= T_OBS)
        viable = sub[sub.viable == True]
        if len(viable) > 0:
            ax.scatter(viable['sigma_m_cm2_g'], viable['r_s_kpc'],
                       s=100, facecolors='none', edgecolors=cfg['color'],
                       linewidths=2.5, zorder=5)

    # t_obs line on colorbar
    cb = plt.colorbar(sc, ax=ax, shrink=0.8, pad=0.02)
    cb.set_label(r'$t_{\rm cross}$ [Gyr]', fontsize=11)
    cb.ax.axhline(T_OBS_GYR, color='red', lw=2, ls='--')
    cb.ax.text(0.5, T_OBS_GYR, r'$t_{\rm obs}$', color='red',
               fontsize=8, transform=cb.ax.get_yaxis_transform())

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$\sigma_T/m_\chi$ [cm$^2$/g]', fontsize=12)
    ax.set_ylabel(r'$r_s$ [kpc]', fontsize=12)
    ax.set_title('Self-consistent B1938+666 viable region\n'
                 r'(circles = viable, $t_{\rm cross} \leq t_{\rm obs} = 6.37$ Gyr)',
                 fontsize=11)
    ax.legend(loc='best', fontsize=9, framealpha=0.9)
    ax.grid(True, alpha=0.3, which='both')

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {out_path}")


def plot_comparison(df, p3_df, out_path):
    """Compare self-consistent vs elastic-approximation viable region."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left: self-consistent (grid scan)
    ax = axes[0]
    for ml, cfg in MODEL_CFG.items():
        sub = df[df.model_label == ml]
        viable = sub[sub.viable == True]
        if len(viable) > 0:
            ax.scatter(viable['sigma_m_cm2_g'], viable['r_s_kpc'],
                       s=40, c=cfg['color'], marker=cfg['marker'],
                       alpha=0.6, label=cfg['label'])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$\sigma_T/m_\chi$ [cm$^2$/g]', fontsize=11)
    ax.set_ylabel(r'$r_s$ [kpc]', fontsize=11)
    ax.set_title('Self-consistent (P5 grid scan)', fontsize=11)
    ax.legend(loc='best', fontsize=9)
    ax.grid(True, alpha=0.3, which='both')

    # Right: elastic approximation (P3 matches)
    ax = axes[1]
    if p3_df is not None:
        p3_viable = p3_df[p3_df.t_ok == True]
        for ml, cfg in MODEL_CFG.items():
            # Map P3 model names to M1/M2/M3
            p3_key = {'M1': 'M1_dark_photon_massive',
                      'M2': 'M2_scalar_phi_massive',
                      'M3': 'const_rdiss_1p05'}[ml]
            sub = p3_viable[p3_viable.model == p3_key]
            if len(sub) > 0:
                ax.scatter(sub['sigma_m_cm2_g'], sub['r_s_kpc'],
                           s=10, c=cfg['color'], marker=cfg['marker'],
                           alpha=0.4, label=cfg['label'])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$\sigma_T/m_\chi$ [cm$^2$/g]', fontsize=11)
    ax.set_ylabel(r'$r_s$ [kpc]', fontsize=11)
    ax.set_title('Elastic approximation (P3 rescaling)', fontsize=11)
    ax.legend(loc='best', fontsize=9)
    ax.grid(True, alpha=0.3, which='both')

    fig.suptitle('B1938+666 viable region: self-consistent vs elastic approximation',
                 fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {out_path}")


def plot_tcross_contours(df, out_path):
    """t_cross contour plot per model."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=True)

    for ax, (ml, cfg) in zip(axes, MODEL_CFG.items()):
        sub = df[df.model_label == ml]
        if len(sub) < 4:
            ax.set_title(f'{ml}: insufficient data')
            continue

        sig = sub['sigma_m_cm2_g'].values
        rs = sub['r_s_kpc'].values
        tc = sub['t_cross_gyr'].values
        tc_plot = np.where(np.isinf(tc), 100.0, tc)
        tc_plot = np.where(tc_plot == 0, 0.001, tc_plot)

        # Grid for contour
        sig_grid = np.logspace(np.log10(sig.min()), np.log10(sig.max()), 50)
        rs_grid = np.logspace(np.log10(rs.min()), np.log10(rs.max()), 50)
        SIG, RS = np.meshgrid(sig_grid, rs_grid)

        points = np.column_stack([np.log10(sig), np.log10(rs)])
        values = np.log10(tc_plot)
        grid_points = np.column_stack([np.log10(SIG.ravel()), np.log10(RS.ravel())])
        TC = griddata(points, values, grid_points, method='linear').reshape(SIG.shape)

        # Contour plot
        levels = np.logspace(-2, 2, 9)
        cs = ax.contourf(SIG, RS, 10**TC, levels=levels, cmap='viridis_r',
                         norm=LogNorm(vmin=0.01, vmax=100))
        ax.contour(SIG, RS, 10**TC, levels=[T_OBS_GYR], colors='red',
                   linewidths=2, linestyles='--')
        ax.scatter(sig, rs, c='white', s=20, edgecolors='k', linewidths=0.5,
                   zorder=5)

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(r'$\sigma_T/m_\chi$ [cm$^2$/g]', fontsize=11)
        if ax == axes[0]:
            ax.set_ylabel(r'$r_s$ [kpc]', fontsize=11)
        ax.set_title(f'{ml}', fontsize=11)
        ax.grid(True, alpha=0.3, which='both')

    cb = fig.colorbar(cs, ax=axes, shrink=0.8, pad=0.02)
    cb.set_label(r'$t_{\rm cross}$ [Gyr]', fontsize=11)

    fig.suptitle(r'$t_{\rm cross}$ contours (red dashed = $t_{\rm obs} = 6.37$ Gyr)',
                 fontsize=12, y=1.02)
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {out_path}")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    df = load_grid()
    print(f"Loaded {len(df)} grid points")
    if len(df) == 0:
        return

    print("\nViable points per model:")
    for ml in ['M1', 'M2', 'M3']:
        sub = df[df.model_label == ml]
        v = sub.viable.sum()
        print(f"  {ml}: {v}/{len(sub)} viable")
        if v > 0:
            viable = sub[sub.viable == True]
            print(f"    sigma/m range: [{viable.sigma_m_cm2_g.min():.4e}, "
                  f"{viable.sigma_m_cm2_g.max():.4e}]")
            print(f"    r_s range: [{viable.r_s_kpc.min():.4f}, "
                  f"{viable.r_s_kpc.max():.4f}]")

    p3_df = load_p3()

    out1 = os.path.join(OUT_DIR, 'P5_viable_region.png')
    plot_viable_region(df, out1)

    out2 = os.path.join(OUT_DIR, 'P5_viable_comparison.png')
    plot_comparison(df, p3_df, out2)

    out3 = os.path.join(OUT_DIR, 'P5_tcross_contours.png')
    plot_tcross_contours(df, out3)


if __name__ == '__main__':
    main()
