"""
Plot the symmetry-breaking offset from the P5 resim results.

Compares the resim mass ratio (M(20pc)/M(90pc) from a fresh physical-parameter
simulation) to the observed value (0.364 +/- 0.022). The offset
Delta = (M_resim - M_obs)/M_obs is the ACTUAL symmetry-breaking correction,
which can be compared to the first-order steady-state estimate from P3.

Outputs:
  figures/P5_resim_offset.png        : Delta vs sigma/m, per model
  figures/P5_resim_correction.png    : sigma/m_corrected vs sigma/m_naive
"""
from __future__ import annotations

import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJ_ROOT = os.path.normpath(os.path.join(_HERE, '..', '..'))
SUMMARY_CSV = os.path.join(_PROJ_ROOT, 'data', 'P5_resim', 'summary.csv')
P3_CSV = os.path.join(_PROJ_ROOT, 'data', 'P3_rescaled_params.csv')
OUT_DIR = os.path.join(_PROJ_ROOT, 'figures')

MASS_RATIO_OBS = 0.36418166238217653
MASS_RATIO_ERR = 0.022

# Model display config
MODEL_CFG = {
    'M1': {'color': 'tab:green', 'marker': 'o',
           'label': r'M1: dark photon ($v^*\!=\!200$ km/s)'},
    'M2': {'color': 'tab:red', 'marker': 's',
           'label': r'M2: scalar $\phi$ ($v^*\!=\!500$ km/s)'},
    'M3': {'color': 'tab:blue', 'marker': '^',
           'label': r'M3: const $r_{\rm diss}\!=\!1.05$ (control)'},
}


def load_resim():
    if not os.path.exists(SUMMARY_CSV):
        print(f"ERROR: {SUMMARY_CSV} not found")
        sys.exit(1)
    df = pd.read_csv(SUMMARY_CSV)
    # Filter out error rows
    df = df[df.get('mass_ratio_resim', pd.Series([np.nan])).notna()].copy()
    return df


def load_p3_first_order():
    """Load the P3 first-order symmetry-breaking estimates for comparison."""
    if not os.path.exists(P3_CSV):
        return None
    df = pd.read_csv(P3_CSV)
    vd = df[df['model'].isin(['M1_dark_photon_massive',
                              'M2_scalar_phi_massive'])].copy()
    vd = vd[vd['valid_regime'] == True]
    if len(vd) == 0:
        return None
    return vd


def plot_delta_vs_sigma(df, out_path):
    """Delta (mass ratio offset) vs sigma/m, per model.

    This is the central result: how much does the naive elastic rescaling
    miss the true dissipative mass ratio?
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    # Observed band
    ax.axhspan(-MASS_RATIO_ERR / MASS_RATIO_OBS * 100,
               +MASS_RATIO_ERR / MASS_RATIO_OBS * 100,
               color='gray', alpha=0.2, label=r'Observed ($\pm 1\sigma$)')

    for ml, cfg in MODEL_CFG.items():
        sub = df[df['model_label'] == ml]
        if len(sub) == 0:
            continue
        sub = sub.sort_values('sigma_m_cm2_g')
        ax.errorbar(sub['sigma_m_cm2_g'], sub['delta_ratio_pct'],
                   yerr=0,  # no formal error bar on resim
                   fmt=cfg['marker'] + '-', color=cfg['color'],
                   markersize=8, capsize=3, lw=1.5,
                   label=cfg['label'])

    ax.axhline(0, color='k', lw=0.8, ls='-')
    ax.set_xscale('log')
    ax.set_xlabel(r'$\sigma_T/m_\chi$ [cm$^2$/g] (naive, from elastic rescale)',
                  fontsize=12)
    ax.set_ylabel(r'$\Delta M_{20}/M_{90}$ [%]  '
                  r'$= (M_{\rm resim} - M_{\rm obs})/M_{\rm obs}$',
                  fontsize=11)
    ax.set_title('P5: Symmetry-breaking offset from per-point resimulation\n'
                 r'(0% = elastic rescaling valid; '
                 r'$\Delta<0$ = over-cooled, $\sigma/m$ overestimated)',
                 fontsize=11)
    ax.legend(loc='best', fontsize=9, framealpha=0.9)
    ax.grid(True, alpha=0.3, which='both')

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {out_path}")


def plot_correction_factor(df, out_path):
    """sigma/m_corrected vs sigma/m_naive.

    The correction factor is derived from Delta:
      - If Delta < 0 (over-cooled): sigma/m_corrected < sigma/m_naive
      - Approximate: sigma/m_corrected ~ sigma/m_naive * sqrt(M_obs/M_resim)
        (since t_cool ~ 1/sigma in LMFP limit, M_ratio ~ sigma^alpha)
    """
    fig, ax = plt.subplots(figsize=(7, 7))

    # Elastic diagonal
    sig_arr = np.logspace(-3, 0, 50)
    ax.plot(sig_arr, sig_arr, 'k--', lw=1.0,
           label='Elastic symmetry (no correction)')

    for ml, cfg in MODEL_CFG.items():
        sub = df[df['model_label'] == ml]
        if len(sub) == 0:
            continue
        sub = sub.sort_values('sigma_m_cm2_g')

        # Approximate correction: sigma_corrected = sigma_naive * (M_obs/M_resim)^0.5
        # (heuristic; the true relation requires a full sigma scan)
        ratio_obs_over_resim = MASS_RATIO_OBS / sub['mass_ratio_resim']
        sigma_corrected = sub['sigma_m_cm2_g'] * np.sqrt(ratio_obs_over_resim)

        ax.plot(sub['sigma_m_cm2_g'], sigma_corrected,
                cfg['marker'] + '-', color=cfg['color'],
                markersize=8, lw=1.5, label=cfg['label'])

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$\sigma_T/m_\chi$ naive [cm$^2$/g] '
                  r'(from elastic rescaling)', fontsize=11)
    ax.set_ylabel(r'$\sigma_T/m_\chi$ corrected [cm$^2$/g] '
                  r'(from resim $\Delta$)', fontsize=11)
    ax.set_title('P5: Symmetry-breaking correction to $\\sigma_T/m_\\chi$\n'
                 '(heuristic: $\\sigma_{\\rm corr} = \\sigma_{\\rm naive} '
                 '\\sqrt{M_{\\rm obs}/M_{\\rm resim}}$)',
                 fontsize=11)
    ax.legend(loc='best', fontsize=9, framealpha=0.9)
    ax.grid(True, alpha=0.3, which='both')

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {out_path}")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    df = load_resim()
    print(f"Loaded {len(df)} resim results")
    if len(df) == 0:
        print("No results to plot")
        return

    print("\nSummary:")
    cols = ['model_label', 'pick_label', 'sigma_m_cm2_g', 't_evo_gyr',
            'mass_ratio_resim', 'delta_ratio_pct']
    cols = [c for c in cols if c in df.columns]
    print(df[cols].to_string(index=False))

    out1 = os.path.join(OUT_DIR, 'P5_resim_offset.png')
    plot_delta_vs_sigma(df, out1)

    out2 = os.path.join(OUT_DIR, 'P5_resim_correction.png')
    plot_correction_factor(df, out2)


if __name__ == '__main__':
    main()
