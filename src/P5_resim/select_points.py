"""
Select 5 representative match points per model from P3 rescaling results.

For each model (M1_dark_photon_massive, M2_scalar_phi_massive, M3_massless_control),
we pick 5 points spanning the σ/m_phys range:
  - P1: best-fit (highest σ/m, closest to B1938 mass ratio)
  - P2: lower boundary (lowest σ/m)
  - P3: upper boundary (highest σ/m)
  - P4: median σ/m
  - P5: 75th percentile σ/m

For each point we save the FULL physical parameters needed for resimulation:
  - r_s_phys_kpc, rho_0_phys_Msun_pc3, sigma_phys_cm2_g, t_evo_phys_gyr
  - plus the original (lambda, mu, snapshot_time_gyr, r2D_rs) for reference

Output: data/P5_resim_points.csv with 15 rows (5 per model × 3 models)
"""
from __future__ import annotations

import os, sys
import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJ_ROOT = os.path.normpath(os.path.join(_HERE, '..', '..'))
P3_CSV = os.path.join(_PROJ_ROOT, 'data', 'P3_rescaled_params.csv')
OUT_CSV = os.path.join(_PROJ_ROOT, 'data', 'P5_resim_points.csv')

# Models to process — M3 corresponds to M3_massless_control which
# in P3 was not run directly; we use the elastic + const_rdiss_1p05
# runs as proxies for "massless / r_diss=const" models, and
# M1_dark_photon_massive and M2_scalar_phi_massive for the velocity-dependent
# massive emission. We additionally use the M1_highconc run for the
# concentration-scan task (handled separately).
TARGET_MODELS = [
    'M1_dark_photon_massive',
    'M2_scalar_phi_massive',
    # For "M3" we pick from elastic runs (r_diss=1) since the massless control
    # is just an elastic model with constant r_diss; both elastic and
    # const_rdiss_1p05 qualify, we use const_rdiss_1p05 as it has finite r_diss.
    'const_rdiss_1p05',
]
MODEL_LABELS = {
    'M1_dark_photon_massive': 'M1',
    'M2_scalar_phi_massive':  'M2',
    'const_rdiss_1p05':       'M3',  # massless control proxy
}


def select_5_points(df_sub: pd.DataFrame) -> pd.DataFrame:
    """Pick 5 representative points spanning the σ/m range.

    Filters out snapshot_idx=0 (initial NFW, t_evo=0) matches, which
    represent trivial coincidences where the initial mass ratio already
    falls in the observed band — these give meaningless 1-step resims.
    """
    # Filter to viable (t_ok=True) with actual evolution (t_evo > 0)
    viable = df_sub[df_sub['t_ok'] & (df_sub['t_evo_gyr'] > 0)]
    viable = viable.sort_values('sigma_m_cm2_g').reset_index(drop=True)
    n_filtered = len(df_sub[df_sub['t_ok']]) - len(viable)
    if n_filtered > 0:
        print(f"  Filtered out {n_filtered} points with t_evo=0 (snapshot_idx=0)")
    if len(viable) < 5:
        print(f"  WARNING: only {len(viable)} viable points with t_evo>0; using all")
        return viable.copy()

    n = len(viable)
    idx_picks = {
        'P1_lower':  0,           # lowest σ/m (lower boundary)
        'P2_p25':    n // 4,      # 25th percentile
        'P3_median': n // 2,     # median
        'P4_p75':    3 * n // 4, # 75th percentile
        'P5_upper':  n - 1,      # highest σ/m (upper boundary)
    }
    picks = []
    for label, idx in idx_picks.items():
        row = viable.iloc[idx].copy()
        row['pick_label'] = label
        picks.append(row)
    return pd.DataFrame(picks)


def main():
    if not os.path.exists(P3_CSV):
        print(f"ERROR: P3 CSV not found at {P3_CSV}")
        sys.exit(1)

    df = pd.read_csv(P3_CSV)
    print(f"Loaded {len(df)} P3 matches from {P3_CSV}")
    print(f"Available models: {df.model.unique()}")

    all_picks = []
    for mk in TARGET_MODELS:
        sub = df[df['model'] == mk]
        if len(sub) == 0:
            print(f"\n[{mk}] no matches; skipping")
            continue
        print(f"\n[{mk}] {len(sub)} matches, "
              f"σ/m range [{sub.sigma_m_cm2_g.min():.4e}, "
              f"{sub.sigma_m_cm2_g.max():.4e}]")
        picks = select_5_points(sub)
        picks['model_label'] = MODEL_LABELS[mk]
        all_picks.append(picks)
        print(f"  Selected 5 points:")
        for _, p in picks.iterrows():
            print(f"    {p.pick_label}: σ/m={p.sigma_m_cm2_g:.4e} cm²/g, "
                  f"t_evo={p.t_evo_gyr:.4e} Gyr, r_s={p.r_s_kpc:.4e} kpc, "
                  f"ρ₀={p.rho0_msun_pc3:.4e}, r2D/rs={p.r2D_rs:.4f}")

    if not all_picks:
        print("ERROR: No match points selected")
        sys.exit(1)

    out = pd.concat(all_picks, ignore_index=True)
    out.to_csv(OUT_CSV, index=False)
    print(f"\nWrote {len(out)} resim points → {OUT_CSV}")


if __name__ == '__main__':
    main()
