#!/usr/bin/env python3
"""Analyze the working reference snapshots to compute projected mass ratios.

The working reference used:
  - Cored NFW IC (r_s=3.6 kpc, rho_0=7.09e-3, core radius ~0.3 kpc)
  - sigma/m = 3.0 cm^2/g (similar to our P1: 3.33)
  - r_diss = 1.05 (slightly dissipative)
  - 1M DM + 500k gas particles
  - Evolved to t=1.0 code units = 0.978 Gyr

This gives us a reference for:
  1. The initial projected mass ratio (t=0)
  2. The final ratio after SIDM evolution (t=1.0)
  3. Comparison with our 3 new runs (sigma/m = 3.33, 0.167, 7.33)

Note: The reference IC is cored NFW, not pure NFW.
The fluid model uses pure NFW, so the absolute ratios will differ.
But the RELATIVE CHANGE (ratio_final / ratio_init) can be compared.
"""
import numpy as np
import struct
import os
import sys

# Add parent dir for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from read_binary_snap import read_snapshot

# Rescaling parameters
LAMBDA = 0.085 / 3.6
MU = (10.0 / 7.09e-3) * LAMBDA**3
T_SCALE = np.sqrt(LAMBDA**3 / MU)

# Rescaled radii in simulation space
R_INNER_SIM = 0.020 / LAMBDA   # 0.847 kpc
R_OUTER_SIM = 0.090 / LAMBDA   # 3.812 kpc

# Fluid model predictions (for pure NFW)
FLUID_PREDICTIONS = {
    'P1_elastic_control': {'ratio_init': 0.1465, 'ratio_final': 0.1205, 't_actual': 0.233},
    'P2_m3_low_sigma':    {'ratio_init': 0.1465, 'ratio_final': 0.1480, 't_actual': 0.040},
    'P3_m3_high_sigma':   {'ratio_init': 0.1465, 'ratio_final': 0.1321, 't_actual': 0.100},
}


def projected_mass_ratio(coords, masses, r_inner, r_outer):
    """Compute M_2D(r_inner)/M_2D(r_outer) averaged over 3 axes."""
    ratios = []
    for axis in ['z', 'y', 'x']:
        if axis == 'z':
            rp = np.sqrt(coords[:, 0]**2 + coords[:, 1]**2)
        elif axis == 'y':
            rp = np.sqrt(coords[:, 0]**2 + coords[:, 2]**2)
        else:
            rp = np.sqrt(coords[:, 1]**2 + coords[:, 2]**2)
        mi = masses[rp <= r_inner].sum()
        mo = masses[rp <= r_outer].sum()
        if mo > 0:
            ratios.append(mi / mo)
    return np.mean(ratios), np.std(ratios)


def analyze_snapshot(snap_file, label=""):
    """Analyze a single snapshot."""
    print(f"\n--- {label}: {os.path.basename(snap_file)} ---")
    try:
        snap = read_snapshot(snap_file, types_to_read=(1, 2))
    except Exception as e:
        print(f"  ERROR: {e}")
        return None

    header = snap['header']
    t_code = header.get('time', 0.0)
    t_sim_gyr = t_code * 0.978
    t_phys = t_sim_gyr * T_SCALE

    results = {}
    for ptype in snap['coords']:
        coords = snap['coords'][ptype]
        if coords is None or len(coords) == 0:
            continue
        masses = snap['masses'][ptype]
        if masses is None:
            # Use header mass
            mass_table = header.get('mass', [0]*6)
            if ptype < len(mass_table) and mass_table[ptype] > 0:
                masses = np.full(len(coords), mass_table[ptype])
            else:
                continue

        n = len(coords)
        print(f"  PartType{ptype}: N={n}, t_code={t_code:.4f}, t_sim={t_sim_gyr:.4f} Gyr, t_phys={t_phys:.4f} Gyr")

        ratio, ratio_std = projected_mass_ratio(coords, masses, R_INNER_SIM, R_OUTER_SIM)
        print(f"    M_2D(<{R_INNER_SIM:.3f})/M_2D(<{R_OUTER_SIM:.3f}) = {ratio:.6f} +/- {ratio_std:.6f}")

        results[f'ratio_ptype{ptype}'] = ratio
        results[f'ratio_std_ptype{ptype}'] = ratio_std
        results[f'n_ptype{ptype}'] = n

    results['t_code'] = t_code
    results['t_sim_gyr'] = t_sim_gyr
    results['t_phys_gyr'] = t_phys
    return results


def main():
    ref_dir = "D:/graverthermal-sidm/data/P5_nbody_verify/ref_snapshots"
    snaps = sorted([f for f in os.listdir(ref_dir) if f.startswith('snapshot_')])

    print("=" * 70)
    print("Working reference snapshot analysis")
    print("=" * 70)
    print(f"R_inner_sim = {R_INNER_SIM:.4f} kpc (from 20 pc)")
    print(f"R_outer_sim = {R_OUTER_SIM:.4f} kpc (from 90 pc)")
    print(f"Found {len(snaps)} snapshots: {snaps}")

    all_results = []
    for snap in snaps:
        path = os.path.join(ref_dir, snap)
        r = analyze_snapshot(path, label="Reference")
        if r:
            r['file'] = snap
            all_results.append(r)

    # Summary
    if len(all_results) >= 2:
        print("\n" + "=" * 70)
        print("Summary: Reference (sigma/m=3.0, r_diss=1.05)")
        print("=" * 70)

        r0 = all_results[0]  # t=0
        r1 = all_results[-1]  # t=1.0

        for pt in [1, 2]:
            key = f'ratio_ptype{pt}'
            if key in r0 and key in r1:
                ratio_init = r0[key]
                ratio_final = r1[key]
                if ratio_init > 0:
                    change = (ratio_final - ratio_init) / ratio_init * 100
                    print(f"\nPartType{pt}:")
                    print(f"  ratio_init  (t=0)   = {ratio_init:.6f}")
                    print(f"  ratio_final (t=1.0) = {ratio_final:.6f}")
                    print(f"  Change = {change:+.2f}%")

        # Compare with fluid model
        print(f"\n--- Comparison with fluid model (P1: sigma/m=3.33, elastic) ---")
        fp = FLUID_PREDICTIONS['P1_elastic_control']
        if 'ratio_ptype1' in r0 and 'ratio_ptype1' in r1:
            # The N-body ratio_init will differ from fluid's 0.1465 because the IC is cored
            # But we can compare the RELATIVE CHANGE
            nb_init = r1.get('ratio_ptype1', 0)
            nb_final = r1.get('ratio_ptype1', 0)

            # Actually, r0 is t=0 and r1 is t=1.0
            nb_init = r0['ratio_ptype1']
            nb_final = r1['ratio_ptype1']
            nb_change = (nb_final - nb_init) / nb_init * 100

            fluid_init = fp['ratio_init']
            fluid_final = fp['ratio_final']
            fluid_change = (fluid_final - fluid_init) / fluid_init * 100

            print(f"  Fluid (pure NFW): init={fluid_init:.4f}, final={fluid_final:.4f}, change={fluid_change:+.2f}%")
            print(f"  N-body (cored):  init={nb_init:.6f}, final={nb_final:.6f}, change={nb_change:+.2f}%")
            print(f"  Note: Absolute ratios differ due to cored vs pure NFW profile")
            print(f"  The relative change comparison is more meaningful")


if __name__ == '__main__':
    main()
