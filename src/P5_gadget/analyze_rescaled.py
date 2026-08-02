#!/usr/bin/env python3
"""Analyze rescaled N-body snapshots and compare with fluid model.

The projected mass ratio M_2D(R_inner)/M_2D(R_outer) is INVARIANT under the
Schmidt et al. 2026 Appendix G rescaling symmetry:
  - Both radii scale by 1/lambda: R_sim = R_phys / lambda
  - Both masses scale by mu: M_sim = mu * M_phys... wait, actually
    M_2D scales as mu (mass) * lambda^2 (length^2 from the surface density
    integral? No - the 2D projected mass integrates the 3D density along
    line of sight, so M_2D ~ rho * length^3 ~ (mu/lambda^3) * lambda^3 = mu)
  - So the ratio M_2D(R_inner)/M_2D(R_outer) is invariant.

We can therefore compute the ratio directly in simulation space and compare
with the fluid model's ratio (also invariant under rescaling).

Usage:
    python analyze_rescaled.py <snapshot_dir>
    python analyze_rescaled.py data/P5_nbody_verify/sim_snapshots/

Reads all snapshot_NNN files in the directory, computes the projected mass
ratio at rescaled radii, and prints a comparison table.
"""
import numpy as np
import struct
import os
import sys
import glob

# Add parent dir for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from read_binary_snap import read_snapshot

# ============================================================
# Rescaling parameters
# ============================================================
LAMBDA = 0.085 / 3.6
MU = (10.0 / 7.09e-3) * LAMBDA**3
T_SCALE = np.sqrt(LAMBDA**3 / MU)  # t_phys = T_SCALE * t_sim

# Rescaled radii in simulation space
R_INNER_SIM = 0.020 / LAMBDA   # 0.847 kpc (from 20 pc)
R_OUTER_SIM = 0.090 / LAMBDA   # 3.812 kpc (from 90 pc)

# Fluid model predictions at N-body time checkpoints (t_code=0, 0.5, 1.0)
# These were computed by compute_fluid_at_nbody_time.py and match the
# actual physical time our N-body runs reach.
# ratio_init = 0.146506 for all points (pure NFW, r_s=0.085 kpc)
# Our IC is CORED NFW (ratio_init ~ 0.084), so we compare RELATIVE CHANGE.
FLUID_PREDICTIONS = {
    'P1_elastic_control': {
        'ratio_init': 0.146506,
        'ratio_t05': 0.140624,   # at t_code=0.5, t_phys=0.013 Gyr
        'ratio_final': 0.142738,  # at t_code=1.0, t_phys=0.026 Gyr
        't_actual': 0.026,
        'rel_change': 0.142738 / 0.146506,  # 0.9742 (-2.58%)
    },
    'P2_m3_low_sigma':    {
        'ratio_init': 0.146506,
        'ratio_t05': 0.147932,
        'ratio_final': 0.148193,
        't_actual': 0.026,
        'rel_change': 0.148193 / 0.146506,  # 1.0115 (+1.15%)
    },
    'P3_m3_high_sigma':   {
        'ratio_init': 0.146506,
        'ratio_t05': 0.145119,
        'ratio_final': 0.141516,
        't_actual': 0.026,
        'rel_change': 0.141516 / 0.146506,  # 0.9659 (-3.41%)
    },
}

# Working reference (sigma/m=3.0, no SIDM) for comparison
# DM ratio: 0.0841 -> 0.0861, +2.36% change over t=0->1.0 code units
REF_REL_CHANGE = 0.0861 / 0.0841  # 1.0238


def projected_mass_ratio(coords, masses, r_inner, r_outer, avg_axes=True):
    """Compute M_2D(r_inner)/M_2D(r_outer) from N-body particles."""
    if avg_axes:
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
    else:
        rp = np.sqrt(coords[:, 0]**2 + coords[:, 1]**2)
        mi = masses[rp <= r_inner].sum()
        mo = masses[rp <= r_outer].sum()
        return (mi / mo if mo > 0 else np.nan), 0.0


def analyze_snapshot(snap_file, point_name=None):
    """Analyze a single snapshot file."""
    print(f"\n--- Analyzing: {snap_file} ---")
    try:
        snap = read_snapshot(snap_file, types_to_read=(1,))
    except Exception as e:
        print(f"  ERROR reading snapshot: {e}")
        return None

    coords = snap['coords'][1]
    masses = snap['masses'][1]
    header = snap['header']

    t_code = header.get('time', 0.0)
    t_sim_gyr = t_code * 0.978  # code units to Gyr
    t_phys = t_sim_gyr * T_SCALE  # map back to physical time

    n = len(coords)
    print(f"  N = {n} particles, t_code = {t_code:.4f}, t_sim = {t_sim_gyr:.4f} Gyr, t_phys = {t_phys:.4f} Gyr")

    ratio, ratio_std = projected_mass_ratio(
        coords, masses, R_INNER_SIM, R_OUTER_SIM, avg_axes=True)

    # 3D spherical for comparison
    r3d = np.sqrt((coords**2).sum(axis=1))
    m3d_in = masses[r3d <= R_INNER_SIM].sum()
    m3d_out = masses[r3d <= R_OUTER_SIM].sum()
    ratio_3d = m3d_in / m3d_out if m3d_out > 0 else np.nan

    print(f"  M_2D(<{R_INNER_SIM:.3f} kpc) / M_2D(<{R_OUTER_SIM:.3f} kpc) = {ratio:.6f} +/- {ratio_std:.6f}")
    print(f"  M_3D(<{R_INNER_SIM:.3f} kpc) / M_3D(<{R_OUTER_SIM:.3f} kpc) = {ratio_3d:.6f}")

    result = {
        'file': snap_file,
        't_code': t_code,
        't_sim_gyr': t_sim_gyr,
        't_phys_gyr': t_phys,
        'ratio_2d': ratio,
        'ratio_2d_std': ratio_std,
        'ratio_3d': ratio_3d,
        'n_particles': n,
    }

    # Note: absolute ratio comparison is NOT valid because our IC (cored NFW)
    # differs from the fluid model (pure NFW). We compare RELATIVE CHANGE instead.
    if point_name and point_name in FLUID_PREDICTIONS:
        fp = FLUID_PREDICTIONS[point_name]
        print(f"\n  Fluid model: ratio {fp['ratio_init']:.4f} -> {fp['ratio_final']:.4f} "
              f"(rel change {fp['rel_change']:.4f}, i.e. {(fp['rel_change']-1)*100:+.2f}%)")
        print(f"  NOTE: Absolute ratios not comparable (cored vs pure NFW IC).")
        print(f"        Use relative-change comparison after snapshot_000 is analyzed.")

    return result


def main():
    if len(sys.argv) < 2:
        print("Usage: analyze_rescaled.py <snapshot_dir> [point_name]")
        print("  e.g. analyze_rescaled.py data/P5_nbody_verify/sim_snapshots/P1 P1_elastic_control")
        sys.exit(1)

    snap_dir = sys.argv[1]
    point_name = sys.argv[2] if len(sys.argv) > 2 else None

    if not os.path.isdir(snap_dir):
        print(f"ERROR: directory not found: {snap_dir}")
        sys.exit(1)

    # Find all snapshot files (binary format: snapshot_000, snapshot_001, ...)
    snap_files = sorted(glob.glob(os.path.join(snap_dir, "snapshot_*")))
    # Also check HDF5 format
    if not snap_files:
        snap_files = sorted(glob.glob(os.path.join(snap_dir, "snap_*.hdf5")))

    if not snap_files:
        print(f"No snapshot files found in {snap_dir}")
        sys.exit(1)

    print("=" * 70)
    print("Rescaling-symmetry N-body verification analysis")
    print("=" * 70)
    print(f"Snapshot directory: {snap_dir}")
    print(f"Point: {point_name or 'unknown'}")
    print(f"R_inner_sim = {R_INNER_SIM:.4f} kpc (from 20 pc)")
    print(f"R_outer_sim = {R_OUTER_SIM:.4f} kpc (from 90 pc)")
    print(f"Found {len(snap_files)} snapshot(s)")

    results = []
    for sf in snap_files:
        r = analyze_snapshot(sf, point_name)
        if r:
            results.append(r)

    # Summary table
    if results:
        print("\n" + "=" * 70)
        print("Summary")
        print("=" * 70)
        print(f"{'Snapshot':<30} {'t_code':>8} {'t_sim(Gyr)':>10} {'t_phys(Gyr)':>11} {'ratio_2d':>10} {'ratio_3d':>10}")
        print("-" * 85)
        for r in results:
            print(f"{os.path.basename(r['file']):<30} {r['t_code']:>8.4f} {r['t_sim_gyr']:>10.4f} "
                  f"{r['t_phys_gyr']:>11.4f} {r['ratio_2d']:>10.6f} {r['ratio_3d']:>10.6f}")

        if point_name and point_name in FLUID_PREDICTIONS:
            fp = FLUID_PREDICTIONS[point_name]
            
            # Find snapshot_000 (t=0) for initial ratio
            init_result = None
            for r in results:
                if r['t_code'] < 1e-8 or 'snapshot_000' in os.path.basename(r['file']):
                    init_result = r
                    break
            
            # Find the final snapshot (largest t_code)
            final_result = max(results, key=lambda r: r['t_code'])
            
            if init_result and final_result and final_result is not init_result:
                nb_rel_change = final_result['ratio_2d'] / init_result['ratio_2d']
                fluid_rel_change = fp['rel_change']
                # Also subtract the reference (no-SIDM) change to isolate SIDM effect
                nb_sidm_effect = nb_rel_change / REF_REL_CHANGE
                
                print(f"\n{'='*70}")
                print(f"RELATIVE CHANGE COMPARISON (corrected for cored vs pure NFW)")
                print(f"{'='*70}")
                print(f"  N-body:")
                print(f"    ratio_init  = {init_result['ratio_2d']:.6f} (snapshot_000)")
                print(f"    ratio_final = {final_result['ratio_2d']:.6f} "
                      f"(t_code={final_result['t_code']:.4f}, t_phys={final_result['t_phys_gyr']:.4f} Gyr)")
                print(f"    rel_change  = {nb_rel_change:.4f} ({(nb_rel_change-1)*100:+.2f}%)")
                print(f"  Fluid model (pure NFW):")
                print(f"    ratio_init  = {fp['ratio_init']:.4f}")
                print(f"    ratio_final = {fp['ratio_final']:.4f} (t={fp['t_actual']:.3f} Gyr)")
                print(f"    rel_change  = {fluid_rel_change:.4f} ({(fluid_rel_change-1)*100:+.2f}%)")
                print(f"  Working reference (no SIDM, sigma/m=3.0):")
                print(f"    rel_change  = {REF_REL_CHANGE:.4f} ({(REF_REL_CHANGE-1)*100:+.2f}%)")
                print(f"\n  N-body SIDM-only effect (after subtracting gravity):")
                print(f"    nb_rel_change / ref_rel_change = {nb_sidm_effect:.4f} "
                      f"({(nb_sidm_effect-1)*100:+.2f}%)")
                
                # Verdict based on relative-change comparison
                # For elastic control (P1): expect ratio to DECREASE (cusp formation)
                # For low sigma (P2): expect little change
                # For high sigma (P3): expect ratio to DECREASE
                delta_rel = abs(nb_rel_change - fluid_rel_change) / abs(fluid_rel_change - 1) if abs(fluid_rel_change - 1) > 0.01 else 0
                print(f"\n  |nb_rel - fluid_rel| / |fluid_rel - 1| = {delta_rel:.3f}")
                verdict = 'PASS' if delta_rel < 0.3 else 'MARGINAL' if delta_rel < 0.6 else 'FAIL'
                print(f"  Verdict: {verdict} (threshold: 0.3 pass, 0.6 marginal)")
            else:
                print(f"\n  Need at least 2 snapshots (init + final) for relative change comparison.")


if __name__ == '__main__':
    main()
