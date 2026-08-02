#!/usr/bin/env python3
"""Analyze the P1 relaxation rerun snapshots and compare with original + fluid.

Reads the three Phase-B snapshots (t=0, 0.5, 1.0) from output_prod/, computes
the projected enclosed mass ratio M_2D(<R_inner)/M_2D(<R_outer) at each time,
and writes a comparison table against:

  - The original P1 run (sim_snapshots/P1_elastic_control/)
  - The fluid-model predictions (fluid_predictions_nbody_time.csv)

The projected ratio is invariant under the Schmidt 2026 Appendix-G rescaling,
so simulation-space radii are used directly:

  R_inner_sim = 0.847 kpc  (= 20 pc / lambda)
  R_outer_sim = 3.812 kpc  (= 90 pc / lambda)

Usage:
    python analyze_relaxation.py
"""
import os
import sys
import struct
import numpy as np

# --- Paths (relative to this file's directory) ---
HERE = os.path.dirname(os.path.abspath(__file__))
DATA_ROOT = os.path.normpath(os.path.join(HERE, '..'))

# Simulation-space projection radii (invariant under rescaling)
R_INNER_SIM = 0.847   # kpc, from 20 pc
R_OUTER_SIM = 3.812   # kpc, from 90 pc

HEADER_SIZE = 256


def read_block(f):
    sz = struct.unpack('I', f.read(4))[0]
    data = f.read(sz)
    sz2 = struct.unpack('I', f.read(4))[0]
    if sz != sz2:
        raise ValueError(f"Block size mismatch: {sz} vs {sz2}")
    return data


def parse_header(hdr_buf):
    fmt = '<6I6d2d2i6I2i4d2i6Ii'
    sz = struct.calcsize(fmt)
    vals = struct.unpack(fmt, hdr_buf[:sz])
    h = {}
    h['npart'] = list(vals[0:6])
    h['mass'] = list(vals[6:12])
    h['time'] = vals[12]
    h['redshift'] = vals[13]
    return h


def read_snapshot(path):
    """Read PartType1 coords/vels/masses from a Gadget4 binary snapshot."""
    with open(path, 'rb') as f:
        hdr_buf = read_block(f)
        h = parse_header(hdr_buf)
        n1 = h['npart'][1]
        m1 = h['mass'][1]

        coords = read_block(f)
        vels = read_block(f)
        ids = read_block(f)

        coords = np.frombuffer(coords, dtype=np.float32).reshape(-1, 3)
        vels   = np.frombuffer(vels,   dtype=np.float32).reshape(-1, 3)
        ids    = np.frombuffer(ids,    dtype=np.uint32)

        # If per-type mass is 0, the Masses block follows.
        if m1 == 0.0 and n1 > 0:
            mass_blk = read_block(f)
            masses = np.frombuffer(mass_blk, dtype=np.float32)
        else:
            masses = np.full(n1, m1, dtype=np.float32)

    return h, coords, vels, ids, masses


def projected_ratio(coords, masses, r_inner, r_outer):
    """Average M_2D(<r_inner)/M_2D(<r_outer) over the three projection axes."""
    ratios = []
    for ax in range(3):
        # Project along axis `ax`: use the other two coordinates
        idx = [i for i in range(3) if i != ax]
        r2d = np.sqrt(coords[:, idx[0]]**2 + coords[:, idx[1]]**2)
        m_inner = masses[r2d <= r_inner].sum()
        m_outer = masses[r2d <= r_outer].sum()
        ratios.append(m_inner / m_outer if m_outer > 0 else np.nan)
    return np.mean(ratios), np.std(ratios), ratios


def main():
    # --- Phase B snapshots ---
    prod_dir = os.path.join(HERE, 'output_prod')
    snaps = sorted(f for f in os.listdir(prod_dir) if f.startswith('snapshot_'))
    print(f"Phase B snapshots found: {snaps}")

    new_ratios = []
    new_times = []
    for snap_name in snaps:
        path = os.path.join(prod_dir, snap_name)
        h, coords, vels, ids, masses = read_snapshot(path)
        mean, std, all_r = projected_ratio(coords, masses, R_INNER_SIM, R_OUTER_SIM)
        new_ratios.append(mean)
        new_times.append(h['time'])
        print(f"  {snap_name}: t={h['time']:.4f}, ratio={mean:.6f} +- {std:.6f}, "
              f"axes={np.array2string(np.array(all_r), precision=6)}")

    # --- Original P1 snapshots (same IC, no relaxation) ---
    orig_dir = os.path.join(DATA_ROOT, 'sim_snapshots', 'P1_elastic_control')
    orig_ratios = []
    orig_times = []
    if os.path.isdir(orig_dir):
        orig_snaps = sorted(f for f in os.listdir(orig_dir) if f.startswith('snapshot_'))
        for snap_name in orig_snaps:
            path = os.path.join(orig_dir, snap_name)
            h, coords, vels, ids, masses = read_snapshot(path)
            mean, std, _ = projected_ratio(coords, masses, R_INNER_SIM, R_OUTER_SIM)
            orig_ratios.append(mean)
            orig_times.append(h['time'])

    # --- Fluid predictions ---
    fluid_csv = os.path.join(DATA_ROOT, 'fluid_predictions_nbody_time.csv')
    fluid_times, fluid_ratios = [], []
    if os.path.exists(fluid_csv):
        with open(fluid_csv) as f:
            header = f.readline()
            for line in f:
                parts = line.strip().split(',')
                # Expect columns: t_code, ratio, ... (tolerant parse)
                try:
                    fluid_times.append(float(parts[0]))
                    fluid_ratios.append(float(parts[1]))
                except (ValueError, IndexError):
                    pass

    # --- Comparison table ---
    print("\n" + "=" * 80)
    print("P1 relaxation rerun: projected mass ratio comparison")
    print("=" * 80)
    print(f"{'t_code':<8} {'New (relaxed)':<16} {'Original':<16} {'Fluid':<16}")
    print("-" * 56)
    for i, t in enumerate(new_times):
        new_v = new_ratios[i] if i < len(new_ratios) else np.nan
        # Match original / fluid by time
        orig_v = np.interp(t, orig_times, orig_ratios) if orig_times else np.nan
        fl_v   = np.interp(t, fluid_times, fluid_ratios) if fluid_times else np.nan
        print(f"{t:<8.3f} {new_v:<16.6f} {orig_v:<16.6f} {fl_v:<16.6f}")

    # --- Relative changes ---
    if len(new_ratios) >= 3:
        print("\nRelative change (t / t0):")
        r0 = new_ratios[0]
        orig_r0 = orig_ratios[0] if orig_ratios else np.nan
        fl_r0   = fluid_ratios[0] if fluid_ratios else np.nan
        for i, t in enumerate(new_times):
            new_rel  = new_ratios[i] / r0 if r0 else np.nan
            orig_rel = (np.interp(t, orig_times, orig_ratios) / orig_r0) if orig_times else np.nan
            fl_rel   = (np.interp(t, fluid_times, fluid_ratios) / fl_r0) if fluid_times else np.nan
            print(f"  t={t:.3f}  new={new_rel:.4f}  orig={orig_rel:.4f}  fluid={fl_rel:.4f}")

        final_new  = new_ratios[-1] / r0 - 1
        final_orig = (np.interp(new_times[-1], orig_times, orig_ratios) / orig_r0 - 1) if orig_times else np.nan
        print(f"\nFinal relative change at t={new_times[-1]:.3f}:")
        print(f"  New (relaxed IC): {final_new*100:+.2f}%")
        print(f"  Original       : {final_orig*100:+.2f}%")
        print(f"  Fluid model   : -2.6% (from final_analysis_report.md)")
        print()
        if abs(final_new) < 0.05:
            print("  -> CONCLUSION: IC disequilibrium was the main cause.")
            print("     Relaxed IC drops <5%, consistent with fluid model.")
            print("     The P3 'pass' conclusion is now on firmer ground.")
        elif abs(final_new) < 0.10:
            print("  -> CONCLUSION: Partial IC effect. Some residual drop remains;")
            print("     P3 verification needs re-examination.")
        else:
            print("  -> CONCLUSION: IC disequilibrium is NOT the main cause.")
            print("     The drop persists with relaxed IC, pointing to a code/")
            print("     physics modelling issue. Revisit the P3 'pass' verdict.")

    # Save to CSV
    out_csv = os.path.join(HERE, 'relaxation_comparison.csv')
    with open(out_csv, 'w') as f:
        f.write('t_code,new_ratio,original_ratio,fluid_ratio\n')
        for i, t in enumerate(new_times):
            new_v = new_ratios[i]
            orig_v = np.interp(t, orig_times, orig_ratios) if orig_times else ''
            fl_v   = np.interp(t, fluid_times, fluid_ratios) if fluid_times else ''
            f.write(f"{t},{new_v},{orig_v},{fl_v}\n")
    print(f"\nComparison saved to {out_csv}")


if __name__ == '__main__':
    main()
