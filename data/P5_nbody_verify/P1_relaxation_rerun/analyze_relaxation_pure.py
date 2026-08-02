#!/usr/bin/env python3
"""Analyze P1 relaxation rerun - pure Python (no numpy).

Computes projected enclosed mass ratio M_2D(<R_inner)/M_2D(<R_outer)
for each Phase B snapshot and compares with original P1 + fluid model.
"""
import struct
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_ROOT = os.path.normpath(os.path.join(HERE, '..'))

R_INNER_SIM = 0.847   # kpc
R_OUTER_SIM = 3.812   # kpc

HEADER_SIZE = 256


def read_block(f):
    sz = struct.unpack('I', f.read(4))[0]
    data = f.read(sz)
    f.read(4)
    return data


def parse_header(hdr_buf):
    fmt = '<6I6d2d2i6I2i4d2i6Ii'
    sz = struct.calcsize(fmt)
    vals = struct.unpack(fmt, hdr_buf[:sz])
    h = {}
    h['npart'] = list(vals[0:6])
    h['mass'] = list(vals[6:12])
    h['time'] = vals[12]
    return h


def read_snapshot(path):
    """Read PartType1 coords/masses from Gadget4 binary snapshot."""
    with open(path, 'rb') as f:
        hdr_buf = read_block(f)
        h = parse_header(hdr_buf)
        n1 = h['npart'][1]
        n2 = h['npart'][2]
        m1 = h['mass'][1]
        m2 = h['mass'][2]

        # Read blocks: Coordinates, Velocities, ParticleIDs
        coords_blk = read_block(f)
        vels_blk = read_block(f)
        ids_blk = read_block(f)

        # Parse coordinates as float32 (N_total, 3)
        n_total = n1 + n2
        import array
        coords = array.array('f')
        coords.frombytes(coords_blk)
        # Reshape: list of (x,y,z) tuples
        coords_list = [(coords[i*3], coords[i*3+1], coords[i*3+2]) for i in range(n_total)]

        # Masses: use header mass table
        masses = []
        for i in range(n_total):
            if i < n1:
                masses.append(m1)
            else:
                masses.append(m2)

    return h, coords_list, masses, n1, n2


def projected_ratio(coords, masses, r_inner, r_outer):
    """Average M_2D(<r_inner)/M_2D(<r_outer) over 3 projection axes."""
    ratios = []
    for ax in range(3):
        # Project along axis ax: use the other two coordinates
        idx = [i for i in range(3) if i != ax]
        m_inner = 0.0
        m_outer = 0.0
        for i, (c) in enumerate(coords):
            r2d = (c[idx[0]]**2 + c[idx[1]]**2) ** 0.5
            if r2d <= r_outer:
                m_outer += masses[i]
                if r2d <= r_inner:
                    m_inner += masses[i]
        if m_outer > 0:
            ratios.append(m_inner / m_outer)
        else:
            ratios.append(float('nan'))
    mean = sum(ratios) / len(ratios)
    return mean, ratios


def main():
    # Phase B snapshots
    prod_dir = os.path.join(HERE, 'output_prod')
    snaps = sorted(f for f in os.listdir(prod_dir) if f.startswith('snapshot_'))
    print(f"Phase B snapshots: {snaps}")

    new_ratios = []
    new_times = []
    for snap_name in snaps:
        path = os.path.join(prod_dir, snap_name)
        print(f"\nReading {snap_name}...", flush=True)
        h, coords, masses, n1, n2 = read_snapshot(path)
        print(f"  time={h['time']:.4f}, n1={n1}, n2={n2}", flush=True)
        mean, all_r = projected_ratio(coords, masses, R_INNER_SIM, R_OUTER_SIM)
        new_ratios.append(mean)
        new_times.append(h['time'])
        print(f"  ratio={mean:.6f}, axes={[round(r,6) for r in all_r]}")

    # Original P1 snapshots
    orig_dir = os.path.join(DATA_ROOT, 'sim_snapshots', 'P1_elastic_control')
    orig_ratios = []
    orig_times = []
    if os.path.isdir(orig_dir):
        orig_snaps = sorted(f for f in os.listdir(orig_dir) if f.startswith('snapshot_'))
        for snap_name in orig_snaps:
            path = os.path.join(orig_dir, snap_name)
            print(f"\nReading original {snap_name}...", flush=True)
            h, coords, masses, n1, n2 = read_snapshot(path)
            mean, all_r = projected_ratio(coords, masses, R_INNER_SIM, R_OUTER_SIM)
            orig_ratios.append(mean)
            orig_times.append(h['time'])
            print(f"  time={h['time']:.4f}, ratio={mean:.6f}")

    # Fluid predictions (from report)
    fluid_data = [
        (0.0, 0.146506),
        (0.5, 0.140624),
        (1.0, 0.142738),
    ]
    fluid_times = [d[0] for d in fluid_data]
    fluid_ratios = [d[1] for d in fluid_data]

    # Comparison table
    print("\n" + "=" * 80)
    print("P1 relaxation rerun: projected mass ratio comparison")
    print("=" * 80)
    print(f"{'t_code':<8} {'New (relaxed)':<16} {'Original':<16} {'Fluid':<16}")
    print("-" * 56)
    for i, t in enumerate(new_times):
        new_v = new_ratios[i]
        # Match original by time
        orig_v = float('nan')
        for j, ot in enumerate(orig_times):
            if abs(ot - t) < 0.01:
                orig_v = orig_ratios[j]
                break
        fl_v = float('nan')
        for j, ft in enumerate(fluid_times):
            if abs(ft - t) < 0.01:
                fl_v = fluid_ratios[j]
                break
        print(f"{t:<8.3f} {new_v:<16.6f} {orig_v:<16.6f} {fl_v:<16.6f}")

    # Relative changes
    if len(new_ratios) >= 3 and len(orig_ratios) >= 3:
        print("\nRelative change (ratio / ratio_t0):")
        r0_new = new_ratios[0]
        r0_orig = orig_ratios[0]
        r0_fl = fluid_ratios[0]
        for i, t in enumerate(new_times):
            new_rel = new_ratios[i] / r0_new if r0_new else float('nan')
            orig_rel = float('nan')
            for j, ot in enumerate(orig_times):
                if abs(ot - t) < 0.01:
                    orig_rel = orig_ratios[j] / r0_orig
                    break
            fl_rel = float('nan')
            for j, ft in enumerate(fluid_times):
                if abs(ft - t) < 0.01:
                    fl_rel = fluid_ratios[j] / r0_fl
                    break
            print(f"  t={t:.3f}  new={new_rel:.4f}  orig={orig_rel:.4f}  fluid={fl_rel:.4f}")

        final_new = new_ratios[-1] / r0_new - 1
        final_orig = orig_ratios[-1] / r0_orig - 1
        final_fl = fluid_ratios[-1] / r0_fl - 1
        print(f"\nFinal relative change at t={new_times[-1]:.3f}:")
        print(f"  New (relaxed IC): {final_new*100:+.2f}%")
        print(f"  Original       : {final_orig*100:+.2f}%")
        print(f"  Fluid model    : {final_fl*100:+.2f}%")
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

    # Save CSV
    out_csv = os.path.join(HERE, 'relaxation_comparison.csv')
    with open(out_csv, 'w') as f:
        f.write('t_code,new_ratio,original_ratio,fluid_ratio\n')
        for i, t in enumerate(new_times):
            new_v = new_ratios[i]
            orig_v = ''
            for j, ot in enumerate(orig_times):
                if abs(ot - t) < 0.01:
                    orig_v = orig_ratios[j]
                    break
            fl_v = ''
            for j, ft in enumerate(fluid_times):
                if abs(ft - t) < 0.01:
                    fl_v = fluid_ratios[j]
                    break
            f.write(f"{t},{new_v},{orig_v},{fl_v}\n")
    print(f"\nComparison saved to {out_csv}")


if __name__ == '__main__':
    main()
