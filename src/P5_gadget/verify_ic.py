#!/usr/bin/env python3
"""Verify IC properties against fluid model expectations.

The IC was generated from halo_ini.h5 (the initial state from the fluid model).
This script checks:
  - The radial mass profile M(<r) matches the fluid model's initial state
  - The mass ratio M(20pc)/M(90pc) matches fluid_predictions.csv's ratio_init
"""
import os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from read_binary_snap import read_snapshot

LOCAL_NBODY_DATA = "D:/graverthermal-sidm/data/P5_nbody_verify"
IC_PATH = os.path.join(LOCAL_NBODY_DATA, "ics/ic.dat")
FLUID_CSV = os.path.join(LOCAL_NBODY_DATA, "fluid_predictions.csv")

R_INNER_KPC = 0.020
R_OUTER_KPC = 0.090


def main():
    print("=" * 60)
    print("IC verification against fluid model")
    print("=" * 60)

    # Load IC
    print(f"\nReading IC: {IC_PATH}")
    snap = read_snapshot(IC_PATH, types_to_read=(1,))
    coords = snap['coords'][1]
    masses = snap['masses'][1]
    r = np.sqrt(np.sum(coords**2, axis=1))
    print(f"N = {len(r)}, M_total = {masses.sum():.4e}")
    print(f"r: min={r.min():.5f}, max={r.max():.5f}, mean={r.mean():.5f} kpc")

    # Compute cumulative mass profile
    radii = np.linspace(0.001, 0.4, 100)
    cum_mass = np.array([masses[r <= rr].sum() for rr in radii])

    # Mass ratio at R_INNER, R_OUTER
    m_inner = masses[r <= R_INNER_KPC].sum()
    m_outer = masses[r <= R_OUTER_KPC].sum()
    ratio = m_inner / m_outer
    print(f"\nM(<{R_INNER_KPC*1000:.1f}pc) = {m_inner:.4e} ({(r <= R_INNER_KPC).sum()} particles)")
    print(f"M(<{R_OUTER_KPC*1000:.1f}pc) = {m_outer:.4e} ({(r <= R_OUTER_KPC).sum()} particles)")
    print(f"Ratio = {ratio:.6f}")

    # Compare to fluid model
    df = pd.read_csv(FLUID_CSV)
    print("\nFluid predictions:")
    for _, row in df.iterrows():
        print(f"  {row['name']:25s}: ratio_init={row['fluid_ratio_init']:.6f}, "
              f"ratio_final={row['fluid_ratio_final']:.6f}")

    # Find P1 (elastic control) - its ratio_init should match IC
    p1 = df[df['name'] == 'P1_elastic_control'].iloc[0]
    print(f"\nP1 ratio_init = {p1['fluid_ratio_init']:.6f}")
    print(f"IC ratio      = {ratio:.6f}")
    delta = (ratio - p1['fluid_ratio_init']) / p1['fluid_ratio_init'] * 100
    print(f"Δ = {delta:+.2f}%")

    if abs(delta) < 1.0:
        print("✓ IC matches fluid initial state (within 1%)")
    else:
        print("✗ WARNING: IC does not match fluid initial state!")
        print(f"  Possible unit mismatch. Check coordinate units.")
        print(f"  r range in IC: {r.min():.4f} - {r.max():.4f}")
        print(f"  Expected: ~0.001 - ~0.3 kpc (1pc - 300pc)")

    # Mass profile table
    print("\nCumulative mass profile:")
    print(f"{'r (kpc)':>10} {'r (pc)':>10} {'M(<r) (Msun)':>15} {'N':>8}")
    print("-" * 50)
    for rr in [0.005, 0.010, 0.020, 0.030, 0.050, 0.090, 0.150, 0.250, 0.340]:
        m = masses[r <= rr].sum()
        n = (r <= rr).sum()
        print(f"{rr:>10.4f} {rr*1000:>10.1f} {m:>15.4e} {n:>8d}")


if __name__ == "__main__":
    main()
