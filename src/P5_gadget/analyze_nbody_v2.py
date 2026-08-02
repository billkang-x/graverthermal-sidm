#!/usr/bin/env python3
"""
Task #42: Compare N-body results against fluid model predictions.

Reads fluid_predictions.csv (already computed) and downloaded N-body
snapshots (binary Gadget2 format with SnapFormat=1, or HDF5 with SnapFormat=3),
computes Δ_Nbody = (ratio_nbody - ratio_fluid) / ratio_fluid * 100%, and
produces a comparison plot + report.
"""
import os, sys, glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Add local module path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from read_binary_snap import read_snapshot

PROJECT_ROOT = "D:/graverthermal-sidm"
LOCAL_NBODY_DATA = os.path.join(PROJECT_ROOT, "data", "P5_nbody_verify")
FLUID_CSV = os.path.join(LOCAL_NBODY_DATA, "fluid_predictions.csv")

# Radial bin definition (in kpc, matching fluid model)
R_INNER_KPC = 0.020   # 20 pc = 0.020 kpc
R_OUTER_KPC = 0.090   # 90 pc = 0.090 kpc


def compute_nbody_ratio(snap_path, r_inner=R_INNER_KPC, r_outer=R_OUTER_KPC):
    """Read a Gadget4 snapshot (binary or HDF5) and compute M(<r_in)/M(<r_out)."""
    print(f"\n--- N-body ratio from {snap_path} ---")
    snap = read_snapshot(snap_path, types_to_read=(1,))
    h = snap['header']
    time = h.get('time', h.get('Time', 0.0))
    print(f"  Snapshot time={time}")

    if 1 not in snap.get('coords', {}):
        raise KeyError("PartType1 not in snapshot - check particle type")
    coords = snap['coords'][1]
    if 'masses' not in snap or 1 not in snap['masses']:
        raise KeyError("Masses not found for PartType1")
    masses = snap['masses'][1]

    r = np.sqrt(np.sum(coords**2, axis=1))
    print(f"  N particles: {len(r)}")
    print(f"  r range: {r.min():.4f} - {r.max():.4f} kpc")
    print(f"  Using r_inner={r_inner} kpc, r_outer={r_outer} kpc")

    m_inner = np.sum(masses[r <= r_inner])
    m_outer = np.sum(masses[r <= r_outer])

    n_inner = (r <= r_inner).sum()
    n_outer = (r <= r_outer).sum()

    if m_outer == 0:
        print(f"  WARNING: m_outer=0, no particles within {r_outer*1000:.1f} pc")
        print(f"  r_max in snapshot: {r.max():.4f}")
        return None, None, None, None

    ratio = m_inner / m_outer
    print(f"  M(<{r_inner*1000:.1f}pc) = {m_inner:.4e} ({n_inner} particles)")
    print(f"  M(<{r_outer*1000:.1f}pc) = {m_outer:.4e} ({n_outer} particles)")
    print(f"  Ratio = {ratio:.6f}")

    com = np.average(coords, axis=0, weights=masses)
    print(f"  Center of mass: {com}")

    return ratio, m_inner, m_outer, time


def find_snapshots(name):
    """Find downloaded snapshot files for a given point name.

    Looks for: snapshot_NNN (binary), snapshot_NNN.hdf5, snapdir_NNN/snapshot_NNN.N.hdf5
    Returns sorted list of paths.
    """
    d = os.path.join(LOCAL_NBODY_DATA, name)
    if not os.path.exists(d):
        return []
    snaps = []
    # Binary/HDF5 directly in dir
    for f in os.listdir(d):
        if f.startswith("snapshot_"):
            snaps.append(os.path.join(d, f))
    # snapdirs
    for sd in os.listdir(d):
        if sd.startswith("snapdir_"):
            sub_d = os.path.join(d, sd)
            for f in os.listdir(sub_d):
                if f.startswith("snapshot_"):
                    snaps.append(os.path.join(sub_d, f))
    return sorted(snaps)


def main():
    print("=" * 70)
    print("Task #42: N-body vs fluid model comparison")
    print("=" * 70)

    # 1. Load fluid predictions
    print("\n[1] Loading fluid predictions...")
    if not os.path.exists(FLUID_CSV):
        print(f"  ERROR: {FLUID_CSV} not found.")
        return
    fluid_df = pd.read_csv(FLUID_CSV)
    cols = [c for c in ['name', 'model_key', 'sigma_m_100', 'fluid_ratio_final']
            if c in fluid_df.columns]
    print(fluid_df[cols].to_string(index=False))

    # 2. Find N-body snapshots
    print("\n[2] Looking for N-body snapshots...")
    nbody_results = {}
    for _, row in fluid_df.iterrows():
        name = row['name']
        snaps = find_snapshots(name)
        if not snaps:
            print(f"  {name}: no snapshots in {os.path.join(LOCAL_NBODY_DATA, name)}")
            continue
        # Use last snapshot (final state)
        last_snap = snaps[-1]
        print(f"  {name}: {len(snaps)} snapshots, using {os.path.basename(last_snap)}")
        try:
            ratio, m_inner, m_outer, t = compute_nbody_ratio(last_snap)
            nbody_results[name] = {
                "snap_file": os.path.basename(last_snap),
                "ratio": ratio,
                "m_inner": m_inner,
                "m_outer": m_outer,
                "time": t,
            }
        except Exception as e:
            print(f"  ERROR reading snapshot: {e}")
            import traceback
            traceback.print_exc()

    # 3. Compute Δ_Nbody
    print("\n[3] Computing discrepancies...")
    results = []
    for _, row in fluid_df.iterrows():
        name = row['name']
        ratio_f = row['fluid_ratio_final']
        if name in nbody_results and nbody_results[name]['ratio'] is not None:
            ratio_n = nbody_results[name]['ratio']
            delta = (ratio_n - ratio_f) / ratio_f * 100.0
            status = "match" if abs(delta) < 10 else "discrepancy"
            print(f"  {name}: fluid={ratio_f:.4f}, nbody={ratio_n:.4f}, Δ={delta:+.2f}% [{status}]")
            results.append({
                "name": name, "model": row.get('model_key', ''),
                "sigma_m_100": row.get('sigma_m_100', ''),
                "fluid_ratio": ratio_f,
                "nbody_ratio": ratio_n, "delta_pct": delta,
                "status": status,
                "nbody_time": nbody_results[name]['time'],
                "nbody_snap": nbody_results[name]['snap_file'],
            })
        else:
            print(f"  {name}: no N-body data")
            results.append({
                "name": name, "model": row.get('model_key', ''),
                "sigma_m_100": row.get('sigma_m_100', ''),
                "fluid_ratio": ratio_f,
                "nbody_ratio": None, "delta_pct": None,
                "status": "no_nbody_data",
                "nbody_time": None,
                "nbody_snap": None,
            })

    # 4. Save CSV
    df = pd.DataFrame(results)
    out_csv = os.path.join(LOCAL_NBODY_DATA, "nbody_vs_fluid.csv")
    df.to_csv(out_csv, index=False)
    print(f"\n[4] Saved table to {out_csv}")

    # 5. Plot
    has_nbody = [r for r in results if r["nbody_ratio"] is not None]
    if has_nbody:
        fig, ax = plt.subplots(figsize=(10, 5))
        names = [r["name"] for r in has_nbody]
        deltas = [r["delta_pct"] for r in has_nbody]
        # Chinese convention: red=positive (overestimate), green=negative
        # but for clarity use green=match, red=discrepancy
        colors = ['green' if abs(d) < 10 else 'red' for d in deltas]
        bars = ax.bar(names, deltas, color=colors, alpha=0.7)
        ax.axhline(0, color='black', linewidth=0.5)
        ax.axhline(10, color='red', linewidth=1, linestyle='--', alpha=0.5, label='±10% threshold')
        ax.axhline(-10, color='red', linewidth=1, linestyle='--', alpha=0.5)
        for bar, d in zip(bars, deltas):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() * (1.01 if d>=0 else 1.05),
                    f'{d:+.1f}%', ha='center', va='bottom' if d >= 0 else 'top', fontsize=10)
        ax.set_ylabel('Δ_Nbody (%)')
        ax.set_title('N-body vs fluid model discrepancy')
        ax.legend()
        plt.xticks(rotation=20, ha='right')
        plt.tight_layout()
        out_png = os.path.join(LOCAL_NBODY_DATA, "nbody_verification.png")
        plt.savefig(out_png, dpi=150)
        print(f"[5] Saved plot to {out_png}")
    else:
        print("[5] No N-body data yet - skipping plot.")

    # 6. Summary
    print("\n[6] Summary")
    print("-" * 70)
    for r in results:
        if r["nbody_ratio"] is not None:
            print(f"  {r['name']:25s}: Δ={r['delta_pct']:+.2f}%  [{r['status']}]")
        else:
            print(f"  {r['name']:25s}: NO N-BODY DATA")
    print("-" * 70)

    n_match = sum(1 for r in results if r["status"] == "match")
    n_disc = sum(1 for r in results if r["status"] == "discrepancy")
    n_missing = sum(1 for r in results if r["status"] == "no_nbody_data")
    print(f"\nMatch: {n_match}, Discrepancy: {n_disc}, Missing: {n_missing}")

    if n_disc > 0:
        print("\n*** DISCREPANCY DETECTED ***")
        print("Points with |Δ| > 10%:")
        for r in results:
            if r["status"] == "discrepancy":
                print(f"  - {r['name']}: Δ={r['delta_pct']:+.2f}%")
        print("\nNext steps: diagnose the cooling rate implementation.")
        print("Likely causes:")
        print("  1. COOLING_PREFACTOR calibration is wrong")
        print("  2. Heat conduction coefficients (K) need adjustment")
        print("  3. Fluid model fundamentally fails at this regime")
    elif n_match > 0 and n_disc == 0 and n_missing == 0:
        print("\n*** ALL POINTS MATCH *** - cooling rate implementation verified.")


if __name__ == "__main__":
    main()
