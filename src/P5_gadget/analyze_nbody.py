"""
Task #42: Compare N-body results against fluid model predictions.

Reads fluid_predictions.csv (already computed) and the downloaded N-body
snapshots, computes Δ_Nbody, and produces a comparison plot + report.
"""
import os
import sys
import numpy as np
import pandas as pd
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PROJECT_ROOT = "D:/graverthermal-sidm"
LOCAL_NBODY_DATA = os.path.join(PROJECT_ROOT, "data", "P5_nbody_verify")
FLUID_CSV = os.path.join(LOCAL_NBODY_DATA, "fluid_predictions.csv")

R_INNER_PC = 20.0
R_OUTER_PC = 90.0


def compute_nbody_ratio(snap_path, r_inner_pc=R_INNER_PC, r_outer_pc=R_OUTER_PC):
    """Read a Gadget4 snapshot and compute the mass ratio M(<r_in)/M(<r_out)."""
    print(f"\n--- N-body ratio from {snap_path} ---")
    with h5py.File(snap_path, 'r') as f:
        h = f['Header']
        time = h.attrs['Time']
        box = h.attrs['BoxSize']
        mass_table = h.attrs['MassTable']
        n_total = h.attrs['NumPart_Total']

        print(f"  Snapshot time={time}, box={box}")
        print(f"  NumPart_Total={n_total}, MassTable={mass_table}")

        if 'PartType1' not in f:
            raise KeyError("PartType1 not in snapshot - check particle type")
        p1 = f['PartType1']
        coords = p1['Coordinates'][:]
        masses = p1['Masses'][:]

        # If coordinates are in code units, convert using UnitLength
        # Check the units via attributes
        # The IC was written with kpc units directly.
        # Gadget4 may rescale to code units (BoxSize-based); need to check.
        # For safety, assume coords are in the same units as IC (kpc).
        r = np.sqrt(np.sum(coords**2, axis=1))

        r_inner_kpc = r_inner_pc / 1000.0
        r_outer_kpc = r_outer_pc / 1000.0

        m_inner = np.sum(masses[r <= r_inner_kpc])
        m_outer = np.sum(masses[r <= r_outer_kpc])

        n_inner = (r <= r_inner_kpc).sum()
        n_outer = (r <= r_outer_kpc).sum()

        if m_outer == 0:
            print(f"  WARNING: m_outer=0, no particles within {r_outer_pc} pc")
            # Try larger radii
            r_max = r.max()
            print(f"  r_max in snapshot: {r_max:.4f}")
            print(f"  (Note: if r_max << {r_outer_pc/1000:.4f} kpc, units are wrong)")
            return None, None, None, None

        ratio = m_inner / m_outer
        print(f"  M(<{r_inner_pc}pc) = {m_inner:.4e} Msun ({n_inner} particles)")
        print(f"  M(<{r_outer_pc}pc) = {m_outer:.4e} Msun ({n_outer} particles)")
        print(f"  Ratio = {ratio:.6f}")

        com = np.average(coords, axis=0, weights=masses)
        print(f"  Center of mass: {com}")

    return ratio, m_inner, m_outer, time


def main():
    print("=" * 70)
    print("Task #42: N-body vs fluid model comparison")
    print("=" * 70)

    # 1. Load fluid predictions
    print("\n[1] Loading fluid predictions...")
    if not os.path.exists(FLUID_CSV):
        print(f"  ERROR: {FLUID_CSV} not found. Run compute_fluid_predictions.py first.")
        return
    fluid_df = pd.read_csv(FLUID_CSV)
    print(fluid_df[['name', 'model_key', 'sigma_m_100', 'fluid_ratio_final']].to_string(index=False))

    # 2. Find N-body snapshots
    print("\n[2] Looking for N-body snapshots...")
    nbody_results = {}
    for _, row in fluid_df.iterrows():
        name = row['name']
        snap_dir = os.path.join(LOCAL_NBODY_DATA, name)
        if not os.path.exists(snap_dir):
            print(f"  {name}: no local data")
            continue
        snaps = sorted([f for f in os.listdir(snap_dir) if f.startswith('snap_') and f.endswith('.hdf5')])
        if not snaps:
            print(f"  {name}: no snapshots in {snap_dir}")
            continue
        last_snap = snaps[-1]
        snap_path = os.path.join(snap_dir, last_snap)
        print(f"  {name}: using {last_snap}")
        try:
            ratio, m_inner, m_outer, t = compute_nbody_ratio(snap_path)
            nbody_results[name] = {
                "snap_file": last_snap,
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
                "name": name, "model": row['model_key'],
                "sigma_m_100": row['sigma_m_100'],
                "fluid_ratio": ratio_f,
                "nbody_ratio": ratio_n, "delta_pct": delta,
                "status": status,
                "nbody_time": nbody_results[name]['time'],
            })
        else:
            print(f"  {name}: no N-body data")
            results.append({
                "name": name, "model": row['model_key'],
                "sigma_m_100": row['sigma_m_100'],
                "fluid_ratio": ratio_f,
                "nbody_ratio": None, "delta_pct": None,
                "status": "no_nbody_data",
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
        colors = ['green' if abs(d) < 10 else 'red' for d in deltas]
        bars = ax.bar(names, deltas, color=colors, alpha=0.7)
        ax.axhline(0, color='black', linewidth=0.5)
        ax.axhline(10, color='red', linewidth=1, linestyle='--', alpha=0.5, label='±10% threshold')
        ax.axhline(-10, color='red', linewidth=1, linestyle='--', alpha=0.5)
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


if __name__ == "__main__":
    main()
