#!/usr/bin/env python3
"""Final analysis: compare N-body results with fluid model predictions.

Reads all snapshots (000, 001, 002) for P1, P2, P3 and produces a comprehensive
comparison table with the fluid model predictions.

Output:
- Console summary table
- data/P5_nbody_verify/final_comparison.csv
"""
import os, sys, numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))

from read_binary_snap import read_snapshot

# Rescaling parameters
LAMBDA = 0.085 / 3.6
R_INNER_SIM = 0.020 / LAMBDA   # 0.847 kpc
R_OUTER_SIM = 0.090 / LAMBDA   # 3.812 kpc
T_SCALE = np.sqrt(LAMBDA**3 / ((10.0/7.09e-3)*LAMBDA**3))

LOCAL_BASE = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), "data", "P5_nbody_verify", "sim_snapshots")

POINTS = ["P1_elastic_control", "P2_m3_low_sigma", "P3_m3_high_sigma"]
SNAP_NAMES = ["snapshot_000", "snapshot_001", "snapshot_002"]

# Fluid model predictions at t_code=0, 0.5, 1.0 (from compute_fluid_at_nbody_time.py)
FLUID = {
    'P1_elastic_control': {0.0: 0.146506, 0.5: 0.140624, 1.0: 0.142738},
    'P2_m3_low_sigma':    {0.0: 0.146506, 0.5: 0.147932, 1.0: 0.148193},
    'P3_m3_high_sigma':   {0.0: 0.146506, 0.5: 0.145119, 1.0: 0.141516},
}

# Reference (no SIDM, gravity only)
REF = {0.0: 0.084074, 0.5: 0.085067, 1.0: 0.086061}  # 0.5 is estimated


def compute_ratio(snap_path):
    """Compute projected mass ratio from a snapshot file."""
    if not os.path.exists(snap_path):
        return None
    
    snap = read_snapshot(snap_path, types_to_read=(1,))
    coords = snap['coords'][1]
    masses = snap['masses'][1]
    h = snap['header']
    t_code = h.get('time', 0.0)
    
    ratios = []
    for axes in [(0,1), (0,2), (1,2)]:  # xy, xz, yz
        rp = np.sqrt(coords[:,axes[0]]**2 + coords[:,axes[1]]**2)
        mi = masses[rp <= R_INNER_SIM].sum()
        mo = masses[rp <= R_OUTER_SIM].sum()
        ratios.append(mi / mo if mo > 0 else float('nan'))
    
    return {
        't_code': t_code,
        'ratio': np.mean(ratios),
        'std': np.std(ratios),
        'ratios': ratios,
        'n_particles': len(coords),
    }


def main():
    print("=" * 90)
    print("FINAL ANALYSIS: N-body vs Fluid Model Verification")
    print("=" * 90)
    print(f"R_inner_sim = {R_INNER_SIM:.4f} kpc, R_outer_sim = {R_OUTER_SIM:.4f} kpc")
    print(f"T_SCALE = {T_SCALE:.6f}")
    
    all_results = []
    
    for point in POINTS:
        print(f"\n{'='*60}")
        print(f"  {point}")
        print(f"{'='*60}")
        
        for snap_name in SNAP_NAMES:
            snap_path = os.path.join(LOCAL_BASE, point, snap_name)
            r = compute_ratio(snap_path)
            
            if r is None:
                print(f"  {snap_name}: NOT FOUND")
                continue
            
            t_code = r['t_code']
            t_phys = t_code * 0.978 * T_SCALE
            
            # Find matching fluid time
            t_key = 0.0 if t_code < 0.01 else (0.5 if abs(t_code - 0.5) < 0.1 else 1.0)
            fluid_ratio = FLUID[point].get(t_key)
            ref_ratio = REF.get(t_key)
            
            print(f"  {snap_name}: t_code={t_code:.4f}, t_phys={t_phys:.6f} Gyr")
            print(f"    N-body ratio = {r['ratio']:.6f} +/- {r['std']:.6f}")
            print(f"      (xy={r['ratios'][0]:.6f}, xz={r['ratios'][1]:.6f}, yz={r['ratios'][2]:.6f})")
            if fluid_ratio:
                print(f"    Fluid ratio  = {fluid_ratio:.6f}")
            if ref_ratio:
                print(f"    Reference    = {ref_ratio:.6f}")
            
            all_results.append({
                'point': point,
                'snapshot': snap_name,
                't_code': t_code,
                't_phys_gyr': t_phys,
                'nb_ratio': r['ratio'],
                'nb_std': r['std'],
                'fluid_ratio': fluid_ratio,
                'ref_ratio': ref_ratio,
            })
    
    # Compute relative changes
    print(f"\n{'='*90}")
    print("RELATIVE CHANGE COMPARISON")
    print(f"{'='*90}")
    
    print(f"\n{'Point':<25} {'t_code':>6} {'NB rel':>10} {'FL rel':>10} {'REF rel':>10} {'NB-FL':>10}")
    print("-" * 75)
    
    for point in POINTS:
        point_results = [r for r in all_results if r['point'] == point]
        if not point_results:
            continue
        
        init = point_results[0]
        nb_init = init['nb_ratio']
        
        for r in point_results:
            t_key = 0.0 if r['t_code'] < 0.01 else (0.5 if abs(r['t_code'] - 0.5) < 0.1 else 1.0)
            
            nb_rel = r['nb_ratio'] / nb_init if nb_init > 0 else float('nan')
            fl_init = FLUID[point][0.0]
            fl_rel = r['fluid_ratio'] / fl_init if r['fluid_ratio'] and fl_init > 0 else float('nan')
            ref_init = REF[0.0]
            ref_rel = r['ref_ratio'] / ref_init if r['ref_ratio'] and ref_init > 0 else float('nan')
            
            diff = nb_rel - fl_rel if not np.isnan(nb_rel) and not np.isnan(fl_rel) else float('nan')
            
            print(f"{point:<25} {r['t_code']:>6.2f} {nb_rel:>10.4f} {fl_rel:>10.4f} {ref_rel:>10.4f} {diff:>+10.4f}")
    
    # SIDM-only effect (remove gravity)
    print(f"\n{'='*90}")
    print("SIDM-ONLY EFFECT (N-body rel / Reference rel)")
    print(f"{'='*90}")
    
    print(f"\n{'Point':<25} {'t_code':>6} {'NB SIDM':>10} {'FL SIDM':>10} {'Verdict':>15}")
    print("-" * 70)
    
    for point in POINTS:
        point_results = [r for r in all_results if r['point'] == point]
        if not point_results:
            continue
        
        init = point_results[0]
        nb_init = init['nb_ratio']
        
        for r in point_results[1:]:  # skip t=0
            t_key = 0.5 if abs(r['t_code'] - 0.5) < 0.1 else 1.0
            
            nb_rel = r['nb_ratio'] / nb_init if nb_init > 0 else float('nan')
            fl_rel = r['fluid_ratio'] / FLUID[point][0.0] if r['fluid_ratio'] else float('nan')
            ref_rel = REF[t_key] / REF[0.0]
            
            nb_sidm = nb_rel / ref_rel if not np.isnan(nb_rel) else float('nan')
            fl_sidm = fl_rel / ref_rel if not np.isnan(fl_rel) else float('nan')
            
            # Verdict
            if not np.isnan(nb_sidm) and not np.isnan(fl_sidm):
                delta = abs(nb_sidm - fl_sidm) / abs(fl_sidm - 1) if abs(fl_sidm - 1) > 0.001 else 0
                verdict = 'PASS' if delta < 0.3 else 'MARGINAL' if delta < 0.6 else 'FAIL'
            else:
                verdict = 'N/A'
            
            print(f"{point:<25} {r['t_code']:>6.2f} {nb_sidm:>10.4f} {fl_sidm:>10.4f} {verdict:>15}")
    
    # Save to CSV
    df = pd.DataFrame(all_results)
    out_csv = os.path.join(os.path.dirname(LOCAL_BASE), "final_comparison.csv")
    df.to_csv(out_csv, index=False)
    print(f"\nSaved to {out_csv}")


if __name__ == '__main__':
    main()
