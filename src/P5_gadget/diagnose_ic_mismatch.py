#!/usr/bin/env python3
"""Diagnose IC vs fluid model mismatch.

Computes the projected mass ratio for:
1. Current IC (N-body, R_MAX=0.34 kpc)
2. NFW analytic with various R_MAX cutoffs
3. Fluid model initial state (halo_ini.h5)

This tells us whether the mismatch is due to:
  (a) R_MAX being too small (missing outer halo)
  (b) Different density profile (NFW vs hydrostatic-adjusted)
  (c) Projection method difference
"""
import os, sys
import numpy as np
import h5py

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from read_binary_snap import read_snapshot
from projected_mass import projected_mass_ratio

R_INNER = 0.020  # kpc = 20 pc
R_OUTER = 0.090  # kpc = 90 pc

# NFW parameters
R_S = 0.085  # kpc
RHO_0 = 10.0  # Msun/pc^3


def rho_nfw(r_pc):
    x = r_pc / (R_S * 1000.0)
    return RHO_0 / (x * (1+x)**2)


def analytic_projected_mass_ratio(r_max_kpc):
    """Compute projected mass ratio for NFW truncated at r_max."""
    trapz = np.trapezoid
    r_max_pc = r_max_kpc * 1000.0
    r_in_pc = R_INNER * 1000.0
    r_out_pc = R_OUTER * 1000.0

    def M_2d(R_target_pc):
        # Part 1: interior (r < R_target)
        r_arr = np.linspace(0.5, R_target_pc, 5001)
        rho_arr = rho_nfw(r_arr)
        M_in = 2 * np.pi * trapz(r_arr * rho_arr, r_arr)
        # Part 2: exterior (r > R_target, up to r_max)
        r_arr2 = np.linspace(R_target_pc * (1+1e-6), r_max_pc, 50001)
        rho_arr2 = rho_nfw(r_arr2)
        integrand = 2 * r_arr2 * rho_arr2 / np.sqrt(r_arr2**2 - R_target_pc**2)
        M_out = R_target_pc * trapz(integrand, r_arr2)
        return M_in + M_out

    m_in = M_2d(r_in_pc)
    m_out = M_2d(r_out_pc)
    return m_in / m_out if m_out > 0 else np.nan, m_in, m_out


def main():
    print("=" * 70)
    print("IC vs fluid model mismatch diagnosis")
    print("=" * 70)

    # 1. Current IC (N-body)
    print("\n[1] Current IC (N-body, R_MAX=0.34 kpc):")
    ic_path = "D:/graverthermal-sidm/data/P5_nbody_verify/ics/ic.dat"
    snap = read_snapshot(ic_path, types_to_read=(1,))
    coords = snap['coords'][1]
    masses = snap['masses'][1]
    r3d = np.sqrt((coords**2).sum(axis=1))
    print(f"  N={len(r3d)}, r_max={r3d.max():.4f} kpc, M_total={masses.sum():.4e}")

    # 3D spherical ratio
    m3d_in = masses[r3d <= R_INNER].sum()
    m3d_out = masses[r3d <= R_OUTER].sum()
    print(f"  3D: M(20pc)={m3d_in:.4e}, M(90pc)={m3d_out:.4e}, ratio={m3d_in/m3d_out:.6f}")

    # 2D projected ratio (3-axis avg)
    ratio_2d, m2d_in, m2d_out = projected_mass_ratio(coords, masses, R_INNER, R_OUTER)
    print(f"  2D: M(20pc)={m2d_in:.4e}, M(90pc)={m2d_out:.4e}, ratio={ratio_2d:.6f}")

    # 2. Analytic NFW with various R_MAX
    print("\n[2] Analytic NFW projected mass ratio vs R_MAX:")
    print(f"  {'R_MAX (kpc)':>12} {'ratio':>10} {'M(20pc)':>12} {'M(90pc)':>12}")
    for r_max in [0.34, 0.5, 1.0, 2.0, 4.25, 10.0, 50.0, 200.0]:
        r, m_in, m_out = analytic_projected_mass_ratio(r_max)
        print(f"  {r_max:>12.2f} {r:>10.6f} {m_in:>12.4e} {m_out:>12.4e}")

    # 3. Fluid model initial state (halo_ini.h5 - raw NFW before hydrostatic adjustment)
    print("\n[3] Fluid model halo_ini.h5 (raw NFW, code units 0.02-50):")
    halo_path = "D:/graverthermal-sidm/data/P5_nbody_verify/_tmp/P1_elastic_control/halo_ini.h5"
    if os.path.exists(halo_path):
        with h5py.File(halo_path, 'r') as f:
            r_code = f['r'][:]
            rho_code = f['rho'][:]
        # Scale: r_s = 0.085 kpc, rho_s = 10 Msun/pc^3
        r_kpc = r_code * R_S
        rho_msun_pc3 = rho_code * RHO_0
        print(f"  r range: {r_kpc.min():.5f} - {r_kpc.max():.5f} kpc")
        print(f"  rho range: {rho_msun_pc3.min():.4e} - {rho_msun_pc3.max():.4e}")

        # Use rescale.py projected_enclosed_mass
        sys.path.insert(0, 'D:/graverthermal-sidm/src/P3_rescaling')
        from rescale import projected_enclosed_mass
        M_in = projected_enclosed_mass(r_kpc, rho_msun_pc3, R_INNER,
                                       r_unit='kpc', rho_unit='Msun_pc3')
        M_out = projected_enclosed_mass(r_kpc, rho_msun_pc3, R_OUTER,
                                        r_unit='kpc', rho_unit='Msun_pc3')
        print(f"  projected_enclosed_mass: M(20pc)={M_in:.4e}, M(90pc)={M_out:.4e}")
        print(f"  Ratio = {M_in/M_out:.6f}")

    # 4. Summary
    print("\n[4] Summary:")
    print(f"  IC 2D ratio:        {ratio_2d:.6f}")
    r_analytic, _, _ = analytic_projected_mass_ratio(4.25)
    print(f"  NFW analytic (R_MAX=4.25 kpc, matching fluid): {r_analytic:.6f}")
    r_analytic_034, _, _ = analytic_projected_mass_ratio(0.34)
    print(f"  NFW analytic (R_MAX=0.34 kpc, matching IC):   {r_analytic_034:.6f}")
    print(f"\n  If IC ratio ≈ NFW(0.34), the IC is correct for its R_MAX.")
    print(f"  The mismatch with fluid is because IC has smaller R_MAX.")
    print(f"  → Need to regenerate IC with R_MAX ≥ 4.25 kpc to match fluid halo.")


if __name__ == "__main__":
    main()
