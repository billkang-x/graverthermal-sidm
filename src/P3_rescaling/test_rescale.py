"""Test the rescaling module with a single snapshot.

This test loads the elastic P2 run, scans r_2D/r_s to find where the projected
mass ratio matches B1938+666 observations, and computes the rescaled parameters
for the matching points.
"""
import sys, os

# Ensure we can import the rescale module (which itself sets up paths to
# external/gravothermalsidm and src/cross_sections).
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import numpy as np
from astropy import units as ut
from astropy import constants as ct

# rescale.py handles all sys.path setup for SourcePy/record and dsidm_models
from rescale import (
    projected_enclosed_mass,
    find_matching_radii,
    compute_rescaling,
    rescale_parameters,
    M_20PC, M_90PC, MASS_RATIO_OBS, MASS_RATIO_ERR,
    RS_SIM_KPC, RHO0_SIM, SIGMA_M_SIM, T_ZOBS_GYR,
)
from SourcePy.record import HaloRecord


def get_scales(halo_ini):
    """Compute dimensionful scales from halo initialization dict."""
    r_s = halo_ini['r_s'] * ut.kpc
    rho_s = halo_ini['rho_s'] * ut.M_sun / ut.pc**3
    scale_r_kpc = r_s.to('kpc').value
    scale_rho = rho_s.to('Msun/pc**3').value
    scale_v = np.sqrt(ct.G * 4 * np.pi * rho_s * r_s**3 / r_s).to('km/s').value
    scale_t = (1.0 / np.sqrt(4.0 * np.pi * rho_s * ct.G)).to('Gyr').value
    return scale_r_kpc, scale_rho, scale_v, scale_t


def main():
    # ---- Load elastic snapshots ----
    dir_data = r'D:/graverthermal-sidm/data/P2_runs/elastic'
    if not os.path.isdir(dir_data):
        print(f"ERROR: P2 elastic data dir not found: {dir_data}")
        return

    halorec = HaloRecord(dir_data)
    list_files, list_times = halorec.glob_pickle_files()
    print(f"Elastic snapshots: {len(list_files)}")
    if len(list_files) == 0:
        print("No snapshots to test.")
        return

    # ---- Get scales ----
    halo_ini, _ = halorec.get_halo_initialization()
    scale_r_kpc, scale_rho, scale_v, scale_t = get_scales(halo_ini)
    print(f"scale_r = {scale_r_kpc:.4f} kpc, scale_rho = {scale_rho:.4e}, "
          f"scale_v = {scale_v:.4f} km/s, scale_t = {scale_t:.6f} Gyr")

    # ---- Test with the last snapshot ----
    data = halorec.get_halo_state_pickled(file_halo=list_files[-1])
    r_kpc = data['r'] * scale_r_kpc
    rho = data['rho'] * scale_rho
    t_gyr = data['t'] * scale_t
    print(f"\nLast snapshot: t = {t_gyr:.4f} Gyr")
    print(f"  r range: {r_kpc.min():.4e} to {r_kpc.max():.4e} kpc")
    print(f"  rho range: {rho.min():.4e} to {rho.max():.4e} Msun/pc^3")

    # ---- Test projected mass at 20 and 90 pc ----
    r_test_20 = 20.0 / 1000.0   # 20 pc in kpc
    r_test_90 = 90.0 / 1000.0   # 90 pc in kpc
    M_2d_20 = projected_enclosed_mass(r_kpc, rho, r_test_20)
    M_2d_90 = projected_enclosed_mass(r_kpc, rho, r_test_90)
    print(f"\n  M_2D(r < 20 pc) = {M_2d_20:.4e} Msun  (obs: {M_20PC:.4e})")
    print(f"  M_2D(r < 90 pc) = {M_2d_90:.4e} Msun  (obs: {M_90PC:.4e})")
    print(f"  Ratio: {M_2d_20/M_2d_90:.4f} (observed: {MASS_RATIO_OBS:.4f})")

    # ---- Scan r_2D/r_s ----
    print("\n--- Scanning r_2D/r_s ---")
    scan = find_matching_radii(r_kpc, rho, RS_SIM_KPC, n_scan=100)
    print(f"  {'r_2D/r_s':<12} {'M_inner':<14} {'M_outer':<14} {'ratio':<10} {'n_sigma':<10}")
    for row in scan[::10]:
        r2D_rs, M_in, M_out, ratio = row
        n_sigma = abs(ratio - MASS_RATIO_OBS) / MASS_RATIO_ERR if ratio > 0 else np.nan
        print(f"  {r2D_rs:<12.4f} {M_in:<14.4e} {M_out:<14.4e} "
              f"{ratio:<10.4f} {n_sigma:<10.2f}")

    # ---- Find matches within 3σ ----
    matches_mask = np.array([
        (row[3] > 0) and (abs(row[3] - MASS_RATIO_OBS) / MASS_RATIO_ERR <= 3)
        for row in scan
    ])
    print(f"\n  Matching points (3σ): {int(np.sum(matches_mask))}")

    if np.any(matches_mask):
        print(f"\n  {'r_2D/r_s':<10} {'λ':<12} {'μ':<12} {'r_s':<10} "
              f"{'ρ_0':<14} {'σ/m':<10} {'t_evo':<10} {'t_ok'}")
        for row in scan[matches_mask][:8]:
            r2D_rs, M_in, M_out, ratio = row
            lam, mu = compute_rescaling(r2D_rs, M_in, M_out)
            params = rescale_parameters(lam, mu, t_gyr)
            t_ok = 'YES' if params['t_evo_gyr'] <= T_ZOBS_GYR else 'no'
            print(f"  {r2D_rs:<10.4f} {lam:<12.4e} {mu:<12.4e} "
                  f"{params['r_s_kpc']:<10.4f} {params['rho0_msun_pc3']:<14.4e} "
                  f"{params['sigma_m_cm2_g']:<10.2f} {params['t_evo_gyr']:<10.2f} {t_ok}")


if __name__ == '__main__':
    main()
