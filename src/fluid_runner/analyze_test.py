"""Analyze the elastic test snapshots we already have."""
import sys, os
sys.path.insert(0, r'D:/graverthermal-sidm/external/gravothermalsidm')

import numpy as np
from astropy import units as ut
from SourcePy.evolve import Halo
from SourcePy.record import HaloRecord

dir_data = r'D:/graverthermal-sidm/data/P2_runs/_test_fast'
halorec = HaloRecord(dir_data)

list_files, list_times = halorec.glob_pickle_files()
print(f"Snapshots: {len(list_files)}")
print(f"Time range: {list_times[0]:.6f} to {list_times[-1]:.6f} (dimensionless)")

# Get the halo initialization to access scales
halo_ini, orig = halorec.get_halo_initialization()
r_s = halo_ini['r_s']  # kpc
rho_s = halo_ini['rho_s']  # Msun/pc^3
sigma_m = halo_ini['sigma_m_with_units']
w = halo_ini['w_units']

# Recompute scales
from astropy import constants as ct
scale_r = r_s * ut.kpc
scale_rho = rho_s * ut.M_sun/ut.pc**3
scale_m = 4*np.pi * scale_rho * scale_r**3
scale_u = ct.G * scale_m / scale_r
scale_v = np.sqrt(scale_u)
scale_t = 1./np.sqrt(4.*np.pi*scale_rho * ct.G)

scale_t_gyr = scale_t.to('Gyr').value
scale_rho_mspc = scale_rho.to('Msun/pc**3').value
scale_v_kms = scale_v.to('km/s').value
scale_r_kpc = scale_r.to('kpc').value

print(f"\nScales:")
print(f"  scale_t = {scale_t_gyr:.6f} Gyr")
print(f"  scale_rho = {scale_rho_mspc:.6e} Msun/pc^3")
print(f"  scale_v = {scale_v_kms:.2f} km/s")
print(f"  scale_r = {scale_r_kpc:.2f} kpc")

# Load and analyze snapshots
times_gyr = []
rho_center = []
v_center = []
for f in list_files:
    data = halorec.get_halo_state_pickled(file_halo=f)
    times_gyr.append(data['t'] * scale_t_gyr)
    rho_center.append(data['rho'][3] * scale_rho_mspc)
    v_center.append(np.sqrt(data['p'][3]/data['rho'][3]) * scale_v_kms)

times_gyr = np.array(times_gyr)
rho_center = np.array(rho_center)
v_center = np.array(v_center)

print(f"\nCentral density evolution (first/last 5):")
for i in [0, 1, 2, 3, 4, -5, -4, -3, -2, -1]:
    if i < len(times_gyr):
        print(f"  t={times_gyr[i]:.6f} Gyr: rho={rho_center[i]:.4e} Msun/pc^3, v={v_center[i]:.2f} km/s")

# Find minimum density (core formation)
idx_min = np.argmin(rho_center)
print(f"\nCore formation (min density):")
print(f"  t_core = {times_gyr[idx_min]:.6f} Gyr")
print(f"  rho_min = {rho_center[idx_min]:.4e} Msun/pc^3")

# Final state
print(f"\nFinal state:")
print(f"  t = {times_gyr[-1]:.6f} Gyr (dimensionless: {list_times[-1]:.4f})")
print(f"  rho = {rho_center[-1]:.4e} Msun/pc^3")
print(f"  v = {v_center[-1]:.2f} km/s")
print(f"  rho/rho0 = {rho_center[-1]/rho_center[0]:.4e}")

# Check if we're in collapse phase
if rho_center[-1] > rho_center[idx_min]:
    print(f"\n  -> Collapse phase detected (rho increasing from minimum)")
    print(f"  -> rho_final/rho_min = {rho_center[-1]/rho_center[idx_min]:.4f}")
