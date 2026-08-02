"""Test with GravothermalSIDM default parameters to validate the code works."""
import sys, os, time
sys.path.insert(0, r'D:/graverthermal-sidm/external/gravothermalsidm')

import numpy as np
from SourcePy.evolve import Halo
from SourcePy.record import HaloRecord

# Use the DEFAULT GravothermalSIDM parameters (from runHaloEvolution.py)
dir_data = r'D:/graverthermal-sidm/data/P2_runs/_test_default'
import shutil
if os.path.exists(dir_data):
    shutil.rmtree(dir_data)

halorec = HaloRecord(dir_data)

# Default: r_s=2.586, rho_s=0.0194, sigma_m=5, w=1, n_shells=400
haloevo = Halo(halorec, flag_hydrostatic_initial=True)

print(f"Default halo:")
print(f"  r_s = {haloevo.r_s} kpc")
print(f"  rho_s = {haloevo.rho_s} Msun/pc^3")
print(f"  sigma_m = {haloevo.sigma_m_with_units} cm^2/g")
print(f"  w = {haloevo.w_units} km/s")
print(f"  n_shells = {haloevo.n_shells}")
print(f"  scale_t (Gyr) = {haloevo.scale_t.to('Gyr').value:.6f}")
print(f"  t_relax (dimensionless) = {haloevo.t_relax:.4f}")
print(f"  t_relax (Gyr) = {(haloevo.t_relax * haloevo.scale_t).to('Gyr').value:.6f}")

# Run with default settings from runHaloEvolution.py
t0 = time.time()
haloevo.evolve_halo(t_end=25.0, save_frequency_rate=10)
elapsed = time.time() - t0

print(f"\nEvolution completed in {elapsed:.1f}s")
print(f"Final t = {haloevo.t:.4f} (dimensionless)")
print(f"Number of conduction steps: {haloevo.n_conduction}")

list_files, list_times = halorec.glob_pickle_files()
print(f"Saved {len(list_files)} snapshots")
print(f"Time range: {list_times[0]:.4f} to {list_times[-1]:.4f}")

# Analyze
times = []
rho_center = []
for f in list_files:
    data = halorec.get_halo_state_pickled(file_halo=f)
    times.append(data['t'])
    rho_center.append(data['rho'][3])

times = np.array(times)
rho_center = np.array(rho_center)

print(f"\nCentral density evolution:")
print(f"  Initial: rho = {rho_center[0]:.4e}")
idx_min = np.argmin(rho_center)
print(f"  Min rho = {rho_center[idx_min]:.4e} at t = {times[idx_min]:.4f}")
print(f"  Final: rho = {rho_center[-1]:.4e} at t = {times[-1]:.4f}")
print(f"  rho_final/rho_initial = {rho_center[-1]/rho_center[0]:.4e}")
