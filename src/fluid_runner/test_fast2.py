"""
Test elastic halo with parameters optimized for fast evolution.

Key insight: The gravothermal evolution has a rescaling symmetry.
Schmidt et al. 2026 Appendix G shows that (lambda, mu) rescaling preserves
the evolution shape. We can use larger sigma_m (faster evolution) and
rescale the times afterward.

For validation, we just need to see the qualitative behavior:
- Core formation (density drops then rises)
- Collapse (density increases sharply)

The absolute times will differ but the shape should match.
"""
import sys, os, time
sys.path.insert(0, r'D:/graverthermal-sidm/external/gravothermalsidm')

import numpy as np
from astropy import units as ut
from SourcePy.evolve import Halo
from SourcePy.record import HaloRecord

dir_data = r'D:/graverthermal-sidm/data/P2_runs/_test_fast2'
import shutil
if os.path.exists(dir_data):
    shutil.rmtree(dir_data)

halorec = HaloRecord(dir_data)

# Use parameters that give fast evolution:
# - Larger sigma_m = more self-interactions = faster evolution
# - Fewer shells = faster per step
# - Use relaxation time step only (no energy criterion)
# - Larger t_epsilon = bigger steps
haloevo = Halo(halorec,
               profile='NFW',
               r_s=3.6,
               rho_s=7.09e-3,
               sigma_m_with_units=500.0,   # 10x larger -> 10x faster
               w_units=100.0,
               n_shells=50,                 # minimal shells
               r_max=50.0,
               r_min=0.05,                  # avoid very dense center
               flag_hydrostatic_initial=True,
               flag_timestep_use_relaxation=True,
               flag_timestep_use_energy=False,   # only relaxation
               )

print(f"Scale time (Gyr): {haloevo.scale_t.to('Gyr').value:.6f}")
print(f"t_relax (dimensionless): {haloevo.t_relax:.6f}")
print(f"t_relax (Gyr): {(haloevo.t_relax * haloevo.scale_t).to('Gyr').value:.6f}")
print(f"sigma_m (dimensionless): {haloevo.sigma_m:.4e}")

# Check time step
haloevo.t_epsilon = 0.1  # 10% per step (very aggressive but OK for testing)
delta_t2 = min(1. / (haloevo.rho * haloevo.v))
print(f"delta_t2 (relaxation): {delta_t2:.6e}")
print(f"delta_t (with t_epsilon=0.1): {delta_t2 * 0.1:.6e}")
print(f"Steps to t=100: {100.0 / (delta_t2 * 0.1):.2e}")

# Run
t0 = time.time()
haloevo.evolve_halo(t_end=100.0,
                    rho_factor_end=1000.0,
                    save_frequency_timing=0.1,
                    t_epsilon=0.1,
                    r_epsilon=1e-10)
elapsed = time.time() - t0

print(f"\nEvolution completed in {elapsed:.1f}s")
print(f"Final t = {haloevo.t:.4f} (dimensionless)")
print(f"Final t = {(haloevo.t * haloevo.scale_t).to('Gyr').value:.4f} Gyr")
print(f"Number of conduction steps: {haloevo.n_conduction}")

list_files, list_times = halorec.glob_pickle_files()
print(f"Saved {len(list_files)} snapshots")

if len(list_files) > 0:
    times_gyr = []
    rho_center = []
    v_center = []
    for f in list_files:
        data = halorec.get_halo_state_pickled(file_halo=f)
        times_gyr.append(data['t'] * haloevo.scale_t.to('Gyr').value)
        rho_center.append(data['rho'][3] * haloevo.scale_rho.to('Msun/pc**3').value)
        v_center.append(np.sqrt(data['p'][3]/data['rho'][3]) * haloevo.scale_v.to('km/s').value)

    times_gyr = np.array(times_gyr)
    rho_center = np.array(rho_center)
    v_center = np.array(v_center)

    print(f"\nCentral density evolution:")
    print(f"  t=0: rho = {rho_center[0]:.4e}, v = {v_center[0]:.2f}")
    idx_min = np.argmin(rho_center)
    print(f"  min rho = {rho_center[idx_min]:.4e} at t = {times_gyr[idx_min]:.4f} Gyr")
    print(f"  final: rho = {rho_center[-1]:.4e} at t = {times_gyr[-1]:.4f} Gyr, v = {v_center[-1]:.2f}")

    # Check if collapse happened
    if rho_center[-1] > rho_center[idx_min] * 2:
        print("  -> COLLAPSE DETECTED!")
        # Find collapse time (when rho crosses back above initial)
        idx_coll = np.searchsorted(rho_center[idx_min:], rho_center[0])
        if idx_coll < len(times_gyr):
            print(f"  -> t_coll ~ {times_gyr[idx_min + idx_coll]:.4f} Gyr")
