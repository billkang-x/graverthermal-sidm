"""Test elastic halo with larger t_epsilon for faster evolution."""
import sys, os, time
sys.path.insert(0, r'D:/graverthermal-sidm/external/gravothermalsidm')

import numpy as np
from astropy import units as ut
from SourcePy.evolve import Halo
from SourcePy.record import HaloRecord

dir_data = r'D:/graverthermal-sidm/data/P2_runs/_test_fast'
import shutil
if os.path.exists(dir_data):
    shutil.rmtree(dir_data)

halorec = HaloRecord(dir_data)
haloevo = Halo(halorec,
               profile='NFW',
               r_s=3.6,
               rho_s=7.09e-3,
               sigma_m_with_units=50.0,
               w_units=100.0,
               n_shells=200,
               r_max=50.0,
               flag_hydrostatic_initial=True)

print(f"Scale time (Gyr): {haloevo.scale_t.to('Gyr').value:.6f}")
print(f"t_relax (dimensionless): {haloevo.t_relax:.4f}")
print(f"t_relax (Gyr): {(haloevo.t_relax * haloevo.scale_t).to('Gyr').value:.4f}")

# Use larger t_epsilon = 1e-2 (1% energy change per step)
# and t_end = 200 (to reach ~10 Gyr)
t0 = time.time()
haloevo.evolve_halo(t_end=200.0,
                    rho_factor_end=1e4,
                    save_frequency_timing=0.1,
                    t_epsilon=1e-2,   # 100x larger
                    r_epsilon=1e-14)
elapsed = time.time() - t0

print(f"\nEvolution completed in {elapsed:.1f}s")
print(f"Final t = {haloevo.t:.4f} (dimensionless)")
print(f"Final t = {(haloevo.t * haloevo.scale_t).to('Gyr').value:.4f} Gyr")
print(f"Number of conduction steps: {haloevo.n_conduction}")

# Extract central density evolution
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
    print(f"  t=0: rho = {rho_center[0]:.4e} Msun/pc^3, v = {v_center[0]:.2f} km/s")
    idx_min = np.argmin(rho_center)
    print(f"  min rho = {rho_center[idx_min]:.4e} at t = {times_gyr[idx_min]:.4f} Gyr")
    print(f"  final rho = {rho_center[-1]:.4e} at t = {times_gyr[-1]:.4f} Gyr")

    # Schmidt et al. 2026 Fig 2 expects for sigma_m=50, r_diss=1.0:
    # t_core ~ 4-5 Gyr, t_coll ~ 8-10 Gyr
    print(f"\n  Schmidt et al. 2026 expects: t_core ~ 4-5 Gyr, t_coll ~ 8-10 Gyr")
