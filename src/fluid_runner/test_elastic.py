"""Quick test: elastic halo evolution to validate against Schmidt et al. 2026."""
import sys, os, time
sys.path.insert(0, r'D:/graverthermal-sidm/external/gravothermalsidm')

import numpy as np
from astropy import units as ut
from astropy import constants as ct
from SourcePy.evolve import Halo
from SourcePy.record import HaloRecord

# Schmidt et al. 2026 setup
RS_KPC = 3.6
RHO0_MSUN_PER_PC3 = 7.09e-3  # = 7.09e6 Msun/kpc^3
SIGMA_M = 50.0  # cm^2/g
W_UNITS = 100.0  # km/s

dir_data = r'D:/graverthermal-sidm/data/P2_runs/_test_elastic'
import shutil
if os.path.exists(dir_data):
    shutil.rmtree(dir_data)

halorec = HaloRecord(dir_data)
haloevo = Halo(halorec,
               profile='NFW',
               r_s=RS_KPC,
               rho_s=RHO0_MSUN_PER_PC3,
               sigma_m_with_units=SIGMA_M,
               w_units=W_UNITS,
               n_shells=200,
               r_max=50.0,
               flag_hydrostatic_initial=True)

print(f"Scale radius: {haloevo.scale_r}")
print(f"Scale density: {haloevo.scale_rho}")
print(f"Scale velocity: {haloevo.scale_v}")
print(f"Scale time (Gyr): {haloevo.scale_t.to('Gyr').value:.4f}")
print(f"sigma_m (dimensionless): {haloevo.sigma_m:.4e}")
print(f"w (dimensionless): {haloevo.w:.4e}")
print(f"t_relax (dimensionless): {haloevo.t_relax:.4e}")
print(f"t_relax (Gyr): {(haloevo.t_relax * haloevo.scale_t).to('Gyr').value:.4f}")

# Run a short evolution
t0 = time.time()
haloevo.evolve_halo(t_end=200.0,  # need ~100-200 for collapse (t_relax=1.47, collapse ~10 Gyr / 0.05 Gyr ~ 200)
                    rho_factor_end=1e4,
                    save_frequency_timing=0.1,
                    t_epsilon=1e-4,
                    r_epsilon=1e-14)
elapsed = time.time() - t0

print(f"\nEvolution completed in {elapsed:.2f}s")
print(f"Final dimensionless time t = {haloevo.t:.4f}")
print(f"Final time (Gyr) = {(haloevo.t * haloevo.scale_t).to('Gyr').value:.4f}")

# Extract central density evolution
list_files, list_times = halorec.glob_pickle_files()
print(f"\nSaved {len(list_files)} snapshots")

if len(list_files) > 0:
    times_gyr = []
    rho_center = []
    for f in list_files:
        data = halorec.get_halo_state_pickled(file_halo=f)
        times_gyr.append(data['t'] * haloevo.scale_t.to('Gyr').value)
        rho_center.append(data['rho'][3] * haloevo.scale_rho.to('Msun/pc**3').value)

    times_gyr = np.array(times_gyr)
    rho_center = np.array(rho_center)

    print(f"\nCentral density evolution:")
    print(f"  t=0: rho = {rho_center[0]:.4e} Msun/pc^3")
    print(f"  min rho = {rho_center.min():.4e} at t = {times_gyr[np.argmin(rho_center)]:.4f} Gyr")
    print(f"  final rho = {rho_center[-1]:.4e} at t = {times_gyr[-1]:.4f} Gyr")

    # Compare with Schmidt et al. 2026 Fig 2:
    # For sigma_T/m=50, r_diss=1.0, t_core ~ 4-5 Gyr, t_coll ~ 8-10 Gyr
    print(f"\n  (Schmidt et al. 2026 expects t_core ~ 4-5 Gyr, t_coll ~ 8-10 Gyr for sigma_m=50)")
