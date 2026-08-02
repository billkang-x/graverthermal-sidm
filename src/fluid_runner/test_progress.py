"""
Custom halo evolution with progress reporting and step limit.

The upstream evolve_halo() can get stuck when time steps become tiny.
This wrapper adds:
  - Progress reporting every N steps
  - A maximum step count to prevent infinite loops
  - Optional save by step count (not by time fraction)
"""
import sys, os, time
sys.path.insert(0, r'D:/graverthermal-sidm/external/gravothermalsidm')

import numpy as np
from astropy import units as ut
from SourcePy.evolve import Halo
from SourcePy.record import HaloRecord


def evolve_with_progress(haloevo, t_end=100.0, rho_factor_end=1000.0,
                          max_steps=50000, save_every=500,
                          t_epsilon=0.1, r_epsilon=1e-10, verbose=True):
    """Evolve halo with progress reporting and step limit."""
    haloevo.t_epsilon = t_epsilon
    haloevo.r_epsilon = r_epsilon

    # Save initial state if t==0
    if haloevo.t == 0:
        if haloevo.flag_hydrostatic_initial:
            haloevo.hydrostatic_adjustment()
        haloevo.save_halo()

    t0 = time.time()
    last_report = time.time()

    for step in range(max_steps):
        # Conduct heat
        try:
            haloevo.conduct_heat()
        except Exception as e:
            print(f"!!! Error at step {step}, t={haloevo.t:.6e}: {e}")
            haloevo.save_halo(prefix=haloevo.record.prefix_debug)
            break

        # Hydrostatic adjustment
        try:
            haloevo.hydrostatic_adjustment()
        except Exception as e:
            print(f"!!! Error at step {step}, t={haloevo.t:.6e}: {e}")
            haloevo.save_halo(prefix=haloevo.record.prefix_debug)
            break

        # Check termination
        current_rho = haloevo.get_central_quantity(haloevo.rho)
        if haloevo.t >= t_end:
            haloevo.save_halo()
            if verbose:
                print(f"Done: reached t={haloevo.t:.4f} (target {t_end})")
            break
        elif current_rho > rho_factor_end * haloevo.rho_center:
            haloevo.save_halo()
            if verbose:
                print(f"Done: rho_center={current_rho:.4e} > {rho_factor_end}*rho0")
            break

        # Save periodically
        if (step + 1) % save_every == 0:
            haloevo.save_halo()

        # Progress report
        if verbose and (time.time() - last_report > 5.0):
            elapsed = time.time() - t0
            rate = (step + 1) / elapsed
            rho_now = current_rho * haloevo.scale_rho.to('Msun/pc**3').value
            t_gyr = haloevo.t * haloevo.scale_t.to('Gyr').value
            print(f"  step {step+1}/{max_steps}: t={haloevo.t:.4f} ({t_gyr:.4f} Gyr), "
                  f"rho_c={rho_now:.4e}, rate={rate:.1f} steps/s")
            last_report = time.time()

    elapsed = time.time() - t0
    print(f"\nEvolution: {haloevo.n_conduction} steps in {elapsed:.1f}s "
          f"({haloevo.n_conduction/elapsed:.1f} steps/s)")
    print(f"Final t = {haloevo.t:.4f} ({(haloevo.t * haloevo.scale_t).to('Gyr').value:.4f} Gyr)")

    return haloevo


# Run test
dir_data = r'D:/graverthermal-sidm/data/P2_runs/_test_progress'
import shutil
if os.path.exists(dir_data):
    shutil.rmtree(dir_data)

halorec = HaloRecord(dir_data)
haloevo = Halo(halorec,
               profile='NFW',
               r_s=3.6,
               rho_s=7.09e-3,
               sigma_m_with_units=500.0,
               w_units=100.0,
               n_shells=50,
               r_max=50.0,
               r_min=0.05,
               flag_hydrostatic_initial=True,
               flag_timestep_use_relaxation=True,
               flag_timestep_use_energy=False,
               )

print(f"Scale time (Gyr): {haloevo.scale_t.to('Gyr').value:.6f}")
print(f"t_relax (dimensionless): {haloevo.t_relax:.6f}")

evolve_with_progress(haloevo, t_end=200.0, rho_factor_end=1000.0,
                     max_steps=100000, save_every=1000,
                     t_epsilon=0.1, r_epsilon=1e-10)

# Analyze
list_files, list_times = halorec.glob_pickle_files()
print(f"\nSnapshots: {len(list_files)}")
if len(list_files) > 0:
    times_gyr = []
    rho_center = []
    for f in list_files:
        data = halorec.get_halo_state_pickled(file_halo=f)
        times_gyr.append(data['t'] * haloevo.scale_t.to('Gyr').value)
        rho_center.append(data['rho'][3] * haloevo.scale_rho.to('Msun/pc**3').value)
    times_gyr = np.array(times_gyr)
    rho_center = np.array(rho_center)

    idx_min = np.argmin(rho_center)
    print(f"t_core = {times_gyr[idx_min]:.4f} Gyr, rho_min = {rho_center[idx_min]:.4e}")
    print(f"Final: t = {times_gyr[-1]:.4f} Gyr, rho = {rho_center[-1]:.4e}")
    print(f"rho_final/rho_initial = {rho_center[-1]/rho_center[0]:.4e}")
