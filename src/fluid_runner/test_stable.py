"""
Evolve elastic halo with controlled parameters.

Use t_epsilon=1e-2 (proven stable in earlier test that reached t=9),
n_shells=100, and both time step criteria.
"""
import sys, os, time
sys.path.insert(0, r'D:/graverthermal-sidm/external/gravothermalsidm')

import numpy as np
from astropy import units as ut
from SourcePy.evolve import Halo
from SourcePy.record import HaloRecord


def evolve_with_checks(haloevo, t_end=200.0, rho_factor_end=1000.0,
                        max_steps=200000, save_every_steps=200,
                        t_epsilon=1e-2, r_epsilon=1e-12, verbose=True):
    """Evolve with NaN checking and step-based saving."""
    haloevo.t_epsilon = t_epsilon
    haloevo.r_epsilon = r_epsilon

    if haloevo.t == 0:
        if haloevo.flag_hydrostatic_initial:
            haloevo.hydrostatic_adjustment()
        haloevo.save_halo()

    t0 = time.time()
    last_report = time.time()
    last_save_t = haloevo.t

    for step in range(max_steps):
        # Conduct heat
        haloevo.conduct_heat()

        # Check for NaN
        if np.any(np.isnan(haloevo.u)) or np.any(np.isnan(haloevo.p)):
            print(f"!!! NaN detected after conduct_heat at step {step}, t={haloevo.t:.6e}")
            break

        # Hydrostatic adjustment
        haloevo.hydrostatic_adjustment()

        # Check for NaN
        if np.any(np.isnan(haloevo.rho)) or np.any(np.isnan(haloevo.r)):
            print(f"!!! NaN detected after hydrostatic at step {step}, t={haloevo.t:.6e}")
            break

        # Check termination
        current_rho = haloevo.get_central_quantity(haloevo.rho)
        if haloevo.t >= t_end:
            haloevo.save_halo()
            print(f"Done: reached t={haloevo.t:.4f} (target {t_end})")
            break
        elif current_rho > rho_factor_end * haloevo.rho_center:
            haloevo.save_halo()
            print(f"Done: rho_center={current_rho:.4e} > {rho_factor_end}*rho0 (collapse)")
            break

        # Save periodically by step count OR by time fraction
        if (step + 1) % save_every_steps == 0:
            haloevo.save_halo()

        # Progress report
        if verbose and (time.time() - last_report > 5.0):
            elapsed = time.time() - t0
            rate = (step + 1) / elapsed
            rho_now = current_rho * haloevo.scale_rho.to('Msun/pc**3').value
            t_gyr = haloevo.t * haloevo.scale_t.to('Gyr').value
            v_now = haloevo.v[3] * haloevo.scale_v.to('km/s').value
            print(f"  step {step+1}: t={haloevo.t:.2f} ({t_gyr:.3f} Gyr), "
                  f"rho_c={rho_now:.3e}, v_c={v_now:.1f} km/s, {rate:.0f} steps/s")
            last_report = time.time()

    elapsed = time.time() - t0
    print(f"\nEvolution: {haloevo.n_conduction} steps in {elapsed:.1f}s")
    print(f"Final t = {haloevo.t:.4f} ({(haloevo.t * haloevo.scale_t).to('Gyr').value:.4f} Gyr)")
    return haloevo


# Run with proven-stable parameters
dir_data = r'D:/graverthermal-sidm/data/P2_runs/_test_stable'
import shutil
if os.path.exists(dir_data):
    shutil.rmtree(dir_data)

halorec = HaloRecord(dir_data)
haloevo = Halo(halorec,
               profile='NFW',
               r_s=3.6,
               rho_s=7.09e-3,
               sigma_m_with_units=50.0,  # Schmidt benchmark
               w_units=100.0,
               n_shells=100,
               r_max=50.0,
               r_min=0.02,
               flag_hydrostatic_initial=True,
               flag_timestep_use_relaxation=True,
               flag_timestep_use_energy=True,
               )

print(f"Scale time (Gyr): {haloevo.scale_t.to('Gyr').value:.6f}")
print(f"t_relax (dimensionless): {haloevo.t_relax:.4f}")
print(f"t_relax (Gyr): {(haloevo.t_relax * haloevo.scale_t).to('Gyr').value:.4f}")

# Use t_epsilon=1e-2 (proven stable) and run to t=200
evolve_with_checks(haloevo, t_end=200.0, rho_factor_end=1000.0,
                   max_steps=500000, save_every_steps=500,
                   t_epsilon=1e-2, r_epsilon=1e-12)

# Analyze
list_files, list_times = halorec.glob_pickle_files()
print(f"\nSnapshots: {len(list_files)}")
if len(list_files) > 0:
    # Filter out the broken file
    valid_files = []
    for f in list_files:
        try:
            data = halorec.get_halo_state_pickled(file_halo=f)
            if data and not np.any(np.isnan(data.get('rho', [np.nan]))):
                valid_files.append(f)
        except:
            pass

    print(f"Valid snapshots: {len(valid_files)}")
    if len(valid_files) > 0:
        times_gyr = []
        rho_center = []
        v_center = []
        for f in valid_files:
            data = halorec.get_halo_state_pickled(file_halo=f)
            times_gyr.append(data['t'] * haloevo.scale_t.to('Gyr').value)
            rho_center.append(data['rho'][3] * haloevo.scale_rho.to('Msun/pc**3').value)
            v_center.append(np.sqrt(data['p'][3]/data['rho'][3]) * haloevo.scale_v.to('km/s').value)
        times_gyr = np.array(times_gyr)
        rho_center = np.array(rho_center)
        v_center = np.array(v_center)

        idx_min = np.argmin(rho_center)
        print(f"\nt_core = {times_gyr[idx_min]:.4f} Gyr, rho_min = {rho_center[idx_min]:.4e}")
        print(f"Final: t = {times_gyr[-1]:.4f} Gyr, rho = {rho_center[-1]:.4e}, v = {v_center[-1]:.2f}")
        print(f"rho_final/rho_initial = {rho_center[-1]/rho_center[0]:.4e}")

        # Save plot
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
        ax1.semilogy(times_gyr, rho_center, 'b-')
        ax1.set_ylabel(r'$\rho_{\rm center}$ [$M_\odot$/pc$^3$]')
        ax1.set_title(r'Elastic: $\sigma_T/m = 50$ cm$^2$/g, $r_{\rm diss}=1.0$')
        ax1.axvline(times_gyr[idx_min], color='r', ls='--', alpha=0.5,
                    label=f't_core={times_gyr[idx_min]:.2f} Gyr')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        ax2.plot(times_gyr, v_center, 'b-')
        ax2.set_xlabel('Time [Gyr]')
        ax2.set_ylabel(r'$\nu_{\rm center}$ [km/s]')
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(r'D:/graverthermal-sidm/figures/test_stable_evolution.png', dpi=100)
        print(f"Plot saved: figures/test_stable_evolution.png")
