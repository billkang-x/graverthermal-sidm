#!/usr/bin/env python3
"""Recompute fluid_ratio_init from scratch, matching compute_fluid_predictions.py exactly."""
import os, sys
import numpy as np

sys.path.insert(0, 'D:/graverthermal-sidm/src/cross_sections')
sys.path.insert(0, 'D:/graverthermal-sidm/src/fluid_runner')
sys.path.insert(0, 'D:/graverthermal-sidm/external/gravothermalsidm')
sys.path.insert(0, 'D:/graverthermal-sidm/src/P3_rescaling')

from astropy import units as ut
from astropy import constants as ct
from dissipative_halo import DissipativeHalo
from SourcePy.record import HaloRecord
from rescale import projected_enclosed_mass

# Build the same evo as P1_elastic_control
out_dir = 'D:/graverthermal-sidm/data/P5_nbody_verify/_tmp_fresh'
os.makedirs(out_dir, exist_ok=True)
# Clean
for f in os.listdir(out_dir):
    if f.endswith('.pickle') or f.endswith('.pkl'):
        os.remove(os.path.join(out_dir, f))

rec = HaloRecord(out_dir)
sigma_m = 0.1
sig_fn = lambda T: np.full(np.atleast_1d(T).shape, sigma_m, dtype=float)
rd_fn = lambda T: np.full(np.atleast_1d(T).shape, 1.0, dtype=float)
evo = DissipativeHalo(rec,
                       sigma_m_eff_callable=sig_fn,
                       rdiss_eff_callable=rd_fn,
                       flag_dissipation=False,
                       profile='NFW', r_s=0.085, rho_s=10.0,
                       sigma_m_with_units=sigma_m, w_units=100.0,
                       n_shells=100, r_max=50.0, r_min=0.02,
                       flag_hydrostatic_initial=True,
                       flag_timestep_use_relaxation=True,
                       flag_timestep_use_energy=True)
evo.t_epsilon = 1e-2
evo.r_epsilon = 1e-12

print(f'scale_r_kpc = {evo.scale_r.to("kpc").value}')
print(f'scale_rho = {evo.scale_rho.to("Msun/pc^3").value}')

# Initial state (matching compute_fluid_predictions.py lines 152-156)
if evo.t == 0:
    if evo.flag_hydrostatic_initial:
        evo.hydrostatic_adjustment()
    evo.save_halo()

# Compute initial ratio (matching compute_fluid_predictions.py lines 158-168)
list_files, _ = rec.glob_pickle_files()
print(f'Pickle files: {list_files}')
data_init = rec.get_halo_state_pickled(file_halo=list_files[-1])
scale_r_kpc = evo.scale_r.to('kpc').value
scale_rho = evo.scale_rho.to('Msun/pc**3').value
r_kpc_init = data_init['r'] * scale_r_kpc
rho_arr_init = data_init['rho'] * scale_rho

print(f'\nr_kpc_init range: {r_kpc_init.min():.6f} - {r_kpc_init.max():.6f}')
print(f'rho range: {rho_arr_init.min():.4e} - {rho_arr_init.max():.4e}')

M_in_init = projected_enclosed_mass(r_kpc_init, rho_arr_init, 0.02,
                                     r_unit='kpc', rho_unit='Msun_pc3')
M_out_init = projected_enclosed_mass(r_kpc_init, rho_arr_init, 0.09,
                                      r_unit='kpc', rho_unit='Msun_pc3')
ratio_init = M_in_init / M_out_init if M_out_init > 0 else np.nan
print(f'\nInitial ratio: {ratio_init:.6f}')
print(f'  M_in(<20pc) = {M_in_init:.6e}')
print(f'  M_out(<90pc) = {M_out_init:.6e}')
print(f'\nCSV fluid_ratio_init = 0.319457')
print(f'CSV fluid_ratio_final = 0.274295')

# Also save the initial state for comparison
import pickle
with open(os.path.join(out_dir, 'time0_fresh.pickle'), 'wb') as f:
    pickle.dump(data_init, f)
print(f'\nSaved fresh initial state to {out_dir}/time0_fresh.pickle')
