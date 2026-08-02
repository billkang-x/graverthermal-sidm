"""Diagnose the time step issue in the elastic halo evolution."""
import sys, os
sys.path.insert(0, r'D:/graverthermal-sidm/external/gravothermalsidm')

import numpy as np
from astropy import units as ut
from astropy import constants as ct
from SourcePy.evolve import Halo
from SourcePy.record import HaloRecord

dir_data = r'D:/graverthermal-sidm/data/P2_runs/_diag'
import shutil
if os.path.exists(dir_data):
    shutil.rmtree(dir_data)

halorec = HaloRecord(dir_data)
haloevo = Halo(halorec,
               profile='NFW',
               r_s=3.6,
               rho_s=7.09e-3,  # Msun/pc^3
               sigma_m_with_units=50.0,  # cm^2/g
               w_units=100.0,  # km/s
               n_shells=200,
               r_max=50.0,
               flag_hydrostatic_initial=True)

print(f"Scale time (Gyr): {haloevo.scale_t.to('Gyr').value:.6f}")
print(f"t_relax (dimensionless): {haloevo.t_relax:.6f}")
print(f"t_relax (Gyr): {(haloevo.t_relax * haloevo.scale_t).to('Gyr').value:.6f}")
print()

# Look at the derived parameters to understand the time step
print("Derived parameters after initialization:")
print(f"  v (1D vel dispersion) range: {haloevo.v.min():.4e} to {haloevo.v.max():.4e}")
print(f"  v in km/s: {haloevo.v.min()*haloevo.scale_v.to('km/s').value:.2f} to {haloevo.v.max()*haloevo.scale_v.to('km/s').value:.2f}")
print(f"  rho range: {haloevo.rho.min():.4e} to {haloevo.rho.max():.4e}")
print(f"  rho (Msun/pc^3): {haloevo.rho.min()*haloevo.scale_rho.to('Msun/pc**3').value:.4e} to {haloevo.rho.max()*haloevo.scale_rho.to('Msun/pc**3').value:.4e}")
print(f"  u range: {haloevo.u.min():.4e} to {haloevo.u.max():.4e}")
print(f"  L range: {haloevo.L.min():.4e} to {haloevo.L.max():.4e}")
print(f"  Kn range: {haloevo.Kn.min():.4e} to {haloevo.Kn.max():.4e}")
print(f"  sigma_m (dimensionless): {haloevo.sigma_m:.4e}")
print()

# Time step components
haloevo.t_epsilon = 1e-4
delta_t1 = min(abs(haloevo.u[0]/(haloevo.L[0]/haloevo.m[0])),
               min(abs(haloevo.u[1:] / ((haloevo.L[1:]-haloevo.L[:-1])/(haloevo.m[1:]-haloevo.m[:-1])))))
delta_t2 = min(1. / (haloevo.rho * haloevo.v))
print(f"  delta_t1 (energy): {delta_t1:.6e}")
print(f"  delta_t2 (relaxation): {delta_t2:.6e}")
print(f"  delta_t (min * t_epsilon): {min(delta_t1, delta_t2) * haloevo.t_epsilon:.6e}")
print(f"  -> in Gyr: {min(delta_t1, delta_t2) * haloevo.t_epsilon * haloevo.scale_t.to('Gyr').value:.6e}")
print()

# Where is the minimum time step?
idx_min_dt2 = np.argmin(1. / (haloevo.rho * haloevo.v))
print(f"  Min delta_t2 at shell {idx_min_dt2}: r={haloevo.r[idx_min_dt2]:.4e}, rho={haloevo.rho[idx_min_dt2]:.4e}, v={haloevo.v[idx_min_dt2]:.4e}")

# The issue: with sigma_m=50 cm^2/g and these halo parameters, the
# dimensionless sigma_m is 0.27, and the relaxation time is short.
# Let's check what sigma_m gives t_relax ~ 100 (so collapse ~ 10 Gyr)
print("\n--- Checking sigma_m needed for t_coll ~ 10 Gyr ---")
target_t_coll_gyr = 10.0
target_t_coll_dimless = target_t_coll_gyr / haloevo.scale_t.to('Gyr').value
print(f"  Target t_coll dimensionless: {target_t_coll_dimless:.2f}")
# t_relax ~ 2/(3 a C sigma_m F v rho_scale) / t_scale
# For elastic, t_coll ~ few * t_relax
# So we need t_relax ~ 20-50 dimensionless
# t_relax = 1.47 for sigma_m=50 -> need sigma_m ~ 50*1.47/30 ~ 2.4
# OR we just need to run to larger t_end
print(f"  Current t_relax = {haloevo.t_relax:.2f} for sigma_m=50")
print(f"  To get t_coll ~ 10 Gyr, need t_end ~ {target_t_coll_dimless:.0f} dimensionless")
print(f"  But each step is ~{min(delta_t1, delta_t2) * haloevo.t_epsilon:.2e} dimensionless")
print(f"  -> Need ~{target_t_coll_dimless / (min(delta_t1, delta_t2) * haloevo.t_epsilon):.2e} steps")
print()

# Alternative: use larger sigma_m (more strongly self-interacting) to speed up
# OR reduce n_shells for testing
print("--- Testing with larger sigma_m=500 and fewer shells ---")
dir_data2 = r'D:/graverthermal-sidm/data/P2_runs/_diag2'
if os.path.exists(dir_data2):
    shutil.rmtree(dir_data2)
halorec2 = HaloRecord(dir_data2)
haloevo2 = Halo(halorec2,
                profile='NFW',
                r_s=3.6,
                rho_s=7.09e-3,
                sigma_m_with_units=500.0,  # 10x larger
                w_units=100.0,
                n_shells=100,  # fewer shells
                r_max=50.0,
                flag_hydrostatic_initial=True)
print(f"  t_relax (dimensionless): {haloevo2.t_relax:.4f}")
print(f"  t_relax (Gyr): {(haloevo2.t_relax * haloevo2.scale_t).to('Gyr').value:.4f}")

haloevo2.t_epsilon = 1e-4
delta_t1_2 = min(abs(haloevo2.u[0]/(haloevo2.L[0]/haloevo2.m[0])),
                 min(abs(haloevo2.u[1:] / ((haloevo2.L[1:]-haloevo2.L[:-1])/(haloevo2.m[1:]-haloevo2.m[:-1])))))
delta_t2_2 = min(1. / (haloevo2.rho * haloevo2.v))
print(f"  delta_t1 (energy): {delta_t1_2:.6e}")
print(f"  delta_t2 (relaxation): {delta_t2_2:.6e}")
print(f"  delta_t: {min(delta_t1_2, delta_t2_2) * haloevo2.t_epsilon:.6e}")
print(f"  Steps to t=10: {10.0 / (min(delta_t1_2, delta_t2_2) * haloevo2.t_epsilon):.2e}")
