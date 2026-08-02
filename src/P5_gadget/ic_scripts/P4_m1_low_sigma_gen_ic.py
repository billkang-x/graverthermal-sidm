#!/usr/bin/env python3
'''Generate NFW initial conditions for P4_m1_low_sigma.

Halo: rho(r) = rho_0 / ((r/r_s) * (1 + r/r_s)^2)
  r_s = 0.085 kpc
  rho_0 = 10.0 Msun/pc^3
  N_particles = 100000
  Box: truncated at r_max = 0.34 kpc

Output: ic.hdf5 in Gadget4 HDF5 snapshot format.
'''
import numpy as np
import h5py
import os

# Constants
G = 4.302e-3   # (pc Msun^-1 (km/s)^2) * (kpc/1000pc) -> kpc Msun^-1 (km/s)^2
G_KPC = 4.302e-3  # kpc Msun^-1 (km/s)^2
MSUN = 1.989e33   # g
KPC = 3.0857e21   # cm
PC = 3.0857e18    # cm
KM = 1e5          # cm

R_S_KPC = 0.085
RHO_0_MSUN_PC3 = 10.0
RHO_0_MSUN_KPC3 = RHO_0_MSUN_PC3 * 1e9   # convert pc^-3 -> kpc^-3
N = 100000
R_MAX_KPC = 0.34

# Mass inside r_max: M(<r) = 4*pi*rho_0*r_s^3 * (ln(1+r/r_s) - r/(r+r_s))
def enclosed_mass(r):
    x = r / R_S_KPC
    return 4 * np.pi * RHO_0_MSUN_KPC3 * R_S_KPC**3 * (np.log1p(x) - x / (1 + x))

M_TOTAL = enclosed_mass(R_MAX_KPC)
M_PARTICLE = M_TOTAL / N
print(f"Halo: r_s={R_S_KPC} kpc, rho_0={RHO_0_MSUN_PC3} Msun/pc^3")
print(f"M_total(<{R_MAX_KPC} kpc) = {M_TOTAL:.4e} Msun, N={N}, m_part={M_PARTICLE:.4e} Msun")

# Sample radii from the NFW profile (inverse-CDF)
# rho(r) * 4*pi*r^2 dr -> dM/dr = 4*pi*rho_0*r_s^3 * r/(r+r_s)^3
# Sample by rejection or by tabulated inverse-CDF
# We'll tabulate CDF(r) and invert numerically
r_grid = np.linspace(1e-6, R_MAX_KPC, 10001)
x_grid = r_grid / R_S_KPC
M_grid = 4 * np.pi * RHO_0_MSUN_KPC3 * R_S_KPC**3 * (np.log1p(x_grid) - x_grid/(1+x_grid))
# CDF normalized to [0,1]
CDF = M_grid / M_grid[-1]

# Sample uniform in [0,1] and invert
u = np.random.uniform(0, 1, N)
r = np.interp(u, CDF, r_grid)

# Sample isotropic direction
phi = np.random.uniform(0, 2*np.pi, N)
cos_theta = np.random.uniform(-1, 1, N)
sin_theta = np.sqrt(1 - cos_theta**2)

x = r * sin_theta * np.cos(phi)
y = r * sin_theta * np.sin(phi)
z = r * cos_theta

# Velocities: circular velocity v_c(r) = sqrt(G*M(r)/r), and isotropic distribution
# for an NFW halo the velocity dispersion profile is more complex (solved by
# solving the Jeans equation). For simplicity use v_circular = sqrt(GM/r) as
# initial tangential velocity, isotropic component ratio sigma_r:sigma_t = 1:2.
M_r = 4 * np.pi * RHO_0_MSUN_KPC3 * R_S_KPC**3 * (np.log1p(r/R_S_KPC) - (r/R_S_KPC)/(1+r/R_S_KPC))
v_circ = np.sqrt(G_KPC * M_r / r)

# Use 3D Maxwellian with sigma = v_circ/sqrt(3) per axis as approximation
# (initial velocities are not critical for SIDM halos - they equilibrate)
sigma_v = v_circ / np.sqrt(3)
vx = np.random.normal(0, sigma_v, N)
vy = np.random.normal(0, sigma_v, N)
vz = np.random.normal(0, sigma_v, N)

print(f"v_circ range: {v_circ.min():.2f} - {v_circ.max():.2f} km/s")

# Write to HDF5 in Gadget4 format
out = "ic.hdf5"
with h5py.File(out, 'w') as f:
    # Header group
    h = f.create_group('Header')
    h.attrs['NumPart_Total'] = np.array([N, 0, 0, 0, 0, 0], dtype=np.uint32)
    h.attrs['NumPart_Total_HighWord'] = np.array([0, 0, 0, 0, 0, 0], dtype=np.uint32)
    h.attrs['NumPart_ThisFile'] = np.array([N, 0, 0, 0, 0, 0], dtype=np.uint32)
    h.attrs['MassTable'] = np.array([M_PARTICLE, 0, 0, 0, 0, 0], dtype=np.float64)
    h.attrs['Time'] = 0.0
    h.attrs['NumFilesPerSnapshot'] = 1
    h.attrs['BoxSize'] = R_MAX_KPC * 2
    h.attrs['Omega0'] = 1.0
    h.attrs['OmegaBaryon'] = 0.0
    h.attrs['OmegaLambda'] = 0.0
    h.attrs['HubbleParam'] = 0.7
    h.attrs['Redshift'] = 0.0
    h.attrs['Flag_Sfr'] = 0
    h.attrs['Flag_Cooling'] = 0
    h.attrs['Flag_StellarAge'] = 0
    h.attrs['Flag_Metals'] = 0
    h.attrs['Flag_Feedback'] = 0
    h.attrs['Flag_DoublePrecision'] = 1

    # Particle group (PartType0 = gas, PartType1 = halo DM)
    # For SIDM, particles go in PartType1 (DM)
    p1 = f.create_group('PartType1')
    p1.create_dataset('ParticleIDs', data=np.arange(1, N+1, dtype=np.uint32))
    p1.create_dataset('Coordinates', data=np.column_stack([x, y, z]).astype(np.float64))
    p1.create_dataset('Velocities', data=np.column_stack([vx, vy, vz]).astype(np.float64))
    p1.create_dataset('Masses', data=np.full(N, M_PARTICLE, dtype=np.float64))
    # Internal energy for DM not needed
    p1.create_dataset('InternalEnergy', data=np.zeros(N, dtype=np.float64))

print(f"Wrote {out}: {os.path.getsize(out)/1e6:.1f} MB")
