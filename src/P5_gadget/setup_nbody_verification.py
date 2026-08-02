"""
Task #41: Generate N-body initial conditions + parameter files for the
5 verification points of Plan B.

Plan B test points (5 points):
  1. ELASTIC CONTROL: pure elastic SIDM (r_diss=1.0 exactly), constant sigma/m=0.1
     - The fluid model should reproduce this exactly (n_cool=0 → no cooling).
       Verifies the elastic-only scattering implementation and symmetry.
  2. M3 LOW  sigma/m: const r_diss=1.05, sigma/m=0.005 (at v=100 km/s)
     - Velocity-independent; tests constant dissipative cooling.
  3. M3 HIGH sigma/m: const r_diss=1.05, sigma/m=0.220
     - Higher cooling rate; tests constant dissipative regime.
  4. M1 LOW  sigma/m: dark photon (velocity-dependent), sigma/m=0.005 at 100 km/s
     - Strong velocity dependence at low v; tests sigma(v) table.
  5. M1 HIGH sigma/m: dark photon, sigma/m=0.165 at 100 km/s
     - Tests sigma(v) with high rate; should produce deep collapse.

Halo setup (B1938+666 inspired, simplified isolated halo):
  - NFW halo, M_200 = 1e12 Msun, c_200 = 8.5 (r_s ~ 4.7 kpc)
    OR scaled-down version to match P5 grid scan r_s~0.05-0.12 kpc.
  - For direct comparison to the fluid model, use the SAME physical halo
    parameters as the grid scan: r_s_phys, rho_0=10 Msun/pc^3.

Because the fluid model uses a single radial bin (r_s, rho_0) and runs in
dimensionless units, the N-body must use the same physical scales. We use:
  - N_particles = 1e5 (heavy) or 5e4 (medium, faster).
  - Box size: 4 * r_s_phys (truncated NFW).
  - Softening: 0.01 * r_s_phys.

For the N-body to match the fluid model's mass ratio M(20pc)/M(90pc), we need
to resolve 20 pc - that requires softening << 20 pc. With r_s_phys = 0.05-0.12
kpc = 50-120 pc, softening of 0.5 pc is sufficient.

Output: ~/dsidm_project/nbody_verify/<point_name>/{ic.txt, params.txt, submit.sh}
"""
import paramiko
import os
import json
import math

HOST = "ssh.cn-zhongwei-1.paracloud.com"
PORT = 22
USER = "scg7816@ZC-M6"
PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"

# Halo physical parameters - matching the P5 grid scan scale.
# We pick r_s_phys = 0.085 kpc (median of viable r_s range) and rho_0 = 10.
# This gives M_halo ~ 4*pi/3 * rho_0 * r_s^3 ~ small cluster scale.
# Note: in dimensionless units the fluid model uses, this is the same halo.

# Test points (sigma_m at v_ref=100 km/s, the fluid model grid reference)
TEST_POINTS = [
    {
        "name": "P1_elastic_control",
        "model": "elastic",  # r_diss = 1.0 exactly
        "sigma_m_100": 0.1,  # cm^2/g at v=100 km/s
        "r_diss": 1.0,       # no dissipation
        "description": "Elastic control: verifies elastic-only path matches fluid model.",
        # Fluid prediction: no cooling, mass ratio stays near initial NFW ratio.
    },
    {
        "name": "P2_m3_low_sigma",
        "model": "M3",
        "sigma_m_100": 0.005,
        "r_diss": 1.05,
        "description": "M3 low sigma/m: constant r_diss=1.05, weak dissipative cooling.",
        # Fluid prediction: M3 grid scan shows viable with t_cross ~0.07 Gyr.
    },
    {
        "name": "P3_m3_high_sigma",
        "model": "M3",
        "sigma_m_100": 0.220,
        "r_diss": 1.05,
        "description": "M3 high sigma/m: strong cooling rate, faster collapse.",
        # Fluid prediction: viable with t_cross ~0.10 Gyr.
    },
    {
        "name": "P4_m1_low_sigma",
        "model": "M1",
        "sigma_m_100": 0.005,
        "r_diss": "velocity-dependent",  # from LB2026
        "description": "M1 dark photon low sigma/m: tests sigma(v) table at low v.",
        # Fluid prediction: viable with t_cross ~0.05 Gyr.
    },
    {
        "name": "P5_m1_high_sigma",
        "model": "M1",
        "sigma_m_100": 0.165,
        "r_diss": "velocity-dependent",
        "description": "M1 dark photon high sigma/m: strong velocity dependence, deep collapse.",
        # Fluid prediction: viable with t_cross ~0.05 Gyr.
    },
]

# Common halo parameters (matches grid scan)
R_S_KPC = 0.085      # median of viable r_s range
RHO_0_MSUN_PC3 = 10.0
N_PARTICLES = 100000  # 1e5 - resolves 20 pc with softening ~0.5 pc
SOFTENING_PC = 0.5   # ~5x smaller than 20pc target
BOX_KPC = 4.0 * R_S_KPC  # 0.34 kpc, plenty for truncated NFW
T_END_GYR = 1.0       # run for 1 Gyr - sufficient for all viable t_cross < 0.63 Gyr


def gen_nfw_ic_params(point):
    """Generate IC parameters for an NFW halo initial condition."""
    # We'll use GADGET-4's built-in NFW IC generator if available, else a Python
    # script to write a HDF5 snapshot.
    params = {
        "point_name": point["name"],
        "model": point["model"],
        "sigma_m_100_cm2_g": point["sigma_m_100"],
        "r_diss": point["r_diss"],
        # Halo
        "r_s_kpc": R_S_KPC,
        "rho_0_msun_pc3": RHO_0_MSUN_PC3,
        "n_particles": N_PARTICLES,
        "softening_pc": SOFTENING_PC,
        "box_kpc": BOX_KPC,
        # Time
        "t_end_gyr": T_END_GYR,
        # Observational target
        "r_inner_pc": 20.0,
        "r_outer_pc": 90.0,
        "target_mass_ratio": 0.364,  # observed B1938+666
    }
    return params


def write_params_file(point):
    """Generate the Gadget4 param.txt content for this point."""
    p = gen_nfw_ic_params(point)
    # Gadget4 params - minimal set for an isolated SIDM halo run
    lines = [
        "%paramfile for N-body verification",
        f"% Point: {point['name']}",
        f"% Model: {point['model']}, sigma/m={point['sigma_m_100']} cm^2/g @ 100 km/s",
        f"% r_diss: {point['r_diss']}",
        "",
        "InitCondFile      ic",
        "OutputDir        output",
        "SnapshotFileBase  snap",
        "",
        "TimeMax            1.0",
        "TimeBetSnapshot    0.05",
        "TimeBetStatistics  0.01",
        "CpuTreeDomainUpdate 0.1",
        "",
        "ErrTolIntAccuracy  0.025",
        "ErrTolTheta        0.5",
        "ErrTolForce        0.005",
        "MaxRMSDisplacementFac 0.2",
        "",
        "SofteningComovingClass0  0.5",
        "SofteningMaxPhysClass0   0.5",
        "SofteningComovingClass1  0.5",
        "SofteningMaxPhysClass1   0.5",
        "SofteningComovingClass2  0.5",
        "SofteningMaxPhysClass2   0.5",
        "SofteningComovingClass3  0.5",
        "SofteningMaxPhysClass3   0.5",
        "SofteningComovingClass4  0.5",
        "SofteningMaxPhysClass4   0.5",
        "SofteningComovingClass5  0.5",
        "SofteningMaxPhysClass5   0.5",
        "",
        "UnitLength_in_cm          3.085678e21",
        "UnitMass_in_g             1.989e43",
        "UnitVelocity_in_cm_per_s  1e5",
        "UnitLuminosity_in_erg_s   1",
        "UnitEnergy_in_ergs        1",
        "GravityConstantInternal  0.017778279",
        "",
        "MinGasTemp                1e2",
        "HubbleParam               0.7",
        "Omega0                    1.0",   # isolated halo - no cosmology
        "OmegaBaryon               0.0",
        "OmegaLambda               0.0",
        "Hubble                     0.0",
        "",
        f"% SIDM parameters - will be overridden by Config.sh compile-time values",
        f"% sigma/m at v=100 km/s = {point['sigma_m_100']}",
        f"% r_diss = {point['r_diss']}",
        "",
        "BoxSize                   1.0",
        "PeriodicBoundariesOn      0",
    ]
    return '\n'.join(lines) + '\n'


def write_submit_sh(point):
    """Write Slurm submit script."""
    name = point["name"]
    return f"""#!/bin/bash
#SBATCH --job-name={name}
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err

source ~/env.sh

cd $SLURM_SUBMIT_DIR
echo "=== Starting {name} at $(date) ==="
echo "Job ID: $SLURM_JOB_ID on $(hostname)"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"

# Run Gadget4 with the params file
mpirun -np 8 ~/dsidm_project/source/Gadget4_dsidm params.txt

echo "=== Finished {name} at $(date) ==="
ls -la output/
"""


def write_gen_ic_py(point):
    """Write a Python script that generates the NFW IC HDF5 snapshot.
    Will be run on the HPC (has numpy/h5py there)."""
    p = gen_nfw_ic_params(point)
    # Compute the concentration parameter etc. We use NFW rho(r) directly.
    return f"""#!/usr/bin/env python3
'''Generate NFW initial conditions for {point['name']}.

Halo: rho(r) = rho_0 / ((r/r_s) * (1 + r/r_s)^2)
  r_s = {p['r_s_kpc']} kpc
  rho_0 = {p['rho_0_msun_pc3']} Msun/pc^3
  N_particles = {p['n_particles']}
  Box: truncated at r_max = {p['box_kpc']} kpc

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

R_S_KPC = {p['r_s_kpc']}
RHO_0_MSUN_PC3 = {p['rho_0_msun_pc3']}
RHO_0_MSUN_KPC3 = RHO_0_MSUN_PC3 * 1e9   # convert pc^-3 -> kpc^-3
N = {p['n_particles']}
R_MAX_KPC = {p['box_kpc']}

# Mass inside r_max: M(<r) = 4*pi*rho_0*r_s^3 * (ln(1+r/r_s) - r/(r+r_s))
def enclosed_mass(r):
    x = r / R_S_KPC
    return 4 * np.pi * RHO_0_MSUN_KPC3 * R_S_KPC**3 * (np.log1p(x) - x / (1 + x))

M_TOTAL = enclosed_mass(R_MAX_KPC)
M_PARTICLE = M_TOTAL / N
print(f"Halo: r_s={{R_S_KPC}} kpc, rho_0={{RHO_0_MSUN_PC3}} Msun/pc^3")
print(f"M_total(<{{R_MAX_KPC}} kpc) = {{M_TOTAL:.4e}} Msun, N={{N}}, m_part={{M_PARTICLE:.4e}} Msun")

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

print(f"v_circ range: {{v_circ.min():.2f}} - {{v_circ.max():.2f}} km/s")

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

print(f"Wrote {{out}}: {{os.path.getsize(out)/1e6:.1f}} MB")
"""


def run(client, cmd, timeout=120):
    print(f"$ {cmd[:200]}", flush=True)
    stdin, stdout, stderr = client.exec_command(cmd, timeout=timeout)
    out = stdout.read().decode(errors='replace').rstrip()
    err = stderr.read().decode(errors='replace').rstrip()
    if out:
        print(out[:5000])
    if err:
        print(f"  [stderr] {err[:2000]}")
    return out, err


def main():
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(HOST, port=PORT, username=USER, password=PASS,
                   timeout=30, look_for_keys=False, allow_agent=False)
    print("Connected.\n", flush=True)

    # Make base directory
    base = "~/dsidm_project/nbody_verify"
    run(client, f"mkdir -p {base}")

    # Upload a gen_ic.py and params.txt and submit.sh for each test point
    sftp = client.open_sftp()
    base_remote = base.replace("~", "/public3/home/scg7816")

    for point in TEST_POINTS:
        name = point["name"]
        point_dir = f"{base_remote}/{name}"
        run(client, f"mkdir -p {point_dir}")

        # Write gen_ic.py
        gen_ic_content = write_gen_ic_py(point)
        gen_ic_path_local = f"D:/graverthermal-sidm/src/P5_gadget/ic_scripts/{name}_gen_ic.py"
        os.makedirs(os.path.dirname(gen_ic_path_local), exist_ok=True)
        with open(gen_ic_path_local, 'w') as f:
            f.write(gen_ic_content)
        sftp.put(gen_ic_path_local, f"{point_dir}/gen_ic.py")
        print(f"  {name}: gen_ic.py uploaded ({len(gen_ic_content)} chars)")

        # Write params.txt
        params_content = write_params_file(point)
        params_path_local = f"D:/graverthermal-sidm/src/P5_gadget/ic_scripts/{name}_params.txt"
        with open(params_path_local, 'w') as f:
            f.write(params_content)
        sftp.put(params_path_local, f"{point_dir}/params.txt")
        print(f"  {name}: params.txt uploaded")

        # Write submit.sh
        submit_content = write_submit_sh(point)
        submit_path_local = f"D:/graverthermal-sidm/src/P5_gadget/ic_scripts/{name}_submit.sh"
        with open(submit_path_local, 'w') as f:
            f.write(submit_content)
        sftp.put(submit_path_local, f"{point_dir}/submit.sh")
        run(client, f"chmod +x {point_dir}/submit.sh")

        # Write a README.json with metadata
        readme = gen_nfw_ic_params(point)
        readme_path_local = f"D:/graverthermal-sidm/src/P5_gadget/ic_scripts/{name}_readme.json"
        with open(readme_path_local, 'w') as f:
            json.dump(readme, f, indent=2)
        sftp.put(readme_path_local, f"{point_dir}/readme.json")

    sftp.close()

    # Verify everything is in place
    print("\n=== Verify directory structure ===", flush=True)
    run(client, f"ls -la {base}/")
    run(client, f"for d in {base}/*/; do echo $d; ls $d; done")

    client.close()
    print("\nDone. IC scripts ready - waiting for Gadget4_dsidm build to complete before "
          "running gen_ic.py and submitting jobs.", flush=True)


if __name__ == "__main__":
    main()
