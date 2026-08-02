"""
Generate NFW ICs LOCALLY (where we have numpy/h5py), then upload to HPC.

We generate ONE identical halo IC file (same NFW halo), and reuse it for all 5
test points because the halo itself is the same — only the SIDM parameters
differ, which are read from the params.txt at runtime.

This is much simpler than the per-point gen_ic.py approach.
"""
import os
import numpy as np
import h5py
import paramiko

HOST = "ssh.cn-zhongwei-1.paracloud.com"
PORT = 22
USER = "scg7816@ZC-M6"
PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"

# Halo parameters (matching setup_nbody_verification.py)
R_S_KPC = 0.085
RHO_0_MSUN_PC3 = 10.0
RHO_0_MSUN_KPC3 = RHO_0_MSUN_PC3 * 1e9  # 1e10 Msun/kpc^3
N_PARTICLES = 100000
R_MAX_KPC = 4.0 * R_S_KPC  # 0.34 kpc
SOFTENING_KPC = 0.0005  # 0.5 pc

# Constants (kpc, Msun, km/s)
G_KPC = 4.302e-3  # kpc Msun^-1 (km/s)^2


def enclosed_mass(r):
    x = r / R_S_KPC
    return 4 * np.pi * RHO_0_MSUN_KPC3 * R_S_KPC**3 * (np.log1p(x) - x / (1 + x))


def gen_nfw_ic(out_path):
    """Generate NFW IC HDF5 file."""
    print(f"Generating NFW IC: N={N_PARTICLES}, r_s={R_S_KPC} kpc, rho_0={RHO_0_MSUN_PC3} Msun/pc^3")

    # Compute total mass
    M_TOTAL = enclosed_mass(R_MAX_KPC)
    M_PARTICLE = M_TOTAL / N_PARTICLES
    print(f"  M_total(<{R_MAX_KPC:.3f} kpc) = {M_TOTAL:.4e} Msun")
    print(f"  m_particle = {M_PARTICLE:.4e} Msun")

    # Sample radii via inverse CDF
    r_grid = np.linspace(1e-6, R_MAX_KPC, 10001)
    x_grid = r_grid / R_S_KPC
    M_grid = 4 * np.pi * RHO_0_MSUN_KPC3 * R_S_KPC**3 * (np.log1p(x_grid) - x_grid / (1 + x_grid))
    CDF = M_grid / M_grid[-1]

    np.random.seed(42)  # reproducible
    u = np.random.uniform(0, 1, N_PARTICLES)
    r = np.interp(u, CDF, r_grid)

    # Isotropic direction
    phi = np.random.uniform(0, 2 * np.pi, N_PARTICLES)
    cos_theta = np.random.uniform(-1, 1, N_PARTICLES)
    sin_theta = np.sqrt(1 - cos_theta**2)
    x = r * sin_theta * np.cos(phi)
    y = r * sin_theta * np.sin(phi)
    z = r * cos_theta

    # Velocities: use circular velocity with small isotropic dispersion
    # (initial conditions matter little for SIDM - they virialize quickly)
    M_r = 4 * np.pi * RHO_0_MSUN_KPC3 * R_S_KPC**3 * (np.log1p(r / R_S_KPC) - (r / R_S_KPC) / (1 + r / R_S_KPC))
    v_circ = np.sqrt(G_KPC * M_r / r)
    sigma_v = v_circ / np.sqrt(3)  # rough isotropic estimate
    vx = np.random.normal(0, sigma_v, N_PARTICLES)
    vy = np.random.normal(0, sigma_v, N_PARTICLES)
    vz = np.random.normal(0, sigma_v, N_PARTICLES)

    print(f"  v_circ range: {v_circ.min():.2f} - {v_circ.max():.2f} km/s")

    # Write HDF5 in Gadget4 format
    with h5py.File(out_path, 'w') as f:
        h = f.create_group('Header')
        h.attrs['NumPart_Total'] = np.array([N_PARTICLES, 0, 0, 0, 0, 0], dtype=np.uint32)
        h.attrs['NumPart_Total_HighWord'] = np.array([0, 0, 0, 0, 0, 0], dtype=np.uint32)
        h.attrs['NumPart_ThisFile'] = np.array([N_PARTICLES, 0, 0, 0, 0, 0], dtype=np.uint32)
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

        # Particles in PartType1 (DM)
        p1 = f.create_group('PartType1')
        p1.create_dataset('ParticleIDs', data=np.arange(1, N_PARTICLES + 1, dtype=np.uint32))
        p1.create_dataset('Coordinates',
                          data=np.column_stack([x, y, z]).astype(np.float64))
        p1.create_dataset('Velocities',
                          data=np.column_stack([vx, vy, vz]).astype(np.float64))
        p1.create_dataset('Masses',
                          data=np.full(N_PARTICLES, M_PARTICLE, dtype=np.float64))
        p1.create_dataset('InternalEnergy',
                          data=np.zeros(N_PARTICLES, dtype=np.float64))

    size_mb = os.path.getsize(out_path) / 1e6
    print(f"  Wrote {out_path}: {size_mb:.1f} MB")
    return M_TOTAL, M_PARTICLE


def upload_to_hpc(local_path, remote_path):
    """Upload file to HPC via SFTP."""
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(HOST, port=PORT, username=USER, password=PASS,
                   timeout=30, look_for_keys=False, allow_agent=False)
    sftp = client.open_sftp()
    print(f"  Uploading {local_path} -> {remote_path}...")
    sftp.put(local_path, remote_path)
    sftp.close()
    client.close()
    print(f"  Uploaded.")


def main():
    print("=" * 70)
    print("Generate NFW IC locally and upload to HPC for all 5 verification points")
    print("=" * 70)

    # 1. Generate one shared IC file (same halo for all 5 points - only SIDM params differ)
    local_dir = "D:/graverthermal-sidm/data/P5_nbody_verify/ics"
    os.makedirs(local_dir, exist_ok=True)
    ic_local = os.path.join(local_dir, "ic.hdf5")

    if not os.path.exists(ic_local):
        M_total, m_part = gen_nfw_ic(ic_local)
    else:
        print(f"IC file already exists: {ic_local}")
        # Read mass from existing file
        with h5py.File(ic_local, 'r') as f:
            M_total = f['Header'].attrs['MassTable'][0] * f['Header'].attrs['NumPart_Total'][0]
            m_part = f['Header'].attrs['MassTable'][0]
        print(f"  M_total = {M_total:.4e}, m_particle = {m_part:.4e}")

    # 2. Upload to each point's directory on the HPC
    print("\nUploading IC to 5 HPC point directories...")
    base_remote = "/public3/home/scg7816/dsidm_project/nbody_verify"
    for point_name in ["P1_elastic_control", "P2_m3_low_sigma", "P3_m3_high_sigma",
                        "P4_m1_low_sigma", "P5_m1_high_sigma"]:
        remote_path = f"{base_remote}/{point_name}/ic.hdf5"
        print(f"\n  {point_name}:")
        upload_to_hpc(ic_local, remote_path)

    print("\nDone. All 5 IC files uploaded. Ready to submit Slurm jobs.")


if __name__ == "__main__":
    main()
