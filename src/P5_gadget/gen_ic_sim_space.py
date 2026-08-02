#!/usr/bin/env python3
"""Generate NFW IC in SIMULATION SPACE for rescaling-symmetry verification.

The physical-space halo (r_s=0.085 kpc, rho_0=10 Msun/pc^3) is too dense for
direct N-body. Using Schmidt et al. 2026 Appendix G rescaling symmetry, we map
to the simulation space (r_s=3.6 kpc, rho_0=7.09e-3 Msun/pc^3), which matches
the working reference HPC run that completed successfully.

The projected mass ratio M_2D(R_inner)/M_2D(R_outer) is INVARIANT under the
rescaling (both radii scale by lambda, both masses scale by mu, so ratio is
preserved). We can therefore:
  1. Run N-body in simulation space
  2. Compute the ratio at rescaled radii R_inner_sim = lambda * R_inner_phys
  3. Compare directly with the fluid model's ratio (also invariant)

Outputs a binary Gadget2 IC (ICFormat=1, GADGET2_HEADER) matching the format
the working reference uses.

Rescaling parameters:
  lambda = r_s_phys / r_s_sim = 0.085 / 3.6 = 0.023611
  mu     = (rho_0_phys / rho_0_sim) * lambda^3 = 0.018565

Physical radii -> Simulation radii (multiply by 1/lambda = 42.353):
  R_inner_phys = 0.020 kpc (20 pc) -> R_inner_sim = 0.847 kpc
  R_outer_phys = 0.090 kpc (90 pc) -> R_outer_sim = 3.812 kpc
  R_MAX_phys  = 0.340 kpc         -> R_MAX_sim  = 14.40 kpc (4*r_s_sim)

Simulation-space halo:
  r_s     = 3.6 kpc
  rho_0   = 7.09e-3 Msun/pc^3 = 7.09e6 Msun/kpc^3
  R_MAX   = 14.4 kpc (4 * r_s)
  M_total = 4*pi*rho_0*r_s^3 * (ln(5) - 4/5) ~ 3.42e8 Msun
  N       = 500000 (matches working reference's DM count)
  m_part  ~ 684 Msun
"""
import numpy as np
import struct
import os

# ============================================================
# Simulation-space halo parameters (matches working reference)
# ============================================================
R_S_KPC = 3.6
RHO_0_MSUN_PC3 = 7.09e-3
RHO_0_MSUN_KPC3 = RHO_0_MSUN_PC3 * 1e9  # 7.09e6 Msun/kpc^3
N_PART = 500000
R_MAX_KPC = 4.0 * R_S_KPC  # 14.4 kpc, matches 4*r_s

# Rescaling (for documentation / verification)
LAMBDA = 0.085 / R_S_KPC
MU = (10.0 / RHO_0_MSUN_PC3) * LAMBDA**3
R_INNER_PHYS = 0.020  # kpc
R_OUTER_PHYS = 0.090  # kpc
R_INNER_SIM = R_INNER_PHYS / LAMBDA   # 0.847 kpc
R_OUTER_SIM = R_OUTER_PHYS / LAMBDA   # 3.812 kpc

# Gadget4 code units (matching working reference params.txt)
# UnitLength = 1 kpc, UnitMass = 1e10 Msun, UnitVelocity = 1 km/s
# -> UnitTime = 1 kpc / 1 km/s = 0.978 Gyr
UNIT_LENGTH_KPC = 1.0
UNIT_mass_MSUN = 1.0e10
UNIT_VELOCITY_KM_S = 1.0
T_UNIT_GYR = 0.978

G_KPC = 4.302e-3  # kpc Msun^-1 (km/s)^2  (G in these units)

# ============================================================
# NFW profile
# ============================================================
def enclosed_mass(r):
    """M(<r) for NFW profile, in Msun. r in kpc."""
    x = r / R_S_KPC
    return 4 * np.pi * RHO_0_MSUN_KPC3 * R_S_KPC**3 * (np.log1p(x) - x / (1 + x))


def nfw_density(r):
    """NFW density in Msun/kpc^3. r in kpc."""
    x = r / R_S_KPC
    return RHO_0_MSUN_KPC3 / (x * (1 + x)**2)


# ============================================================
# Sample particle positions from NFW
# ============================================================
def sample_positions(n, r_max, seed=42):
    """Sample n positions uniformly from NFW profile within r_max."""
    rng = np.random.default_rng(seed)
    # Build CDF
    r_grid = np.linspace(1e-6, r_max, 20001)
    x_grid = r_grid / R_S_KPC
    M_grid = 4 * np.pi * RHO_0_MSUN_KPC3 * R_S_KPC**3 * (np.log1p(x_grid) - x_grid / (1 + x_grid))
    CDF = M_grid / M_grid[-1]
    # Invert
    u = rng.uniform(0, 1, n)
    r = np.interp(u, CDF, r_grid)
    # Isotropic direction
    phi = rng.uniform(0, 2 * np.pi, n)
    cos_theta = rng.uniform(-1, 1, n)
    sin_theta = np.sqrt(1 - cos_theta**2)
    x = r * sin_theta * np.cos(phi)
    y = r * sin_theta * np.sin(phi)
    z = r * cos_theta
    return np.column_stack([x, y, z]), r


def sample_velocities(r, seed=43):
    """Sample velocities from an isotropic Gaussian with sigma = v_circ/sqrt(3).

    For an NFW halo the exact solution of the Jeans equation would give a more
    accurate anisotropic dispersion profile, but for SIDM halos the initial
    velocities equilibrate quickly via scattering. We use the simple v_circ
    estimate matching the prior IC generator.
    """
    rng = np.random.default_rng(seed)
    M_r = 4 * np.pi * RHO_0_MSUN_KPC3 * R_S_KPC**3 * (
        np.log1p(r / R_S_KPC) - (r / R_S_KPC) / (1 + r / R_S_KPC))
    # Guard against r=0
    M_r = np.maximum(M_r, 1.0)
    r_safe = np.maximum(r, 1e-6)
    v_circ = np.sqrt(G_KPC * M_r / r_safe)
    sigma_v = v_circ / np.sqrt(3)
    vx = rng.normal(0, 1, len(r)) * sigma_v
    vy = rng.normal(0, 1, len(r)) * sigma_v
    vz = rng.normal(0, 1, len(r)) * sigma_v
    return np.column_stack([vx, vy, vz]), v_circ


# ============================================================
# Write Gadget2 binary IC (SnapFormat=1, GADGET2_HEADER)
# ============================================================
def make_gadget2_header(npart, mass_table, time, redshift, boxsize,
                        omega0, omega_l, hubble, num_files=1):
    """Build a 256-byte Gadget2 header. NTYPES=3 but arrays padded to 6."""
    npart6 = np.zeros(6, dtype=np.uint32)
    npart6[:3] = npart[:3]
    mass6 = np.zeros(6, dtype=np.float64)
    mass6[:3] = mass_table[:3]

    fields = [
        ('6I', npart6),
        ('6d', mass6),
        ('d', time),
        ('d', redshift),
        ('i', 0),  # flag_sfr
        ('i', 0),  # flag_feedback
        ('6I', npart6.copy()),  # npartTotal
        ('i', 0),  # flag_cooling
        ('i', num_files),
        ('d', boxsize),
        ('d', omega0),
        ('d', omega_l),
        ('d', hubble),
        ('i', 0),  # flag_stellarage
        ('i', 0),  # flag_metals
        ('6I', np.zeros(6, dtype=np.uint32)),  # npartTotalHighWord
        ('i', 0),  # flag_entropy
    ]
    buf = b''
    for fmt, data in fields:
        buf += struct.pack(fmt, *np.asarray(data).ravel())
    if len(buf) > 256:
        raise RuntimeError(f"Header too big: {len(buf)} bytes")
    buf += b'\x00' * (256 - len(buf))
    return buf


def write_ic_binary(filename, coords, vels, particle_ids, mass_per_particle,
                    boxsize, time=0.0):
    """Write Gadget2 binary IC (ICFormat=1).

    All particles are PartType1 (DM). Masses are set in the header (constant
    per type), so no Masses block is needed.
    """
    n = len(coords)
    npart = np.array([0, n, 0], dtype=np.uint32)
    mass_table = np.array([0.0, float(mass_per_particle), 0.0], dtype=np.float64)

    print(f"\nWriting binary IC: {filename}")
    print(f"  N = {n} particles (PartType1)")
    print(f"  Mass per particle = {mass_per_particle:.4f} Msun")
    print(f"  BoxSize = {boxsize:.4f} kpc")

    with open(filename, 'wb') as f:
        # Header block
        hdr = make_gadget2_header(npart, mass_table, time, 0.0, boxsize,
                                  0.0, 0.0, 1.0)
        f.write(struct.pack('I', 256))
        f.write(hdr)
        f.write(struct.pack('I', 256))

        # Block helper
        def write_block(name, data_bytes):
            nblk = len(data_bytes)
            f.write(struct.pack('I', nblk))
            f.write(data_bytes)
            f.write(struct.pack('I', nblk))
            print(f"  {name}: {nblk} bytes")

        # Coordinates (N, 3) float32
        coords_f32 = coords.astype(np.float32)
        write_block("Coordinates", coords_f32.tobytes())

        # Velocities (N, 3) float32
        vels_f32 = vels.astype(np.float32)
        write_block("Velocities", vels_f32.tobytes())

        # ParticleIDs (N,) uint32
        pids_u32 = particle_ids.astype(np.uint32)
        write_block("ParticleIDs", pids_u32.tobytes())

    fsize = os.path.getsize(filename)
    print(f"  Total file size: {fsize/1e6:.2f} MB")
    return fsize


# ============================================================
# Main
# ============================================================
def main():
    print("=" * 70)
    print("Simulation-space NFW IC generator (rescaling symmetry)")
    print("=" * 70)
    print(f"\nHalo: r_s = {R_S_KPC} kpc, rho_0 = {RHO_0_MSUN_PC3} Msun/pc^3")
    print(f"R_MAX = {R_MAX_KPC} kpc (4 * r_s)")
    print(f"N = {N_PART} particles")
    print(f"\nRescaling: lambda = {LAMBDA:.6f}, mu = {MU:.6f}")
    print(f"  R_inner_sim = {R_INNER_SIM:.4f} kpc  (from 20 pc)")
    print(f"  R_outer_sim = {R_OUTER_SIM:.4f} kpc  (from 90 pc)")

    M_total = enclosed_mass(R_MAX_KPC)
    m_part = M_total / N_PART
    print(f"\nM_total(<{R_MAX_KPC} kpc) = {M_total:.4e} Msun")
    print(f"Mass per particle = {m_part:.4f} Msun")

    # Sample positions
    print("\nSampling positions...")
    coords, r = sample_positions(N_PART, R_MAX_KPC, seed=42)
    print(f"  r range: [{r.min():.4f}, {r.max():.4f}] kpc")
    print(f"  x range: [{coords[:,0].min():.4f}, {coords[:,0].max():.4f}]")

    # Sample velocities
    print("Sampling velocities...")
    vels, v_circ = sample_velocities(r, seed=43)
    print(f"  v_circ range: [{v_circ.min():.2f}, {v_circ.max():.2f}] km/s")
    print(f"  |v| range: [{np.linalg.norm(vels, axis=1).min():.2f}, "
          f"{np.linalg.norm(vels, axis=1).max():.2f}] km/s")

    # Particle IDs
    pids = np.arange(1, N_PART + 1, dtype=np.uint32)

    # Write binary IC
    out_dir = "D:/graverthermal-sidm/data/P5_nbody_verify/ics_sim_space"
    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, "ic.dat")
    boxsize = 2 * R_MAX_KPC  # 28.8 kpc
    write_ic_binary(out_file, coords, vels, pids, m_part, boxsize, time=0.0)

    # ---- Verify: compute projected mass ratio in simulation space ----
    print("\n" + "=" * 70)
    print("Verification: projected mass ratio in simulation space")
    print("=" * 70)
    # Use unit-mass particles (since all have same m_part)
    masses = np.full(N_PART, m_part)

    # Project onto xy plane (z is line of sight)
    r_proj_xy = np.sqrt(coords[:, 0]**2 + coords[:, 1]**2)
    m_inner_xy = masses[r_proj_xy <= R_INNER_SIM].sum()
    m_outer_xy = masses[r_proj_xy <= R_OUTER_SIM].sum()
    ratio_xy = m_inner_xy / m_outer_xy if m_outer_xy > 0 else np.nan

    # Average over 3 axes
    r_proj_xz = np.sqrt(coords[:, 0]**2 + coords[:, 2]**2)
    r_proj_yz = np.sqrt(coords[:, 1]**2 + coords[:, 2]**2)
    m_inner_xz = masses[r_proj_xz <= R_INNER_SIM].sum()
    m_outer_xz = masses[r_proj_xz <= R_OUTER_SIM].sum()
    m_inner_yz = masses[r_proj_yz <= R_INNER_SIM].sum()
    m_outer_yz = masses[r_proj_yz <= R_OUTER_SIM].sum()
    ratio_avg = ((m_inner_xy / m_outer_xy) +
                 (m_inner_xz / m_outer_xz) +
                 (m_inner_yz / m_outer_yz)) / 3.0

    print(f"  R_inner_sim = {R_INNER_SIM:.4f} kpc")
    print(f"  R_outer_sim = {R_OUTER_SIM:.4f} kpc")
    print(f"  M_2D(<R_inner, xy)   = {m_inner_xy:.4e} Msun  (N={np.sum(r_proj_xy <= R_INNER_SIM)})")
    print(f"  M_2D(<R_outer, xy)   = {m_outer_xy:.4e} Msun  (N={np.sum(r_proj_xy <= R_OUTER_SIM)})")
    print(f"  Ratio (xy)           = {ratio_xy:.6f}")
    print(f"  Ratio (avg 3 axes)   = {ratio_avg:.6f}")
    print(f"\n  Fluid model ratio_init = 0.1465 (target)")
    print(f"  N-body IC ratio       = {ratio_avg:.4f}")
    print(f"  Difference            = {(ratio_avg - 0.1465)/0.1465*100:.2f}%")

    # Also verify the density profile matches NFW
    print("\nDensity profile check:")
    r_bins = np.array([0.1, 0.5, 1.0, 2.0, 3.6, 5.0, 7.2, 10.0, 14.4])
    r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
    for i in range(len(r_bins) - 1):
        mask = (r >= r_bins[i]) & (r < r_bins[i+1])
        n_in = np.sum(mask)
        vol = (4/3) * np.pi * (r_bins[i+1]**3 - r_bins[i]**3)
        rho_nbody = n_in * m_part / vol
        rho_nfw = nfw_density(r_centers[i])
        print(f"  r=[{r_bins[i]:.1f}-{r_bins[i+1]:.1f}] kpc: "
              f"rho_nbody={rho_nbody:.3e}, rho_nfw={rho_nfw:.3e}, "
              f"ratio={rho_nbody/rho_nfw:.3f}")

    print(f"\nIC file written to: {out_file}")
    print("Done.")


if __name__ == '__main__':
    main()
