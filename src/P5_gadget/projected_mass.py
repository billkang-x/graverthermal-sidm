#!/usr/bin/env python3
"""
Compute the 2D PROJECTED enclosed mass from N-body particles, matching the
fluid model's `projected_enclosed_mass` definition.

For N-body, the projected mass within projected radius R_2d is computed by
projecting particles onto a plane (say xy) and counting mass within R_2d
in that plane.

  M_2D(R) = Σ_i m_i  for all particles with sqrt(x_i^2 + y_i^2) <= R

This naturally includes the line-of-sight integration. We should average over
3 orthogonal projections (xy, yz, xz) for better statistics.

This must match the fluid model's projected_enclosed_mass, which integrates
ρ(r) along the line of sight analytically.
"""
import numpy as np


def projected_mass_2d(coords, masses, r_2d, axis='z'):
    """Compute projected enclosed mass within projected radius r_2d.

    Projects particles onto the plane perpendicular to `axis`, then sums
    masses within cylindrical radius r_2d.

    Args:
        coords: (N, 3) array of particle positions
        masses: (N,) array of particle masses
        r_2d: projected radius
        axis: 'x', 'y', or 'z' - the line-of-sight axis

    Returns:
        M_2d (scalar): total mass within projected radius r_2d
    """
    if axis == 'z':
        r_proj = np.sqrt(coords[:, 0]**2 + coords[:, 1]**2)
    elif axis == 'y':
        r_proj = np.sqrt(coords[:, 0]**2 + coords[:, 2]**2)
    elif axis == 'x':
        r_proj = np.sqrt(coords[:, 1]**2 + coords[:, 2]**2)
    else:
        raise ValueError(f"Unknown axis: {axis}")
    return masses[r_proj <= r_2d].sum()


def projected_mass_2d_avg(coords, masses, r_2d):
    """Average projected mass over 3 orthogonal projections.

    For an isotropic halo this gives the best estimate of M_2D(r_2d).
    """
    m_z = projected_mass_2d(coords, masses, r_2d, axis='z')
    m_y = projected_mass_2d(coords, masses, r_2d, axis='y')
    m_x = projected_mass_2d(coords, masses, r_2d, axis='x')
    return (m_z + m_y + m_x) / 3.0


def projected_mass_ratio(coords, masses, r_inner, r_outer, avg_axes=True):
    """Compute M_2D(r_inner) / M_2D(r_outer).

    Args:
        coords: (N, 3) positions
        masses: (N,) masses
        r_inner: inner projected radius (e.g., 0.020 kpc = 20 pc)
        r_outer: outer projected radius (e.g., 0.090 kpc = 90 pc)
        avg_axes: if True, average over 3 projections

    Returns:
        ratio (scalar), m_inner, m_outer
    """
    if avg_axes:
        m_inner = projected_mass_2d_avg(coords, masses, r_inner)
        m_outer = projected_mass_2d_avg(coords, masses, r_outer)
    else:
        m_inner = projected_mass_2d(coords, masses, r_inner, axis='z')
        m_outer = projected_mass_2d(coords, masses, r_outer, axis='z')
    if m_outer == 0:
        return np.nan, m_inner, m_outer
    return m_inner / m_outer, m_inner, m_outer


if __name__ == '__main__':
    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from read_binary_snap import read_snapshot

    if len(sys.argv) < 2:
        print("Usage: projected_mass.py <snapshot_file>")
        sys.exit(1)

    snap = read_snapshot(sys.argv[1], types_to_read=(1,))
    coords = snap['coords'][1]
    masses = snap['masses'][1]
    print(f"N = {len(coords)}, M_total = {masses.sum():.4e}")

    r_inner = 0.020  # 20 pc
    r_outer = 0.090  # 90 pc
    ratio, m_in, m_out = projected_mass_ratio(coords, masses, r_inner, r_outer)
    print(f"\nProjected mass (avg of 3 axes):")
    print(f"  M_2D(<{r_inner*1000:.0f}pc) = {m_in:.4e}")
    print(f"  M_2D(<{r_outer*1000:.0f}pc) = {m_out:.4e}")
    print(f"  Ratio = {ratio:.6f}")

    # Also compute 3D spherical for comparison
    r3d = np.sqrt((coords**2).sum(axis=1))
    m3d_in = masses[r3d <= r_inner].sum()
    m3d_out = masses[r3d <= r_outer].sum()
    print(f"\n3D spherical mass (for comparison):")
    print(f"  M_3D(<{r_inner*1000:.0f}pc) = {m3d_in:.4e}")
    print(f"  M_3D(<{r_outer*1000:.0f}pc) = {m3d_out:.4e}")
    print(f"  Ratio = {m3d_in/m3d_out:.6f}")
