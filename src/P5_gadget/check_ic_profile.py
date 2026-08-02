#!/usr/bin/env python3
"""Check the actual radial density profile of the IC particles vs NFW."""
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from read_binary_snap import read_snapshot

R_S = 0.085  # kpc
RHO_0 = 10.0  # Msun/pc^3

def rho_nfw(r_kpc):
    x = r_kpc / R_S
    return RHO_0 / (x * (1+x)**2)

ic_path = "D:/graverthermal-sidm/data/P5_nbody_verify/ics/ic.dat"
snap = read_snapshot(ic_path, types_to_read=(1,))
coords = snap['coords'][1]
masses = snap['masses'][1]
r = np.sqrt((coords**2).sum(axis=1))
N = len(r)
M_total = masses.sum()
m_part = M_total / N
print(f"N={N}, M_total={M_total:.4e}, m_part={m_part:.4e}")
print(f"r range: {r.min():.5f} - {r.max():.5f} kpc")

# Compute density profile in log bins
r_bins = np.logspace(np.log10(r.min()), np.log10(r.max()), 51)
r_centers = np.sqrt(r_bins[:-1] * r_bins[1:])  # geometric mean
dr = np.diff(r_bins)
rho_ic = np.zeros(len(r_centers))
for i in range(len(r_centers)):
    mask = (r >= r_bins[i]) & (r < r_bins[i+1])
    n_in_bin = mask.sum()
    if n_in_bin > 0:
        M_in_bin = masses[mask].sum()
        V_shell = (4/3) * np.pi * (r_bins[i+1]**3 - r_bins[i]**3)
        rho_ic[i] = M_in_bin / V_shell  # Msun/kpc^3
        rho_ic_msun_pc3 = rho_ic[i] / 1e9
    else:
        rho_ic[i] = np.nan

# Convert to Msun/pc^3 for comparison
rho_ic_pc3 = rho_ic / 1e9
rho_nfw_pc3 = rho_nfw(r_centers)

print("\nRadial density profile:")
print(f"{'r (kpc)':>10} {'r (pc)':>10} {'rho_IC (Msun/pc3)':>20} {'rho_NFW (Msun/pc3)':>20} {'ratio':>10}")
for i in range(len(r_centers)):
    if not np.isnan(rho_ic_pc3[i]):
        ratio = rho_ic_pc3[i] / rho_nfw_pc3[i] if rho_nfw_pc3[i] > 0 else np.nan
        print(f"{r_centers[i]:>10.5f} {r_centers[i]*1000:>10.1f} {rho_ic_pc3[i]:>20.4e} {rho_nfw_pc3[i]:>20.4e} {ratio:>10.4f}")

# Plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
ax1.loglog(r_centers, rho_ic_pc3, 'bo-', label='IC (N-body)', markersize=4)
ax1.loglog(r_centers, rho_nfw_pc3, 'r--', label='NFW analytic', linewidth=2)
ax1.set_xlabel('r (kpc)')
ax1.set_ylabel('rho (Msun/pc^3)')
ax1.set_title('IC density profile vs NFW')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Cumulative mass
r_sorted = np.sort(r)
m_sorted = masses[np.argsort(r)]
M_cum = np.cumsum(m_sorted)
ax2.semilogx(r_sorted, M_cum / M_total, 'b-', label='IC cumulative M(<r)/M_total')
# NFW cumulative
r_nfw = np.logspace(np.log10(r.min()), np.log10(r.max()), 200)
x = r_nfw / R_S
M_nfw = np.log1p(x) - x/(1+x)
M_nfw /= M_nfw[-1]
ax2.semilogx(r_nfw, M_nfw, 'r--', label='NFW cumulative', linewidth=2)
ax2.set_xlabel('r (kpc)')
ax2.set_ylabel('M(<r) / M_total')
ax2.set_title('Cumulative mass profile')
ax2.legend()
ax2.grid(True, alpha=0.3)

plt.tight_layout()
out_png = "D:/graverthermal-sidm/data/P5_nbody_verify/ic_density_profile.png"
plt.savefig(out_png, dpi=150)
print(f"\nSaved plot to {out_png}")
