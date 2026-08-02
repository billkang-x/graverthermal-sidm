"""P1 figure: sigma_T(v) and r_diss(v) for three benchmark dSIDM models.

Demonstrates the key physical point: massive emission introduces a
velocity scale v* = sqrt(2 m_mediator / mu) that breaks the (lambda, mu)
rescaling symmetry of Schmidt et al. 2026.  Massless emission (control)
recovers their constant-r_diss ansatz.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from cross_sections.dsidm_models import benchmark_models, sigma_T_born, r_diss

plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 12,
    'legend.fontsize': 10,
    'figure.figsize': (12, 5),
    'figure.dpi': 120,
})

models = benchmark_models()
v = np.logspace(1, 3.5, 200)  # 10 to ~3000 km/s

# Characteristic velocity scales
v_dwarf = 50     # km/s
v_MW = 200
v_cluster = 1000

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

colors = {'M1_dark_photon_massive': '#993C1D',
          'M2_scalar_phi_massive': '#185FA5',
          'M3_massless_control': '#0F6E56'}
labels = {'M1_dark_photon_massive': r'M1: dark photon $V$ (massive, $m_V=1$ MeV)',
          'M2_scalar_phi_massive': r'M2: scalar $\phi$ (massive, $m_\phi=1$ MeV)',
          'M3_massless_control': r'M3: massless $V$ (control)'}

for name, p in models.items():
    sig = sigma_T_born(v, p)
    rd = r_diss(v, p)
    ax1.loglog(v, sig, color=colors[name], label=labels[name], lw=2)
    ax2.semilogx(v, rd, color=colors[name], label=labels[name], lw=2)

# Observational constraint bands
# Dwarf galaxies: sigma/m ~ 0.1-50 cm^2/g at v~50 km/s
ax1.axvspan(10, 100, alpha=0.08, color='gray', label='dwarf scale')
ax1.axvspan(100, 500, alpha=0.05, color='gray', label='MW scale')
ax1.axvspan(500, 3000, alpha=0.03, color='gray', label='cluster scale')

# Reference lines for sigma constraints
ax1.axhline(1.25, color='gray', ls=':', alpha=0.5, lw=1)
ax1.axhline(0.1, color='gray', ls=':', alpha=0.5, lw=1)
ax1.text(2000, 1.4, 'Bullet Cluster', fontsize=9, color='gray')
ax1.text(2000, 0.13, 'cluster cores', fontsize=9, color='gray')

ax1.set_xlabel(r'relative velocity $v$ [km/s]')
ax1.set_ylabel(r'$\sigma_T / m_\chi$ [cm$^2$/g]')
ax1.set_title('Transfer cross section (velocity-dependent)')
ax1.set_xlim(10, 3500)
ax1.set_ylim(1e-3, 1e4)
ax1.legend(loc='lower left', fontsize=9)
ax1.grid(True, alpha=0.3, which='both')

# r_diss panel: mark Schmidt et al. explored range
ax2.axhspan(1.01, 1.3, alpha=0.06, color='#993C1D', label='Schmidt et al. range')
ax2.axvline(v_dwarf, color='gray', ls='--', alpha=0.4)
ax2.axvline(v_MW, color='gray', ls='--', alpha=0.4)
ax2.axvline(v_cluster, color='gray', ls='--', alpha=0.4)
ax2.text(v_dwarf*1.1, 1.001, 'dwarf', fontsize=9, color='gray')
ax2.text(v_MW*1.1, 1.001, 'MW', fontsize=9, color='gray')
ax2.text(v_cluster*1.1, 1.001, 'cluster', fontsize=9, color='gray')

# Mark the symmetry-breaking velocity scale for massive models
for name in ['M1_dark_photon_massive', 'M2_scalar_phi_massive']:
    p = models[name]
    v_star = np.sqrt(2 * p.m_mediator / p.mu) / (1.0/299792.458)
    ax2.axvline(v_star, color=colors[name], ls=':', alpha=0.5, lw=1)
    ax2.text(v_star*1.05, 1.045, r'$v_*$', fontsize=10, color=colors[name])

ax2.set_xlabel(r'relative velocity $v$ [km/s]')
ax2.set_ylabel(r'$r_{\rm diss}(v)$')
ax2.set_title('Dissipation parameter (symmetry-breaking)')
ax2.set_xlim(10, 3500)
ax2.set_ylim(0.998, 1.07)
ax2.legend(loc='upper left', fontsize=9)
ax2.grid(True, alpha=0.3, which='both')

plt.suptitle('P1: Three dSIDM models — velocity dependence breaks rescaling symmetry',
             fontsize=13, y=1.02)
plt.tight_layout()

out = os.path.join(os.path.dirname(__file__), '..', '..', 'figures', 'P1_sigma_r_diss.png')
plt.savefig(out, bbox_inches='tight', dpi=150)
print(f'Saved: {out}')

# Also dump numerical values to CSV
import csv
csv_path = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'P1_cross_sections.csv')
with open(csv_path, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['v_km_s'] + [f'{n}_sigma_cm2g' for n in models] + [f'{n}_r_diss' for n in models])
    for vi in v:
        row = [vi]
        for n, p in models.items():
            row.append(sigma_T_born(np.array([vi]), p)[0])
        for n, p in models.items():
            row.append(r_diss(np.array([vi]), p)[0])
        w.writerow(row)
print(f'Saved: {csv_path}')
