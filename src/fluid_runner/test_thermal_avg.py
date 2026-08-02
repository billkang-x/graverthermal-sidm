"""Test thermal averaging with the velocity-dependent DSIDM models from P1."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'cross_sections'))
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from dsidm_models import benchmark_models, sigma_T_born, r_diss
from thermal_avg import effective_sigma_m_and_rdiss

models = benchmark_models()

# Temperature axis: T = v_1d^2, covering dwarf (v~30 km/s) to cluster (v~1000 km/s) scales
T_axis = np.logspace(1.5, 6.5, 40)  # (km/s)^2  ~ 30 to 3000 km/s 1D

print("=" * 80)
print("Thermal averages for the three benchmark DSIDM models")
print("=" * 80)
print(f"{'Model':<30} {'T [km^2/s^2]':<14} {'<sig_m> [cm^2/g]':<20} {'<r_diss>':<10}")
print("-" * 80)

for name, p in models.items():
    # Build callables that take v in km/s and return sigma_T/m (cm^2/g) and r_diss
    sig_fn = lambda v, p=p: sigma_T_born(np.atleast_1d(v), p)
    rd_fn = lambda v, p=p: r_diss(np.atleast_1d(v), p)

    sigma_m_eff, rdiss_eff, sig_arr, rd_arr = effective_sigma_m_and_rdiss(
        sig_fn, rd_fn, T_axis
    )

    # Report at representative astrophysical temperatures
    for T_rep, label in [(1e3, 'dwarf~30'), (1e4, 'dwarf~100'), (1e5, 'MW~300'), (1e6, 'cluster~1000')]:
        sm = sigma_m_eff(np.array([T_rep]))[0]
        rd = rdiss_eff(np.array([T_rep]))[0]
        print(f"{name:<30} {T_rep:<14.2e} {sm:<20.4e} {rd:<10.4f}")
    print()
