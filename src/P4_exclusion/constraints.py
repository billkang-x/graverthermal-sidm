"""
Observational constraints on SIDM cross sections.

Each constraint is characterized by:
- A characteristic velocity scale v_char (km/s) at which σ_T/m is probed
- An upper bound (or range) on σ_T/m at that velocity (cm²/g)
- A label, color, and reference

Sources (canonical, used widely in the SIDM exclusion-plot literature):
- Bullet Cluster (Markevitch et al. 2002; Randall et al. 2008):
    σ_T/m < 1.25 cm²/g at v ~ 3000 km/s (90% CL)
    Often quoted as the cleanest cluster-scale constraint.
- Cluster cores (Kaplinghat, Tulin, Yu 2016; Andrade et al. 2022):
    σ_T/m ≲ 1 cm²/g at v ~ 1000-1500 km/s
    From halo shapes and core densities in massive clusters.
- Dwarf cores / too-late-to-form issues (Read et al. 2018; Sokolenko et al. 2018;
    Kahlhoefer et al. 2019):
    σ_T/m ~ 1-50 cm²/g preferred at v ~ 30-100 km/s to explain cores
    Upper bound ~ 50 cm²/g from dwarfs that remain cuspy.
- B1938+666 (Vegetti et al. 2026; interpreted by Schmidt et al. 2026):
    Requires σ_T/m up to ~7.8 cm²/g at v ~ 50-200 km/s
    Lower bound from requiring fast-enough collapse.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Optional


@dataclass
class ObservationalConstraint:
    name: str                # short name
    label: str               # plot label (LaTeX)
    v_char_kms: float        # characteristic velocity (km/s)
    v_lo_kms: float          # lower velocity of applicability
    v_hi_kms: float          # upper velocity of applicability
    sigma_upper: Optional[float]  # upper bound σ/m (cm²/g), None if not upper-bounded
    sigma_lower: Optional[float]  # lower bound σ/m (cm²/g), None if not lower-bounded
    color: str
    linestyle: str
    ref: str


CONSTRAINTS = [
    ObservationalConstraint(
        name='bullet',
        label=r'Bullet Cluster',
        v_char_kms=3000., v_lo_kms=2000., v_hi_kms=4500.,
        sigma_upper=1.25, sigma_lower=None,
        color='darkred', linestyle='-',
        ref='Randall et al. 2008',
    ),
    ObservationalConstraint(
        name='cluster_cores',
        label=r'Cluster cores',
        v_char_kms=1200., v_lo_kms=800., v_hi_kms=1800.,
        sigma_upper=1.0, sigma_lower=None,
        color='purple', linestyle='-',
        ref='Kaplinghat et al. 2016',
    ),
    ObservationalConstraint(
        name='dwarf_cores',
        label=r'Dwarf cores',
        v_char_kms=50., v_lo_kms=20., v_hi_kms=120.,
        sigma_upper=50., sigma_lower=1.0,
        color='teal', linestyle='-',
        ref='Read et al. 2018',
    ),
    ObservationalConstraint(
        name='b1938',
        label=r'B1938+666',
        v_char_kms=80., v_lo_kms=30., v_hi_kms=200.,
        sigma_upper=7.8, sigma_lower=0.07,
        color='orange', linestyle='-',
        ref='Schmidt et al. 2026',
    ),
]


def get_constraint(name: str) -> ObservationalConstraint:
    for c in CONSTRAINTS:
        if c.name == name:
            return c
    raise KeyError(f"Unknown constraint: {name}")


if __name__ == '__main__':
    print("Observational constraints:")
    print(f"  {'name':<15} {'v_char [km/s]':<15} "
          f"{'σ/m upper':<12} {'σ/m lower':<12} {'label'}")
    for c in CONSTRAINTS:
        sup = f"{c.sigma_upper}" if c.sigma_upper is not None else "—"
        slo = f"{c.sigma_lower}" if c.sigma_lower is not None else "—"
        print(f"  {c.name:<15} {c.v_char_kms:<15.0f} "
              f"{sup:<12} {slo:<12} {c.label}  ({c.ref})")
