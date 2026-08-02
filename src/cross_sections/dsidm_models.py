"""
Velocity-dependent sigma_T(v) and r_diss(v) for dissipative SIDM models.

Based on Lankester-Broche & Pradler 2026 (JCAP 03, 034; arXiv:2509.12317),
"Towards a theory of dissipative Dark Matter I: the Born limit".

Six dSIDM scenarios unified via effective coupling a_j:
    1. chi_1 chi_2 -> chi_1 chi_2 V    (Dirac fermion, vector mediator+emission)
    2. S_1 S_2 -> S_1 S_2 V            (complex scalar, vector)
    3. chi_1 chi_2 -> chi_1 chi_2 phi  (Dirac fermion, scalar)
    4. chi~_1 chi~_2 -> chi~_1 chi~_2 phi  (Majorana fermion, scalar)
    5. S_1 S_2 -> S_1 S_2 phi          (complex scalar, scalar)
    6. S~_1 S~_2 -> S~_1 S~_2 phi      (real scalar, scalar)

Key equations (all from LB2026):
    Eq (2.6):   sigma_T(v_i)  [elastic, Born]
    Eq (3.21):  <|M_D|^2>     [dipole emission amplitude, massless]
    Eq (3.26):  C_phi = 1 - m_phi^2/omega^2     [massive scalar correction]
    Eq (3.27):  C_V   = 1 + m_V^2/(2 omega^2)   [massive vector correction]
    Eq (3.28):  omega d sigma / d omega           [energy-differential, dipole, long-range]
    Eq (5.6):   epsilon_dot (dipole, long-range, massless)  ~ T^{1/2}
    Eq (5.8):   epsilon_dot (dipole, short-range, massless) ~ T^{5/2}

Mapping to Schmidt et al. 2026:
    sigma_T/m_chi      <- Eq (2.6) / m_chi
    r_diss(v) - 1      <- epsilon_dot / (elastic scattering rate * <k0>)

The elastic scattering rate per particle ~ n <sigma_T v>.
The mean radiated energy <k0> per collision is obtained from
    <k0> = (1/sigma_T) * integral[omega * (d sigma / d omega) d omega]
which for dipole long-range emission yields <k0> ~ mu v_i^2 (i.e., a fixed
fraction of the CM kinetic energy in the massless case, consistent with
Schmidt et al.'s constant-r_diss ansatz).  Massive emission introduces a
velocity-dependent <k0>(v), which breaks the rescaling symmetry.

Units:
    velocities in km/s
    masses in GeV
    sigma_T/m_chi in cm^2/g
    r_diss dimensionless
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from typing import Callable, Optional

# Physical constants (c = hbar = 1 convention; conversions for output units)
# 1 GeV^-2 = 0.3894e-27 cm^2
# 1 GeV = 1.7827e-24 g
# 1 km/s = 1/299792.458 (natural units, c=1)
GEV2_TO_CM2G = 0.3894e-27 / 1.7827e-24  # GeV^-2 -> cm^2/g conversion factor
KM_S_TO_NAT = 1.0 / 299792.458           # km/s -> natural units (c=1)
K_BOLTZMANN = 8.617333e-14               # GeV/K (not used; T expressed as energy)
# Provisional amplitude retained only for the legacy phenomenological
# ``r_diss`` compatibility closure.  Quantitative cooling uses
# ``cross_sections.emission_kernel`` instead.
DEFAULT_C0 = 0.05


@dataclass
class DSIDMParameters:
    """Parameters for a single dSIDM model realization.

    The six scenarios of LB2026 are unified through the effective coupling a_j.
    For identical-particle SIDM (m1 = m2 = m_chi), dipole emission vanishes
    (charge-to-mass ratios equal) and quadrupole dominates.  For the realistic
    case of SIDM halos we consider self-scattering of identical DM particles,
    so quadrupole is the relevant channel -- BUT Schmidt et al.'s r_diss
    parametrization is general; we keep dipole here for the distinguishable
    case (e.g. two-component atomic/mirror DM) and note quadrupole scaling.

    Parameters
    ----------
    model : str
        One of {'chi-V', 'S-V', 'chi-phi', 'chi_tilde-phi', 'S-phi', 'S_tilde-phi'}.
    m_chi : float
        DM mass in GeV.
    m_mediator : float
        Light mediator/emitted particle mass in GeV (m_V or m_phi).
        For massless emission set to 0.
    m_mediator_heavy : float, optional
        Heavy (short-range) mediator mass m_{phi',V'} in GeV.
        If None, long-range mediation is assumed.
    alpha_D : float
        Dark fine-structure constant (vector) or effective scalar coupling.
        For vector: g = sqrt(4*pi*alpha_D).
        For scalar: y (Yukawa) or A/(2m) (trilinear/mass).
    mediation : str
        'long' or 'short'.
    emission_type : str
        'massless' or 'massive'.
    """
    model: str
    m_chi: float
    m_mediator: float
    alpha_D: float
    m_mediator_heavy: Optional[float] = None
    mediation: str = 'long'
    emission_type: str = 'massless'

    def __post_init__(self):
        assert self.model in {
            'chi-V', 'S-V', 'chi-phi', 'chi_tilde-phi', 'S-phi', 'S_tilde-phi'
        }, f"Unknown model: {self.model}"
        assert self.mediation in ('long', 'short')
        assert self.emission_type in ('massless', 'massive')

    @property
    def g_eff(self) -> float:
        """Effective coupling a_j for identical particles (a1 = a2 = a)."""
        if self.model in ('chi-V', 'S-V'):
            return np.sqrt(4 * np.pi * self.alpha_D)
        elif self.model in ('chi-phi', 'chi_tilde-phi'):
            # Yukawa; alpha_D reinterpreted as y^2/(4*pi)
            return np.sqrt(4 * np.pi * self.alpha_D)
        elif self.model in ('S-phi', 'S_tilde-phi'):
            # A/(2m) = sqrt(4*pi*alpha_D)  (parametric mapping)
            return np.sqrt(4 * np.pi * self.alpha_D)

    @property
    def f_phiV(self) -> float:
        """Polarization factor: 1 for scalar, 2 for vector emission."""
        return 2.0 if 'V' in self.model else 1.0

    @property
    def mu(self) -> float:
        """Reduced mass for identical particles: mu = m_chi/2."""
        return self.m_chi / 2.0


def sigma_T_born(v_km_s: np.ndarray, p: DSIDMParameters) -> np.ndarray:
    """Eq. (2.6) of LB2026: Born elastic transfer cross section.

    sigma_T = (a^4) / (8 pi mu^2 v^4) * [ln(1 + 4 mu^2 v^2 / m^2) - ...]

    For identical particles a1 = a2 = a = g_eff.

    Mediator mass m is the one that mediates the ELASTIC 2->2 scattering:
        long-range  : m = m_mediator (light, same as emitted particle)
        short-range : m = m_mediator_heavy (heavy, distinct from emitted)

    Returns sigma_T / m_chi in cm^2/g.
    """
    v = np.asarray(v_km_s, dtype=float) * KM_S_TO_NAT  # natural units (c=1)
    a = p.g_eff
    mu = p.mu

    if p.mediation == 'long':
        m = p.m_mediator
    else:
        m = p.m_mediator_heavy

    if m is None or m == 0:
        # Massless mediator: Rutherford, log-divergent; regulate with a soft
        # mediator-to-reduced-mass ratio of 1e-8.
        m_reg = 1e-8 * mu
        x = 4 * mu**2 * v**2 / m_reg**2
    else:
        x = 4 * mu**2 * v**2 / m**2

    x = np.clip(x, 1e-30, None)
    bracket = np.log1p(x) - x / (1 + x)

    sigma_nat = a**4 / (8 * np.pi * mu**2 * v**4) * bracket
    sigma_over_m = sigma_nat * GEV2_TO_CM2G / p.m_chi
    return sigma_over_m


def born_expansion_parameter(p: DSIDMParameters) -> float:
    """Return alpha_D * mu / m_med for the elastic interaction.

    The Yukawa Born expression requires this parameter to be much smaller
    than one.  A value below one is used only as a lenient plotting mask; it
    should not be read as a precision-error guarantee.
    """
    mediator_mass = (
        p.m_mediator if p.mediation == "long" else p.m_mediator_heavy
    )
    if mediator_mass is None or mediator_mass <= 0:
        return np.inf
    return float(p.alpha_D * p.mu / mediator_mass)


def is_born_valid(p: DSIDMParameters, threshold: float = 1.0) -> bool:
    """Return whether a parameter point passes the lenient Born mask."""
    if threshold <= 0:
        raise ValueError("threshold must be positive")
    return born_expansion_parameter(p) < threshold


def dipole_emission_factor(p: DSIDMParameters) -> float:
    """F_D for identical particles vanishes (charge-to-mass ratios equal).

    For distinguishable particles (e.g. m1 != m2 in two-component DM),
    F_D = a^2 * (a/m1 - a/m2).  Here we expose the structure for the
    distinguishable case and return 0 for identical (quadrupole-dominated).

    In practice, Schmidt et al.'s r_diss parametrization is model-agnostic;
    we treat r_diss as derived from the energy loss rate, which for identical
    particles uses the quadrupole channel (T^{3/2} scaling, eq 5.6 second term).
    """
    # Identical particles: dipole cancels
    return 0.0


def energy_loss_rate_long_range(T: float, n: float, p: DSIDMParameters) -> float:
    """Eq. (5.6) of LB2026, dipole+quadrupole, long-range mediation, massless.

    For identical particles (dipole=0), only quadrupole survives:
        eps_dot = n^2 * F_Q^2 * (24/45) * (mu*T)^{3/2} * (4*alpha+beta)
                  / (6 pi^{7/2} sqrt(2) * f_phiV)

    Here T is the temperature in GeV (= nu^2, 1D velocity dispersion in nat units),
    n is number density in GeV^3.

    Returns energy loss rate per volume in GeV^4 (natural units).
    """
    # For identical fermion DM (chi-V, chi-phi), quadrupole dominates.
    # 4*alpha + beta is O(1) model-dependent; use fiducial value 1.
    model_factor = 1.0  # placeholder; to be refined per model from Sec 3.5
    F_Q_eff = p.g_eff**2 * p.m_chi  # parametric F_Q ~ a^2 * m for identical

    eps = (n**2 / p.f_phiV) * (1.0 / (6 * np.pi**3.5 * np.sqrt(2))) * \
         (F_Q_eff**2 * (24.0/45.0) * (p.mu * T)**1.5 * model_factor)
    return eps


def mean_radiated_energy(v_km_s: float, p: DSIDMParameters) -> float:
    """<k0> per collision.

    For massless emission: <k0> = c * mu * v^2 / 2  (fixed fraction of CM KE).
        Coefficient c ~ O(0.1-1) from integrating eq (3.28).
    For massive emission: <k0>(v) is velocity-dependent; below the threshold
        v* = sqrt(2 m_mediator / mu), emission is exponentially suppressed.
        This is the symmetry-breaking term.
    """
    v = v_km_s * KM_S_TO_NAT
    mu = p.mu
    KE_cm = 0.5 * mu * v**2

    if p.emission_type == 'massless' or p.m_mediator == 0:
        # Fixed fraction (coefficient from dipole integral ~ 0.3-0.5)
        return 0.4 * KE_cm
    else:
        # Massive emission: phase-space suppression factor
        # eps_dot(massive) / eps_dot(massless) ~ exp(-m/T) for m >> T
        # Per-collision: <k0> ~ KE_cm * (1 - m/KE_cm)^{3/2} for KE_cm > m
        m = p.m_mediator
        if KE_cm <= m:
            return 1e-30 * KE_cm  # exponentially suppressed
        x = m / KE_cm
        # Approximate form (full requires numerical integration of eq 3.28 + 3.26/3.27)
        suppression = (1 - x)**1.5 * np.exp(-x)
        return 0.4 * KE_cm * suppression


def r_diss(v_km_s: np.ndarray, p: DSIDMParameters) -> np.ndarray:
    """Effective dissipation parameter r_diss(v) as defined by Schmidt et al. 2026.

    r_diss - 1 = (energy loss rate) / (elastic scattering rate * <k0>)
               = (eps_dot / n) / (n <sigma_T v> * <k0>)

    For massless emission: r_diss - 1 ~ const (recovers Schmidt et al. ansatz).
    For massive emission: r_diss(v) - 1 acquires velocity dependence via
        the suppression factor, breaking the (lambda, mu) rescaling symmetry.

    ====================================================================
    Legacy parametrization of the amplitude C0
    ====================================================================
    This function is retained as an explicit phenomenological compatibility
    closure.  The preferred quantitative path now evaluates the long-range
    identical-fermion differential emission kernel in emission_kernel.py.

    LB2026 Eq. (5.14a) gives the volumetric energy-loss rate for chi-V
    (long-range mediation, massless emission):
        eps_dot = [(44 - 32*pi^2) / (m_chi^{5/2} * pi^{7/2})] *
                  n^2 * T^{3/2} * (5/7) * g^6

    with T the temperature (= 1D velocity dispersion in natural units), n the
    number density, and g = sqrt(4*pi*alpha_D) the dark fine-structure
    coupling.

    Following the derivation:
        r_diss - 1 = eps_dot / [n^2 <sigma_T v> <k0>]
    with <sigma_T v> the thermally-averaged transfer cross section and
    <k0> the mean radiated energy per collision. For massless emission,
    <k0> ~ c_k * mu * v^2 with c_k ~ 0.4 (LB2026 Sec 3.5.2 / Sec 4),
    and the velocity dependence cancels, yielding a constant r_diss - 1.

    Substituting Eq. (2.6) for sigma_T, Eq. (5.14a) for eps_dot, and
    <k0> = c_k * mu * v^2, after angular and MB averaging:

        r_diss - 1 (massless) = N_quad * g^6 * T^{3/2} / m_chi^{5/2}
                                 / (sigma_T_v * c_k * mu * v^2)

    where N_quad = |44 - 32*pi^2| * (5/7) / pi^{7/2} is the LB2026
    prefactor for chi-V (and analogous factors for other channels; see
    LB2026 Table 2 for the ratios). For the chi-phi channel, the
    prefactor differs by the ratio 144/(...) per LB2026 Eq. (5.14b).

    The present project uses the explicit fiducial value C0 = 0.05 for all
    benchmarks.  It is a provisional closure, not a calibration to an
    elastic core-collapse time and not a validated prediction of the
    channel-dependent Lagrangian parameters.  C0_from_LB2026() is retained
    only as a diagnostic; it is not used by the preferred M1/M2 cooling path.

    ====================================================================
    Massive-emission suppression (LB2026 Sec 5.2.2)
    ====================================================================
    For massive emission, the current implementation uses a Boltzmann
    threshold ansatz.  The full phase-space integral of Eqs. (3.26)-(3.28)
    is evaluated directly by emission_kernel.py for the supported channels;
    this threshold ansatz remains a diagnostic for comparison only.

    Model dependence enters via the C_V (vector, Eq 3.27) vs C_phi
    (scalar, Eq 3.26) correction factors, which differ in shape:
        C_V   = 1 + m_V^2 / (2 omega^2)      (vector: enhanced at low omega)
        C_phi = 1 - m_phi^2 / omega^2         (scalar: suppressed at low omega,
                                               vanishes at threshold)
    These modify the pre-exponential factor and introduce a model-dependent
    power-law prefactor on top of the Boltzmann suppression.

    Returns r_diss (dimensionless).
    """
    v = np.asarray(v_km_s, dtype=float)
    v_ref = 100.0  # km/s calibration point
    # C0 is the asymptotic (massless, high-v) value of r_diss - 1.
    # It is currently a provisional phenomenological closure shared by the
    # benchmarks. The channel-dependent value must come from the full
    # differential-emission calculation before quantitative claims are made.
    C0 = DEFAULT_C0

    if p.emission_type == 'massless' or p.m_mediator == 0:
        # Massless: velocity-independent (recovers Schmidt et al. constant r_diss)
        # This is the control case; rescaling symmetry holds.
        return np.full_like(v, 1.0 + C0)

    # Massive emission: compute T = mu v^2 / 2 (CM kinetic energy in GeV)
    T_at_v = 0.5 * p.mu * (v * KM_S_TO_NAT)**2  # GeV
    m = p.m_mediator  # GeV

    # Boltzmann suppression (validated by LB2026 Sec 5.2.2)
    # Clip the exponent to avoid overflow in exp
    boltzmann = np.exp(-np.clip(m / T_at_v, 0, 700))

    # Model-dependent pre-exponential factor from C_V / C_phi (LB2026 Eqs 3.26, 3.27)
    # For a typical emitted energy omega ~ T (CM kinetic energy per particle),
    # the correction factors modify the effective amplitude:
    #   Vector (C_V):   enhancement at low omega (1 + m_V^2/(2 omega^2))
    #   Scalar (C_phi): suppression at low omega (1 - m_phi^2/omega^2), vanishes at threshold
    # We evaluate at omega ~ T (the characteristic CM energy):
    omega = np.clip(T_at_v, 1e-30, None)
    if 'V' in p.model:
        # Vector emission: C_V = 1 + m_V^2/(2 omega^2)
        # This ENHANCES emission near threshold (longitudinal mode contribution)
        model_factor = 1.0 + m**2 / (2.0 * omega**2)
    else:
        # Scalar emission: C_phi = 1 - m_phi^2/omega^2
        # This SUPPRESSES emission near threshold (phase space closure)
        # and vanishes at threshold (omega = m), as required by kinematics
        model_factor = np.clip(1.0 - m**2 / omega**2, 0.0, None)

    # r_diss - 1 = C0 * model_factor * boltzmann
    # At v >> v*: boltzmann -> 1, model_factor -> 1, recovers massless limit
    # At v << v*: boltzmann -> 0, r_diss -> 1 (emission suppressed)
    rdiss_minus_1 = C0 * model_factor * boltzmann
    return 1.0 + rdiss_minus_1


def C0_from_LB2026(p: DSIDMParameters, c_k: float = 0.4) -> float:
    """Estimate the asymptotic (massless) r_diss - 1 from LB2026 Eq. (5.14).

    For the chi-V channel (long-range mediation, identical fermions),
    LB2026 Eq. (5.14a) gives the volumetric energy-loss rate as:

        eps_dot = (44 - 32*pi^2) * n^2 * T^{3/2} / (m_chi^{5/2} * pi^{7/2})
                  * (5/7) * g^6

    (rendering of the published formula; the PDF extraction loses some
    structure, see LB2026 Sec. 5.2.1 for the exact expression).

    And r_diss - 1 = eps_dot / (n^2 * <sigma_T v> * <k0>),
    with <k0> = c_k * mu * v^2 the mean radiated energy per collision
    (LB2026 Sec 3.5.2; c_k ~ 0.4) and <sigma_T v> the thermally-averaged
    transfer cross section.

    IMPORTANT CAVEAT:
        (1) The PDF extraction of Eq. (5.14) loses some subscripts and
            prefactors. A direct dimensional substitution gives C0 of the
            wrong order of magnitude (~1e-16 rather than ~0.05), suggesting
            that either the mass scale in the denominator is the LIGHT
            mediator mass (not m_chi) or that a factor of m_chi has been
            absorbed into the coupling definition. A definitive check
            requires comparing against the numerical results in LB2026
            Fig. 9, not the analytical formula alone.
        (2) For other channels (chi-phi, S-V, etc.), the prefactor differs
            by the ratios given in LB2026 Table 2.
        (3) The function returns the order-of-magnitude estimate; the
            The ABSOLUTE value used in r_diss() is the provisional
            DEFAULT_C0 = 0.05 and is not derived from first principles here.

    Args:
        p: DSIDMParameters instance.
        c_k: mean radiated energy coefficient (~0.4 default, LB2026 Sec 3.5.2).

    Returns:
        C0 = r_diss - 1 in the massless limit at v = 100 km/s.
        ORDER-OF-MAGNITUDE diagnostic only; it is not used to set DEFAULT_C0.
    """
    g = p.g_eff
    mu = p.mu
    m_chi = p.m_chi
    # LB2026 Eq. (5.14a) prefactor for chi-V, long-range, massless.
    N_quad = abs(44.0 - 32.0 * np.pi**2) * (5.0 / 7.0) / np.pi**3.5

    vref = 100.0  # km/s
    vref_nat = vref * KM_S_TO_NAT
    T_at_vref = 0.5 * mu * vref_nat**2  # GeV

    sigma_T_vref = sigma_T_born(np.array([vref]), p)[0]  # cm^2/g
    sigma_T_nat = sigma_T_vref / GEV2_TO_CM2G  # GeV^-2
    sigma_v_nat = sigma_T_nat * vref_nat  # GeV^-1
    k0 = c_k * mu * vref_nat**2  # GeV

    # r_diss - 1 = N_quad * g^6 * T^{3/2} / m_chi^{5/2} / (sigma_v_nat * k0)
    C0 = N_quad * g**6 * T_at_vref**1.5 / m_chi**2.5 / (sigma_v_nat * k0)
    return float(C0)


# ============================================================
# Convenience: parameter sets for the three benchmark scenarios
# ============================================================

def benchmark_models() -> dict:
    """Return three benchmark dSIDM parameter sets for P1 comparison.

    Calibration:
        Mediator mass chosen so v* = sqrt(2 m_mediator / mu) falls in the
        astrophysically relevant range (100-1000 km/s).
        For m_chi = 10 GeV, mu = 5 GeV:
            v* = 200 km/s  ->  m = 1.1 keV = 1.1e-6 GeV
            v* = 500 km/s  ->  m = 6.95 keV = 6.95e-6 GeV

    M1: Dark photon bremsstrahlung, m_V = 1.1 keV (v*~200 km/s)
        - Emission suppressed below v*: symmetry broken at dwarf scale
        - Efficient above v*: MW and cluster scales
    M2: Scalar emission, m_phi = 6.95 keV (v*~500 km/s)
        - Different suppression profile (C_phi vs C_V)
        - Symmetry broken only above MW scale
    M3: Massless vector emission (control)
        - Constant r_diss; rescaling symmetry holds (Schmidt et al. limit)
    """
    return {
        'M1_dark_photon_massive': DSIDMParameters(
            # m_V = 1.1 keV, v* ~ 200 km/s
            model='chi-V', m_chi=10.0, m_mediator=1.1e-6,
            alpha_D=1e-4, mediation='long', emission_type='massive'
        ),
        'M2_scalar_phi_massive': DSIDMParameters(
            # m_phi = 6.95 keV, v* ~ 500 km/s
            model='chi-phi', m_chi=10.0, m_mediator=6.95e-6,
            alpha_D=1e-4, mediation='long', emission_type='massive'
        ),
        'M3_massless_control': DSIDMParameters(
            # Use a very light but finite elastic mediator (10 eV) to regulate
            # the Rutherford divergence. The emitted vector is treated as
            # massless, so no emission threshold is assigned to this control.
            model='chi-V', m_chi=10.0, m_mediator=1e-8,
            alpha_D=1e-4, mediation='long', emission_type='massless'
        ),
    }


if __name__ == '__main__':
    # Quick sanity check
    models = benchmark_models()
    v = np.logspace(0, 3.5, 100)  # 10 to ~3000 km/s

    for name, p in models.items():
        sig = sigma_T_born(v, p)
        rd = r_diss(v, p)
        print(f"\n{name}:")
        print(f"  v=100 km/s: sigma_T/m = {sigma_T_born(np.array([100.0]), p)[0]:.3e} cm^2/g, "
              f"r_diss = {r_diss(np.array([100.0]), p)[0]:.4f}")
        print(f"  v=1000 km/s: sigma_T/m = {sigma_T_born(np.array([1000.0]), p)[0]:.3e} cm^2/g, "
              f"r_diss = {r_diss(np.array([1000.0]), p)[0]:.4f}")
