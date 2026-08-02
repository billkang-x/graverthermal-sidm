"""
Dissipative gravothermal halo evolution.

This module subclasses the GravothermalSIDM `Halo` class to add:
  1. Velocity-dependent elastic cross section: sigma_m(T) where T = v^2
     (1D velocity dispersion in code units). Implemented via thermal averaging
     of the particle-physics sigma_T(v) over a Maxwell-Boltzmann distribution.
  2. Dissipative cooling rate per unit mass, using the Schmidt et al.
     fluid closure:
        C_cool = (8/sqrt(pi)) (sigma_T/m_chi) rho nu^3 (r_diss - 1).
     This enters the energy equation as an extra delta_uc:
        delta_uc -= delta_t * C_cool
     (per Schmidt et al. 2026 Eq. 19, Appendix D).

The cooling rate has the parametric form (in code units).  For velocity-
dependent models, elastic transport still uses collision-weighted averages.
The preferred cooling path supplies an energy-weighted effective cross section
computed from the differential emission kernel; the older collision-weighted
``r_diss`` path remains only for compatibility.

For velocity-dependent sigma_m(v) and r_diss(v), the elastic transport uses:
    sigma_m  -> <sigma_T v>/m_chi / <v>   (effective cross section per mass)
    r_diss   -> <r_diss * sigma_T v> / <sigma_T v>   (effective dissipation)

The cooling term may additionally receive an energy-weighted effective cross
section from ``thermal_avg.effective_cooling_sigma_m``.  If it is not supplied,
the collision-weighted ``r_diss`` closure remains available as a compatibility
approximation.

USAGE:
    from dissipative_halo import DissipativeHalo
    halo = DissipativeHalo(record,
                           sigma_m_eff_callable=...,
                           rdiss_eff_callable=...,
                           # plus standard Halo kwargs: r_s, rho_s, etc.
                           )
    halo.evolve_halo(t_end=...)
"""

from __future__ import annotations

import os, sys
import numpy as np
from functools import partial

# Import the upstream Halo class
_HERE = os.path.dirname(os.path.abspath(__file__))
_EXT = os.path.normpath(os.path.join(_HERE, '..', '..', 'external', 'gravothermalsidm'))
if _EXT not in sys.path:
    sys.path.insert(0, _EXT)

from SourcePy.evolve import Halo as _UpstreamHalo, TDMA_solver
from cooling import specific_cooling_rate


class DissipativeHalo(_UpstreamHalo):
    """Gravothermal halo with velocity-dependent sigma_T(v) and r_diss(v).

    Two extensions over the upstream Halo:
      (a) sigma_m is no longer a constant; it is the thermal-averaged
          <sigma_T v>/m_chi / <v> evaluated at the local temperature T = v^2.
      (b) An extra cooling term is added to the energy equation in
          conduct_heat() representing dissipative energy loss.

    The user supplies callables:
        sigma_m_eff(T_km2_s2) -> effective sigma_m in cm^2/g
        rdiss_eff(T_km2_s2)   -> effective r_diss (dimensionless)
        cooling_sigma_m_eff(T_km2_s2) -> energy-weighted cooling sigma_m in cm^2/g
    where T_km2_s2 is the 1D velocity dispersion squared in (km/s)^2. These
    are built by thermal_avg.py and converted to code units internally.

    Additional kwargs:
        flag_dissipation: bool, default True
            If True, apply the cooling term. If False, only velocity-dependent
            elastic scattering is used (useful for elastic-only control).
        dissipation_prefactor: float, default 1.0
            Explicit sensitivity factor multiplying the fluid closure. It is
            not calibrated to an elastic core-collapse time.
    """

    def __init__(self, record, *, sigma_m_eff_callable=None,
                 rdiss_eff_callable=None, cooling_sigma_m_eff_callable=None,
                 flag_dissipation=True,
                 dissipation_prefactor=1.0, **kwargs):
        # Save callables BEFORE super().__init__ because update_derived_parameters
        # is called there and we override it to use these.
        self._sigma_m_eff_fn = sigma_m_eff_callable
        self._rdiss_eff_fn = rdiss_eff_callable
        self._cooling_sigma_m_eff_fn = cooling_sigma_m_eff_callable
        self.flag_dissipation = flag_dissipation
        self.dissipation_prefactor = dissipation_prefactor

        # cooling luminosity array (energy loss per unit time at each shell boundary)
        # We will allocate after n_shells is known.
        self.L_cool = None

        super().__init__(record, **kwargs)

    # ------------------------------------------------------------------
    # Override: initialize_elastic_scattering returns the F(v) for the
    # *shape* of the cross section. Since we handle sigma_m velocity
    # dependence directly via the thermal-average callable, we set F=1
    # (constant shape) and let update_derived_parameters apply the T-dependent
    # sigma_m.
    # ------------------------------------------------------------------
    def initialize_elastic_scattering(self, model_name):
        # We override to always return constant shape=1; the velocity dependence
        # is injected via sigma_m(T) in update_derived_parameters.
        # BUT: if the user passed 'constant' (default), we keep the upstream behavior.
        if self._sigma_m_eff_fn is None:
            return super().initialize_elastic_scattering(model_name)
        # Velocity-dependent case: F(v) = 1, sigma_m itself carries the dependence.
        return lambda x: 1.0

    # ------------------------------------------------------------------
    # Override: update_derived_parameters
    # Replace constant self.sigma_m with a T-dependent array when a callable
    # is provided. Also compute the cooling luminosity L_cool.
    # ------------------------------------------------------------------
    def update_derived_parameters(self):
        # specific energy and 1D velocity dispersion
        self.u = (3. / 2.) * self.p / self.rho
        self.v = np.sqrt(self.p / self.rho)

        # ---- velocity-dependent sigma_m ----
        if self._sigma_m_eff_fn is not None:
            # v is in dimensionless code units. Convert to dimensionful km/s
            # via scale_v, then take T = v_1d^2 in (km/s)^2 for the callable.
            # BUT: the callable was built with T in (km/s)^2 because sigma_T(v)
            # takes v in km/s. So we need v in km/s:
            v_km_s = self.v * (self.scale_v.to('km/s').value
                               if hasattr(self.scale_v, 'to') else float(self.scale_v.value))
            # T = v_1d^2  in (km/s)^2
            T_km2_s2 = v_km_s ** 2
            # Evaluate sigma_m_eff in cm^2/g, convert to dimensionless code units
            from astropy import units as ut
            sigma_m_with_units = self._sigma_m_eff_fn(T_km2_s2) * ut.cm**2 / ut.g
            self.sigma_m_arr = (sigma_m_with_units).to_value(self.scale_sigma_m)
        else:
            # constant sigma_m (upstream behavior)
            self.sigma_m_arr = np.full(self.n_shells, self.sigma_m)

        # ---- effective thermal conductivity (LMFP + SMFP interpolated) ----
        # Kinv_smfp = sigma_m(v) * F_smfp(v/w) / (b v)
        # Kinv_lmfp = 1 / (a C v p sigma_m(v) F_lmfp(v/w))
        # With F=1 (we set it above), and sigma_m -> sigma_m_arr:
        v_safe = np.clip(self.v, 1e-30, None)
        self.Kinv_smfp = self.sigma_m_arr * self.F_elastic_smfp(self.v / self.w) / (self.b * v_safe)
        self.Kinv_lmfp = 1.0 / (self.a * self.C * v_safe * self.p * self.sigma_m_arr *
                                self.F_elastic_lmfp(self.v / self.w))
        Keff = 1.0 / (self.Kinv_smfp + self.Kinv_lmfp)

        # luminosity (heat flux * 4 pi r^2)
        self.L[1:-1] = -self.r[1:-1] * self.r[1:-1] * (Keff[1:-1] + Keff[2:]) / 2. * \
                       (self.u[2:] - self.u[1:-1]) / ((self.r[2:] - self.r[:-2]) / 2.)
        self.L[0] = -self.r[0] * self.r[0] * (Keff[0] + Keff[1]) / 2. * \
                    (self.u[1] - self.u[0]) / (self.r[1] / 2.)
        self.L[-1] = 0.0

        # ---- dissipative cooling rate ----
        # The executable definition lives in cooling.specific_cooling_rate.
        # The historical derivation below is retained only for provenance;
        # the executable contract is the tested function in cooling.py.
        # Per Schmidt et al. 2026 Eq. 19, the volumetric cooling rate is
        #   C_vol(ρ, ν) = (8√π/3) (σ_T/m)(r_diss - 1) ρ² ν³
        # where ν is the 1D velocity dispersion (our self.v).
        #
        # We need the cooling rate PER UNIT MASS (energy/mass/time), obtained
        # by dividing the volumetric rate by ρ:
        #   C_cool = (8√π/3) (σ_T/m)(r_diss - 1) ρ ν³
        #
        # Equivalent derivation from microphysics:
        #   - Scattering rate per particle: Γ = (ρ/mχ)(σ/m)ν
        #   - Mean radiated energy per collision: <k0> = c0 × ½μν²  (c0 ≈ 0.4)
        #   - For identical particles μ = mχ/2, so <k0> = c0 × (mχ/4)ν²
        #   - Energy loss per unit mass = Γ × <k0> / mχ
        #                               = (ρ/mχ)(σ/m)ν × c0(mχ/4)ν² / mχ
        #                               = (c0/4)(σ/m)ρ ν³
        # Matching the two: (8√π/3) ≈ 4.73, while c0/4 = 0.1.
        # These differ because the scattering rate <σv> for a Yukawa/Born
        # cross section is NOT simply σ_T × v — the MB average introduces
        # a factor of √(8/π) and the Yukawa shape gives an additional O(1)
        # correction. The full relation is:
        #   C_cool = [8√π/3] × (σ/m)(r_diss-1) ρ ν³ × [1/⟨σ_T v⟩_MB_norm]
        # We absorb the O(1) normalization into the effective sigma_m_arr
        # (which is already the MB-averaged ⟨σ_T v⟩/⟨v⟩).
        # The dimensionless prefactor in code units is:
        #   prefactor = (8√π/3) / √(8/π) = (8√π/3) × √(π/8) = (8π/3)/√8 = π√2/3 ≈ 1.481
        # But this assumes ⟨v⟩ = √(8ν²/π) (3D mean speed), and our sigma_m_arr
        # already includes this averaging. So the consistent prefactor is:
        #   prefactor = (8√π/3) × √(π/8) / √(8/π) = (8√π/3) × (π/8) = π√π/3 ≈ 1.860
        # For practical purposes, we use the microphysics derivation:
        #   C_cool = (c0/4) × (σ/m)(r_diss-1) ρ ν³  with c0=0.4 → prefactor = 0.1
        # The current implementation uses the explicit source coefficient in
        # cooling.py; dissipation_prefactor is an optional user-supplied scale,
        # not a calibration to an elastic core-collapse time.
        if self.flag_dissipation and (
            self._cooling_sigma_m_eff_fn is not None or self._rdiss_eff_fn is not None
        ):
            v_km_s = self.v * (self.scale_v.to('km/s').value
                               if hasattr(self.scale_v, 'to') else float(self.scale_v.value))
            T_km2_s2 = v_km_s ** 2
            if self._cooling_sigma_m_eff_fn is not None:
                from astropy import units as ut
                from cooling import specific_cooling_rate_from_moment
                cooling_sigma = self._cooling_sigma_m_eff_fn(T_km2_s2)
                cooling_sigma_units = cooling_sigma * ut.cm**2 / ut.g
                cooling_sigma_code = cooling_sigma_units.to_value(self.scale_sigma_m)
                self.cooling_sigma_m_arr = cooling_sigma_code
                self.C_cool = specific_cooling_rate_from_moment(
                    cooling_sigma_code,
                    self.rho,
                    v_safe,
                    self.dissipation_prefactor,
                )
            else:
                # Compatibility path for callers that have not yet supplied
                # the exact energy-weighted cooling moment.
                rdiss_arr = self._rdiss_eff_fn(T_km2_s2)
                self.C_cool = specific_cooling_rate(
                    self.sigma_m_arr,
                    self.rho,
                    v_safe,
                    rdiss_arr,
                    self.dissipation_prefactor,
                )
            # L_cool is the cooling "luminosity" at each shell boundary,
            # analogously to L (which is heat luminosity). For the energy
            # update we use C_cool directly (per-shell rate).
            if self.L_cool is None:
                self.L_cool = np.empty(self.n_shells, dtype=np.float64)
            # C_cool is already an energy loss per unit mass per unit time.
            self.L_cool[:] = self.C_cool
        else:
            if self.L_cool is None:
                self.L_cool = np.zeros(self.n_shells, dtype=np.float64)
            else:
                self.L_cool[:] = 0.0

        # Knudsen number (use local sigma_m_arr)
        self.Kn = 1.0 / (np.sqrt(np.clip(self.p, 1e-30, None)) * self.sigma_m_arr)

    # ------------------------------------------------------------------
    # Override: conduct_heat to add the cooling term
    # ------------------------------------------------------------------
    def conduct_heat(self):
        # update conduction counter
        self.n_conduction += 1

        # determine minimum time step
        delta_t = self.get_timestep()

        # heat conduction contribution (upstream)
        self.delta_uc[1:] = -delta_t * ((self.L[1:] - self.L[:-1]) /
                                        (self.m[1:] - self.m[:-1]))
        self.delta_uc[0] = -delta_t * (self.L[0] / self.m[0])

        # dissipative cooling contribution
        # C_cool is a specific energy loss rate, so it is subtracted directly.
        if self.flag_dissipation and (
            self._cooling_sigma_m_eff_fn is not None or self._rdiss_eff_fn is not None
        ):
            self.delta_uc -= delta_t * self.C_cool

        # pressure change (adiabatic: delta_p/p = delta_u/u)
        delta_pc = self.p * self.delta_uc / np.clip(self.u, 1e-30, None)

        # update variables
        self.p += delta_pc
        self.u += self.delta_uc

        # update time
        self.t_before = self.t
        self.t += delta_t

    # ------------------------------------------------------------------
    # Override: get_timestep to include a cooling time criterion
    # ------------------------------------------------------------------
    def get_timestep(self):
        # upstream criteria
        delta_t1 = min(abs(self.u[0] / (self.L[0] / self.m[0])),
                       min(abs(self.u[1:] /
                               ((self.L[1:] - self.L[:-1]) /
                                (self.m[1:] - self.m[:-1])))))
        delta_t2 = min(1.0 / (self.rho * self.v))

        # cooling time criterion: delta_t * C_cool < epsilon * u
        if self.flag_dissipation and (
            self._cooling_sigma_m_eff_fn is not None or self._rdiss_eff_fn is not None
        ) and hasattr(self, 'C_cool'):
            u_safe = np.clip(self.u, 1e-30, None)
            cool_rates = self.C_cool  # rate = (energy/mass)/time
            delta_t3_arr = u_safe / np.clip(cool_rates, 1e-30, None)
            # Only consider shells where cooling is significant
            mask = cool_rates > 1e-30 * u_safe
            if np.any(mask):
                delta_t3 = np.min(delta_t3_arr[mask])
            else:
                delta_t3 = np.inf
        else:
            delta_t3 = np.inf

        # determine minimum time step
        if self.flag_timestep_use_relaxation and self.flag_timestep_use_energy:
            delta_t = min(delta_t1, delta_t2, delta_t3)
        elif self.flag_timestep_use_relaxation:
            delta_t = min(delta_t2, delta_t3)
        elif self.flag_timestep_use_energy:
            delta_t = min(delta_t1, delta_t3)
        else:
            delta_t = delta_t3

        return self.t_epsilon * delta_t
