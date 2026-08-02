"""
Maxwell-Boltzmann thermal averaging for velocity-dependent cross sections.

The gravothermal fluid equations use local (thermally-averaged) scattering rates.
For velocity-dependent sigma_T(v) and r_diss(v), the relevant averages are:

    <sigma_T v>      = int dv  f(v; T) * sigma_T(v) * v
    <r_diss sigma_T v> = int dv  f(v; T) * r_diss(v) * sigma_T(v) * v

where f(v; T) is the relative-velocity Maxwell-Boltzmann distribution.

Convention in the gravothermal code:
    T = nu^2  (one-dimensional velocity variance)
    For two identical particles, the relative velocity has component variance 2*T.

We use the standard SIDM MB relative-speed distribution
    f(v|T) = v^2 / (2*sqrt(pi)*T^(3/2)) * exp(-v^2 / (4T))
normalized so int_0^inf f(v) dv = 1.

NOTE on units: The thermal averages should be expressed in the SAME units as
the constant sigma_m used by GravothermalSIDM. That is:
    sigma_m = sigma_T / m_chi    [cm^2 / g]
    v       = relative speed     [in code units = v_dimensionless * w_units]

To keep things simple we work in DIMENSIONLESS form: the user passes
sigma_T(v)/m_chi as a callable returning cm^2/g for v in km/s, and we evaluate
the integral in km/s. The final <sigma_T v> / m_chi is then converted to
code units by the caller (Halo wrapper) using scale_sigma_m and scale_v.

Usage:
    from thermal_avg import make_thermal_interpolators
    sig_avg, rsig_avg = make_thermal_interpolators(sigma_v_callable, rdiss_v_callable,
                                                    v_axis_km_s, T_axis_km2_s2)
    # sig_avg(T)  -> <sigma_T v>/m_chi at temperature T (km^2/s * cm^2/g)
    # rsig_avg(T) -> <r_diss * sigma_T v>/m_chi at temperature T
    # r_diss_eff(T) = rsig_avg(T) / sig_avg(T)
"""

from __future__ import annotations

import numpy as np

from cooling import COOLING_PREFACTOR

try:
    from scipy.integrate import quad as _scipy_quad
except ImportError:  # The core regression tests only require NumPy.
    _scipy_quad = None


def _integrate_1d(integrand, lower: float, upper: float) -> float:
    """Integrate a smooth one-dimensional thermal kernel."""
    if _scipy_quad is not None:
        def scalar_integrand(value):
            result = np.asarray(integrand(value), dtype=float)
            if result.size != 1:
                raise ValueError("quadrature integrand must return one value for scalar input")
            return float(result.reshape(-1)[0])

        value, _ = _scipy_quad(
            scalar_integrand, lower, upper, limit=200, epsrel=1e-7
        )
        return float(value)

    grid = np.linspace(lower, upper, 16_385)
    try:
        values = np.asarray(integrand(grid), dtype=float)
        if values.shape != grid.shape:
            raise ValueError
    except (TypeError, ValueError):
        values = np.array([integrand(float(v)) for v in grid], dtype=float)
    integrate = getattr(np, "trapezoid", None)
    if integrate is None:
        integrate = np.trapz
    return float(integrate(values, grid))


def _linear_interp_with_extrapolation(x, y, x_new):
    """One-dimensional linear interpolation matching SciPy extrapolation."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x_new = np.asarray(x_new, dtype=float)
    out = np.interp(x_new, x, y)
    if x.size < 2:
        return np.full_like(x_new, y[0])
    left = x_new < x[0]
    right = x_new > x[-1]
    out[left] = y[0] + (x_new[left] - x[0]) * (y[1] - y[0]) / (x[1] - x[0])
    out[right] = y[-1] + (x_new[right] - x[-1]) * (y[-1] - y[-2]) / (x[-1] - x[-2])
    return out


# ----------------------------------------------------------------------
# Maxwell-Boltzmann relative-speed distribution (3D, identical particles)
# ----------------------------------------------------------------------
def mb_pdf(v: float | np.ndarray, T: float) -> np.ndarray:
    """3D relative-speed Maxwell-Boltzmann PDF.

    f(v) = v^2 / (2*sqrt(pi)*T^{3/2}) * exp(-v^2 / (4T))

    Normalized so int_0^inf f(v) dv = 1.
    v in same units as sqrt(T).
    """
    v = np.asarray(v, dtype=float)
    if T <= 0:
        out = np.zeros_like(v)
        return out
    # Relative speed of two identical Maxwellian populations: each relative
    # component has variance 2*T when T is the one-particle 1D variance.
    return v**2 / (2.0 * np.sqrt(np.pi) * T**1.5) * np.exp(-v**2 / (4.0 * T))


def mean_relative_speed(T: float | np.ndarray) -> np.ndarray:
    """Mean relative speed for the Maxwell PDF used above."""
    return 4.0 * np.sqrt(np.asarray(T, dtype=float) / np.pi)


# ----------------------------------------------------------------------
# Single-point thermal average
# ----------------------------------------------------------------------
def thermal_avg_sigma_v(
    sigma_v_callable,
    T: float,
    v_max_factor: float = 12.0,
) -> float:
    """<sigma_T v> at temperature T (1D velocity dispersion squared).

    sigma_v_callable(v) -> sigma_T/m_chi in cm^2/g, v in km/s
    T in (km/s)^2  (i.e. T = v_1d^2)
    Returns  <sigma_T v>/m_chi in (cm^2/g)*(km/s).
    """
    if T <= 0:
        return 0.0
    sigma_max = np.sqrt(T)
    v_max = max(v_max_factor * sigma_max, 1.0)

    integrand = lambda v: sigma_v_callable(v) * v * mb_pdf(v, T)
    return _integrate_1d(integrand, 1e-3, v_max)


def thermal_avg_rsig_v(
    sigma_v_callable,
    rdiss_v_callable,
    T: float,
    v_max_factor: float = 12.0,
) -> float:
    """<r_diss * sigma_T v> at temperature T."""
    if T <= 0:
        return 0.0
    sigma_max = np.sqrt(T)
    v_max = max(v_max_factor * sigma_max, 1.0)

    integrand = lambda v: rdiss_v_callable(v) * sigma_v_callable(v) * v * mb_pdf(v, T)
    return _integrate_1d(integrand, 1e-3, v_max)


def thermal_avg_dissipation_moment(
    sigma_v_callable,
    rdiss_v_callable,
    T: float,
    v_max_factor: float = 12.0,
) -> float:
    """Return the exact energy-weighted cooling moment for one component.

    For identical particles the self-scattering cooling closure is

        C_vol / rho = rho * (1/4) * <sigma_m(v) [r(v)-1] v^3>,

    where the average uses the relative-speed Maxwell PDF.  The factor 1/4
    is the pair-counting/identical-particle factor in Schmidt et al. App. A.
    The returned value has units ``(cm^2/g) * (km/s)^3``.
    """
    if T <= 0:
        return 0.0
    sigma_max = np.sqrt(T)
    v_max = max(v_max_factor * sigma_max, 1.0)

    integrand = lambda v: (
        sigma_v_callable(v)
        * (rdiss_v_callable(v) - 1.0)
        * v**3
        * mb_pdf(v, T)
    )
    return 0.25 * _integrate_1d(integrand, 1e-3, v_max)


# ----------------------------------------------------------------------
# Tabulated thermal averages over a (v, T) grid -> interpolators
# ----------------------------------------------------------------------
def make_thermal_interpolators(
    sigma_v_callable,
    rdiss_v_callable,
    T_axis_km2_s2: np.ndarray,
    n_v_quad: int = 320,
    v_max_factor: float = 12.0,
):
    """Precompute <sigma v> and <r*sigma v> on a temperature grid.

    Returns two callables sig_avg(T), rsig_avg(T) that can be evaluated
    at arbitrary T (linear-interpolated in log T).
    """
    T_axis = np.asarray(T_axis_km2_s2, dtype=float)
    sig_arr = np.empty_like(T_axis)
    rsig_arr = np.empty_like(T_axis)

    for i, T in enumerate(T_axis):
        sig_arr[i] = thermal_avg_sigma_v(sigma_v_callable, T, v_max_factor)
        rsig_arr[i] = thermal_avg_rsig_v(sigma_v_callable, rdiss_v_callable, T, v_max_factor)

    # Build log-log interpolators (T > 0)
    logT = np.log10(np.clip(T_axis, 1e-6, None))
    logsig = np.log10(np.clip(sig_arr, 1e-300, None))
    logrsig = np.log10(np.clip(rsig_arr, 1e-300, None))

    def sig_avg(T):
        T = np.atleast_1d(np.asarray(T, dtype=float))
        out = 10.0 ** _linear_interp_with_extrapolation(
            logT, logsig, np.log10(np.clip(T, 1e-6, None))
        )
        return out

    def rsig_avg(T):
        T = np.atleast_1d(np.asarray(T, dtype=float))
        out = 10.0 ** _linear_interp_with_extrapolation(
            logT, logrsig, np.log10(np.clip(T, 1e-6, None))
        )
        return out

    return sig_avg, rsig_avg, sig_arr, rsig_arr


# ----------------------------------------------------------------------
# Convenience wrapper: effective sigma_m(T) and r_diss(T) for the Halo class
# ----------------------------------------------------------------------
def effective_sigma_m_and_rdiss(
    sigma_v_callable,
    rdiss_v_callable,
    T_axis_km2_s2: np.ndarray,
):
    """Return callables sigma_m_eff(T) and r_diss_eff(T) for the Halo wrapper.

    sigma_m_eff(T)  = <sigma_T v>/m_chi  / <v>   (effective cross-section per mass at T)
                    = <sigma_T v>/m_chi  / (4 sqrt(T/pi))

    r_diss_eff(T)    = <r_diss sigma_T v> / <sigma_T v>
    """
    sig_avg, rsig_avg, sig_arr, rsig_arr = make_thermal_interpolators(
        sigma_v_callable, rdiss_v_callable, T_axis_km2_s2
    )

    # Mean relative speed for two identical Maxwellian populations.
    vbar = mean_relative_speed(T_axis_km2_s2)

    sigma_m_arr = sig_arr / np.clip(vbar, 1e-10, None)
    rdiss_arr = rsig_arr / np.clip(sig_arr, 1e-300, None)

    logT = np.log10(np.clip(T_axis_km2_s2, 1e-6, None))
    logsig = np.log10(np.clip(sigma_m_arr, 1e-300, None))
    logrd = np.log10(np.clip(rdiss_arr, 1e-10, None))

    def sigma_m_eff(T):
        T = np.atleast_1d(np.asarray(T, dtype=float))
        return 10.0 ** _linear_interp_with_extrapolation(
            logT, logsig, np.log10(np.clip(T, 1e-6, None))
        )

    def rdiss_eff(T):
        T = np.atleast_1d(np.asarray(T, dtype=float))
        return 10.0 ** _linear_interp_with_extrapolation(
            logT, logrd, np.log10(np.clip(T, 1e-6, None))
        )

    return sigma_m_eff, rdiss_eff, sigma_m_arr, rdiss_arr


def effective_cooling_sigma_m(
    sigma_v_callable,
    rdiss_v_callable,
    T_axis_km2_s2: np.ndarray,
    v_max_factor: float = 12.0,
):
    """Return the energy-weighted effective cooling cross section.

    ``sigma_cool_eff(T)`` is defined so that the exact moment can be passed
    to the fluid closure as

        C_vol / rho = (8/sqrt(pi)) * rho * sigma_cool_eff(T) * T^(3/2).

    For constant ``sigma_T/m`` and constant ``r_diss``, this returns exactly
    ``sigma_T/m * (r_diss - 1)``.  Unlike ``r_diss_eff``, it retains the
    velocity weighting required by the radiated-energy rate.
    """
    T_axis = np.asarray(T_axis_km2_s2, dtype=float)
    if T_axis.ndim != 1 or T_axis.size < 2 or np.any(T_axis <= 0):
        raise ValueError("T_axis_km2_s2 must be a 1D array with at least two positive values")

    moment_arr = np.array([
        thermal_avg_dissipation_moment(
            sigma_v_callable, rdiss_v_callable, T, v_max_factor
        )
        for T in T_axis
    ])
    denominator = COOLING_PREFACTOR * T_axis**1.5
    cooling_arr = np.divide(
        moment_arr,
        denominator,
        out=np.zeros_like(moment_arr),
        where=denominator > 0,
    )
    logT = np.log10(T_axis)

    def sigma_cool_eff(T):
        T = np.atleast_1d(np.asarray(T, dtype=float))
        values = _linear_interp_with_extrapolation(
            logT, cooling_arr, np.log10(np.clip(T, T_axis[0], T_axis[-1]))
        )
        return np.clip(values, 0.0, None)

    return sigma_cool_eff, cooling_arr


def thermal_avg_cooling_kernel_moment(
    cooling_sigma_m_callable,
    T: float,
    v_max_factor: float = 12.0,
) -> float:
    """Return ``(1/4) <sigma_cool_m(v) v^3>`` for a direct microphysical kernel."""
    if T <= 0:
        return 0.0
    sigma_max = np.sqrt(T)
    v_max = max(v_max_factor * sigma_max, 1.0)
    integrand = lambda v: (
        cooling_sigma_m_callable(v) * v**3 * mb_pdf(v, T)
    )
    return 0.25 * _integrate_1d(integrand, 1e-3, v_max)


def effective_cooling_sigma_m_from_kernel(
    cooling_sigma_m_callable,
    T_axis_km2_s2: np.ndarray,
    v_max_factor: float = 12.0,
):
    """Build the exact cooling effective cross section from a direct kernel.

    The callable should return the per-speed combination
    ``2 Q(v)/(m_chi^2 v^2)`` in ``cm^2/g``, where ``Q`` is the emitted-energy
    weighted cross section.  This function performs the remaining relative-
    speed Maxwell average and normalizes it to the Schmidt et al. closure.
    """
    T_axis = np.asarray(T_axis_km2_s2, dtype=float)
    if T_axis.ndim != 1 or T_axis.size < 2 or np.any(T_axis <= 0):
        raise ValueError("T_axis_km2_s2 must be a 1D array with at least two positive values")
    moment_arr = np.array([
        thermal_avg_cooling_kernel_moment(
            cooling_sigma_m_callable, T, v_max_factor
        )
        for T in T_axis
    ])
    denominator = COOLING_PREFACTOR * T_axis**1.5
    cooling_arr = np.divide(
        moment_arr,
        denominator,
        out=np.zeros_like(moment_arr),
        where=denominator > 0,
    )
    logT = np.log10(T_axis)

    def sigma_cool_eff(T):
        T = np.atleast_1d(np.asarray(T, dtype=float))
        values = _linear_interp_with_extrapolation(
            logT, cooling_arr, np.log10(np.clip(T, T_axis[0], T_axis[-1]))
        )
        return np.clip(values, 0.0, None)

    return sigma_cool_eff, cooling_arr


if __name__ == '__main__':
    # Self-test: constant sigma and r_diss should give back the constants.
    sig0 = 5.0  # cm^2/g
    rd0 = 1.05

    sigma_fn = lambda v: np.full_like(np.asarray(v, dtype=float), sig0)
    rdiss_fn = lambda v: np.full_like(np.asarray(v, dtype=float), rd0)

    T_axis = np.logspace(1, 6, 30)  # (km/s)^2
    sig_m, rd_eff, sa, ra = effective_sigma_m_and_rdiss(sigma_fn, rdiss_fn, T_axis)

    print("Constant-sigma test (should recover sig0=5, rd0=1.05):")
    print(f"  sigma_m_eff at T=1e3: {sig_m(np.array([1e3]))[0]:.4f} (expect ~{sig0:.4f})")
    print(f"  sigma_m_eff at T=1e5: {sig_m(np.array([1e5]))[0]:.4f} (expect ~{sig0:.4f})")
    print(f"  rdiss_eff at T=1e3:   {rd_eff(np.array([1e3]))[0]:.4f} (expect ~{rd0:.4f})")
    print(f"  rdiss_eff at T=1e5:   {rd_eff(np.array([1e5]))[0]:.4f} (expect ~{rd0:.4f})")
