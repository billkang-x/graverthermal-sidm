"""Regression tests for the corrected fluid-physics closures."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src" / "fluid_runner"))

from cooling import (
    COOLING_PREFACTOR,
    specific_cooling_rate,
    specific_cooling_rate_from_moment,
)
sys.path.insert(0, str(ROOT / "src" / "cross_sections"))
from dsidm_models import (
    benchmark_models,
    born_expansion_parameter,
    is_born_valid,
    sigma_T_born,
)
from born_valid_scan import (
    parameters_at_eta,
    parameters_at_mass_threshold,
    parameters_for_target_sigma,
    parameters_for_target_sigma_at_threshold,
    threshold_velocity_km_s,
)
from emission_kernel import (
    microscopic_cooling_sigma_m,
    radiated_energy_cross_section,
)
from thermal_avg import (
    effective_cooling_sigma_m,
    mb_pdf,
    mean_relative_speed,
    thermal_avg_dissipation_moment,
    thermal_avg_sigma_v,
)
sys.path.insert(0, str(ROOT / "src" / "P3_rescaling"))
from rescale import (
    M_20PC,
    M_90PC,
    absolute_mass_fit,
    projected_enclosed_mass,
)
from snapshot_diagnostics import best_mass_fit, _log_interp_profile, diagnose_snapshot
sys.path.insert(0, str(ROOT / "src" / "P6_parameter_scan"))
from run_parameter_scan import _model_point, _parse_float_list
from run_direct_time_scan import _nearest_snapshot, _parse_times
from audit_physical_activity import parse_point_id
from audit_nfw_cosmology import isolated_nfw_extrapolation, nfw_mass_factor
from compare_direct_controls import compare_frames


class RelativeSpeedTests(unittest.TestCase):
    def test_pdf_normalization_and_mean(self):
        temperature = 125.0
        velocity = np.linspace(0.0, 16.0 * np.sqrt(temperature), 200_001)
        pdf = mb_pdf(velocity, temperature)
        integrate = getattr(np, "trapezoid", None)
        if integrate is None:
            integrate = np.trapz
        normalization = integrate(pdf, velocity)
        measured_mean = integrate(velocity * pdf, velocity)
        self.assertAlmostEqual(normalization, 1.0, places=7)
        self.assertAlmostEqual(
            measured_mean,
            float(mean_relative_speed(temperature)),
            places=6,
        )

    def test_constant_cross_section_is_preserved(self):
        sigma_m = 3.25
        temperature = 400.0
        sigma_v = thermal_avg_sigma_v(
            lambda velocity: np.full_like(np.asarray(velocity), sigma_m),
            temperature,
        )
        recovered = sigma_v / float(mean_relative_speed(temperature))
        self.assertAlmostEqual(recovered, sigma_m, places=5)

    def test_length_one_model_output_is_accepted(self):
        temperature = 225.0
        sigma_v = thermal_avg_sigma_v(
            lambda velocity: np.full(np.atleast_1d(velocity).shape, 2.5),
            temperature,
        )
        recovered = sigma_v / float(mean_relative_speed(temperature))
        self.assertAlmostEqual(recovered, 2.5, places=5)


class CoolingTests(unittest.TestCase):
    def test_source_coefficient_and_elastic_limit(self):
        self.assertAlmostEqual(COOLING_PREFACTOR, 8.0 / np.sqrt(np.pi))
        elastic = specific_cooling_rate(5.0, 2.0, 3.0, 1.0)
        self.assertEqual(float(elastic), 0.0)

    def test_specific_rate_has_nu_cubed_scaling(self):
        base = specific_cooling_rate(2.0, 3.0, 4.0, 1.05)
        doubled_nu = specific_cooling_rate(2.0, 3.0, 8.0, 1.05)
        self.assertAlmostEqual(float(doubled_nu / base), 8.0)

    def test_constant_model_fractional_cooling_is_rescaling_invariant(self):
        sigma_m, rho, nu, time, rdiss = 2.0, 3.0, 4.0, 5.0, 1.05
        u = 1.5 * nu**2
        invariant = specific_cooling_rate(sigma_m, rho, nu, rdiss) * time / u

        length_scale = 7.0
        mass_scale = 11.0
        sigma_scaled = sigma_m * length_scale**2 / mass_scale
        rho_scaled = rho * mass_scale / length_scale**3
        nu_scaled = nu * np.sqrt(mass_scale / length_scale)
        time_scaled = time * np.sqrt(length_scale**3 / mass_scale)
        u_scaled = 1.5 * nu_scaled**2
        scaled = (
            specific_cooling_rate(sigma_scaled, rho_scaled, nu_scaled, rdiss)
            * time_scaled
            / u_scaled
        )
        self.assertAlmostEqual(float(scaled / invariant), 1.0, places=12)

    def test_energy_weighted_moment_recovers_constant_closure(self):
        sigma_m, rdiss, temperature = 3.25, 1.05, 400.0
        sigma_fn = lambda velocity: np.full_like(np.asarray(velocity, dtype=float), sigma_m)
        rdiss_fn = lambda velocity: np.full_like(np.asarray(velocity, dtype=float), rdiss)
        moment = thermal_avg_dissipation_moment(sigma_fn, rdiss_fn, temperature)
        expected = COOLING_PREFACTOR * sigma_m * (rdiss - 1.0) * temperature**1.5
        self.assertAlmostEqual(moment, expected, places=5)

        cooling_eff, _ = effective_cooling_sigma_m(
            sigma_fn, rdiss_fn, np.logspace(1, 5, 24)
        )
        self.assertAlmostEqual(
            float(cooling_eff(np.array([temperature]))[0]),
            sigma_m * (rdiss - 1.0),
            places=5,
        )

    def test_specific_rate_from_moment_matches_source_closure(self):
        sigma_cool, rho, nu = 0.25, 2.0, 4.0
        from_moment = specific_cooling_rate_from_moment(sigma_cool, rho, nu)
        direct = specific_cooling_rate(sigma_cool, rho, nu, 2.0)
        self.assertAlmostEqual(float(from_moment), float(direct), places=12)


class BornMaskTests(unittest.TestCase):
    def test_fiducial_massive_points_are_flagged(self):
        models = benchmark_models()
        self.assertGreater(born_expansion_parameter(models["M1_dark_photon_massive"]), 1.0)
        self.assertGreater(born_expansion_parameter(models["M2_scalar_phi_massive"]), 1.0)
        self.assertFalse(is_born_valid(models["M1_dark_photon_massive"]))

    def test_lowering_alpha_can_enter_lenient_mask(self):
        model = benchmark_models()["M1_dark_photon_massive"]
        model.alpha_D *= 1.0e-4
        self.assertTrue(is_born_valid(model))

    def test_parameterization_hits_requested_eta(self):
        model = benchmark_models()["M2_scalar_phi_massive"]
        controlled = parameters_at_eta(model, 0.3)
        self.assertAlmostEqual(born_expansion_parameter(controlled), 0.3, places=12)
        self.assertTrue(is_born_valid(controlled))

    def test_mass_hierarchy_preserves_threshold_and_eta(self):
        model = benchmark_models()["M1_dark_photon_massive"]
        threshold = threshold_velocity_km_s(model)
        changed = parameters_at_mass_threshold(model, 1.5, threshold, 0.1)
        self.assertAlmostEqual(born_expansion_parameter(changed), 0.1, places=12)
        self.assertAlmostEqual(
            threshold_velocity_km_s(changed), threshold, places=10
        )

    def test_target_sigma_solver_preserves_born_control(self):
        model = benchmark_models()["M2_scalar_phi_massive"]
        solved = parameters_for_target_sigma(model, 0.1, 1.0, 100.0)
        self.assertAlmostEqual(born_expansion_parameter(solved), 0.1, places=12)
        sigma = sigma_T_born(np.array([100.0]), solved)[0]
        self.assertAlmostEqual(float(sigma), 1.0, places=10)

    def test_target_sigma_solver_accepts_requested_threshold(self):
        model = benchmark_models()["M1_dark_photon_massive"]
        solved = parameters_for_target_sigma_at_threshold(model, 0.1, 1.0, 20.0)
        self.assertAlmostEqual(born_expansion_parameter(solved), 0.1, places=12)
        self.assertAlmostEqual(threshold_velocity_km_s(solved), 20.0, places=10)
        sigma = sigma_T_born(np.array([100.0]), solved)[0]
        self.assertAlmostEqual(float(sigma), 1.0, places=10)


class B1938MassFitTests(unittest.TestCase):
    def test_shell_projection_recovers_uniform_sphere_mass(self):
        radius_pc = np.array([1.0, 2.0, 3.0])
        density = np.full(3, 2.5)
        measured = projected_enclosed_mass(
            radius_pc, density, 3.0, r_unit="pc", rho_unit="Msun_pc3"
        )
        expected = (4.0 * np.pi / 3.0) * 3.0**3 * 2.5
        self.assertAlmostEqual(measured, expected, places=10)

    def test_shell_projection_matches_uniform_cylinder_formula(self):
        radius_pc = np.array([1.0, 2.0, 3.0])
        density = np.full(3, 2.5)
        projected_radius = 1.2
        measured = projected_enclosed_mass(
            radius_pc, density, projected_radius,
            r_unit="pc", rho_unit="Msun_pc3",
        )
        expected_volume = (4.0 * np.pi / 3.0) * (
            3.0**3 - (3.0**2 - projected_radius**2) ** 1.5
        )
        self.assertAlmostEqual(measured, expected_volume * 2.5, places=10)

    def test_exact_common_mass_rescaling_has_zero_chi2(self):
        fit = absolute_mass_fit(M_20PC / 2.5, M_90PC / 2.5)
        self.assertAlmostEqual(fit["mu"], 2.5, places=12)
        self.assertAlmostEqual(fit["chi2"], 0.0, places=12)
        self.assertAlmostEqual(fit["M_inner_phys"], M_20PC, places=6)
        self.assertAlmostEqual(fit["M_outer_phys"], M_90PC, places=6)

    def test_ratio_mismatch_is_penalized_by_absolute_fit(self):
        fit = absolute_mass_fit(M_20PC, 0.8 * M_90PC)
        self.assertGreater(fit["chi2"], 0.0)
        self.assertNotEqual(fit["residual_inner_sigma"], 0.0)
        self.assertNotEqual(fit["residual_outer_sigma"], 0.0)

    def test_best_snapshot_diagnostic_selects_minimum_chi2(self):
        rows = [
            {"snapshot_idx": 0, "mass_chi2": 4.0},
            {"snapshot_idx": 1, "mass_chi2": 0.25},
        ]
        self.assertEqual(best_mass_fit(rows)["snapshot_idx"], 1)
        self.assertIsNone(best_mass_fit([]))

    def test_log_profile_interpolation_preserves_power_law(self):
        radius = np.array([0.1, 1.0, 10.0])
        profile = radius ** -2
        self.assertAlmostEqual(_log_interp_profile(radius, profile, 2.0), 0.25)

    def test_snapshot_mass_fit_refines_beyond_radius_grid(self):
        radius = np.logspace(-3, 2, 192)
        density = 1.0 / (radius * (1.0 + radius) ** 2)
        data = {"r": radius, "rho": density, "t": 0.0}
        coarse = diagnose_snapshot(
            data, 1.0, 1.0, 1.0, r_s_kpc=1.0, n_scan=16
        )
        refined = min(coarse, key=lambda row: row["mass_chi2"])
        self.assertLess(refined["mass_chi2"], 1e-10)


class ParameterScanTests(unittest.TestCase):
    def test_scan_point_preserves_target_sigma_and_threshold(self):
        point = _model_point("M1_dark_photon_massive", 0.1, 1.0, 20.0, 1.0)
        self.assertTrue(point["born_valid"])
        self.assertAlmostEqual(point["eta_B"], 0.1, places=12)
        self.assertAlmostEqual(point["sigma_m_100_kms"], 1.0, places=10)
        self.assertAlmostEqual(point["threshold_actual_kms"], 20.0, places=10)

    def test_scan_value_parser_accepts_benchmark_threshold(self):
        self.assertEqual(_parse_float_list("0.1, 1, 10"), [0.1, 1.0, 10.0])
        self.assertEqual(_parse_float_list("none,20", allow_none=True), [None, 20.0])

    def test_direct_time_parser_sorts_and_deduplicates(self):
        self.assertEqual(_parse_times("0.5, 0.1, 0.5"), [0.1, 0.5])
        self.assertEqual(_parse_times("0,0.01", allow_initial=True), [0.0, 0.01])

    def test_physical_activity_parser_decodes_final_point_id(self):
        self.assertEqual(
            parse_point_id("M1_eta0p1_sigma0p1_vstar100"),
            (0.1, 0.1, 100.0),
        )

    def test_direct_time_scan_selects_nearest_saved_snapshot(self):
        import pandas as pd

        frame = pd.DataFrame({
            "snapshot_idx": [0, 0, 1, 1, 2, 2],
            "snapshot_time_gyr": [0.0, 0.0, 0.1, 0.1, 0.3, 0.3],
        })
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "diagnostics.csv"
            frame.to_csv(path, index=False)
            snapshot_idx, time_gyr = _nearest_snapshot(path, 0.24)
            evolved_idx, evolved_time = _nearest_snapshot(path, 0.0)
            initial_idx, initial_time = _nearest_snapshot(path, 0.0, include_initial=True)
        self.assertEqual(snapshot_idx, 2)
        self.assertAlmostEqual(time_gyr, 0.3)
        self.assertEqual(evolved_idx, 1)
        self.assertAlmostEqual(evolved_time, 0.1)
        self.assertEqual(initial_idx, 0)
        self.assertAlmostEqual(initial_time, 0.0)

    def test_direct_control_comparison_matches_point_and_time(self):
        import pandas as pd

        common = {
            "scan_point_id": ["M1_eta0p1_sigma1_vstar5"],
            "requested_source_time_gyr": [0.1],
            "status": ["complete"],
            "direct_mass_ratio": [0.36],
            "direct_vmax_kms": [5.3],
            "direct_steps": [100],
        }
        dissipative = pd.DataFrame({
            **common,
            "direct_mass_chi2": [0.02],
            "direct_max_cooling_code": [1e-28],
        })
        dissipative = pd.concat([
            dissipative,
            dissipative.assign(
                requested_source_time_gyr=0.05,
                direct_mass_chi2=0.019,
            ),
        ], ignore_index=True)
        elastic = pd.DataFrame({
            **common,
            "direct_mass_chi2": [0.02],
            "direct_max_cooling_code": [0.0],
        })
        comparison = compare_frames(dissipative, elastic).iloc[0]
        self.assertEqual(comparison["delta_chi2"], 0.0)
        self.assertEqual(comparison["delta_mass_ratio"], 0.0)
        self.assertTrue(comparison["identical_mass_observables"])

    def test_isolated_nfw_extrapolation_closes_overdensity_identity(self):
        result = isolated_nfw_extrapolation(10.0, 50.0, 0.881)
        concentration = result["c_delta_isolated_extrapolation"]
        delta_c = (
            200.0 / 3.0 * concentration**3 / nfw_mass_factor(concentration)
        )
        recovered_rho_s = delta_c * result["rho_critical_msun_pc3"]
        self.assertAlmostEqual(recovered_rho_s, 50.0, places=8)
        self.assertGreater(result["r_delta_over_modeled_rmax"], 1.0)


class EmissionKernelTests(unittest.TestCase):
    def test_massless_kernel_has_expected_velocity_scaling(self):
        model = benchmark_models()["M3_massless_control"]
        q100 = radiated_energy_cross_section(100.0, model)
        q200 = radiated_energy_cross_section(200.0, model)
        self.assertGreater(q100, 0.0)
        self.assertAlmostEqual(q200 / q100, 4.0, places=5)
        cooling = microscopic_cooling_sigma_m(np.array([50.0, 200.0]), model)
        self.assertAlmostEqual(float(cooling[1] / cooling[0]), 1.0, places=5)

    def test_massive_kernel_enforces_emission_threshold(self):
        model = benchmark_models()["M1_dark_photon_massive"]
        self.assertEqual(radiated_energy_cross_section(100.0, model), 0.0)
        self.assertGreater(radiated_energy_cross_section(1000.0, model), 0.0)

    def test_massless_kernel_matches_published_rate(self):
        model = benchmark_models()["M3_massless_control"]
        m_chi = model.m_chi
        temperature = m_chi * (100.0 / 299792.458) ** 2
        relative_variance = 2.0 * temperature / m_chi
        nodes, weights = np.polynomial.legendre.leggauss(64)
        v_nat = 0.5 * (nodes + 1.0) * 0.02
        weights = 0.5 * 0.02 * weights
        relative_pdf = (
            np.sqrt(2.0 / np.pi)
            * v_nat**2
            / relative_variance**1.5
            * np.exp(-v_nat**2 / (2.0 * relative_variance))
        )
        q = np.array([
            radiated_energy_cross_section(v / (1.0 / 299792.458), model)
            for v in v_nat
        ])
        numerical_rate = 0.5 * np.sum(weights * relative_pdf * v_nat * q)
        analytic_rate = (
            (44.0 - 3.0 * np.pi**2)
            / (m_chi**2.5 * np.pi**3.5)
            * temperature**1.5
            * model.g_eff**6
            * (5.0 / 144.0)
        )
        self.assertAlmostEqual(numerical_rate / analytic_rate, 1.0, places=3)


if __name__ == "__main__":
    unittest.main()
