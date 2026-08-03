import sys
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src" / "P7_lens_joint"))

from uvfits_data import decode_aips_baseline, weighted_parallel_hand_average
from visibility_likelihood import (
    VisibilityNoiseModel,
    band_averaged_point_visibility,
    classify_noise_bins,
    complex_gaussian_log_likelihood,
    subchannel_frequencies_hz,
)
from host_tidal_prior import (
    FreePseudoJaffeSummary,
    PublishedTidalReference,
    minimum_current_radius_pc,
    tidal_radius_power_law,
    weighted_quantile,
)
from lens_forward_model import (
    BandIntegratedNufft,
    ImageGrid,
    SourceGrid,
    bilinear_source_operator,
    marginalized_source_log_likelihood_explicit,
    pseudo_jaffe_projected_mass,
)


class UVFitsTests(unittest.TestCase):
    def test_standard_aips_baseline_decoding(self):
        first, second = decode_aips_baseline(np.array([258.0, 517.0, 5651.0]))
        np.testing.assert_array_equal(first, np.array([1, 2, 22]))
        np.testing.assert_array_equal(second, np.array([2, 5, 19]))

    def test_parallel_hand_average_uses_inverse_variance(self):
        visibility, weight = weighted_parallel_hand_average(
            np.array([[2.0, 4.0]]),
            np.array([[1.0, 3.0]]),
            np.array([[1.0, 3.0]]),
        )
        self.assertAlmostEqual(visibility[0].real, 3.5)
        self.assertAlmostEqual(visibility[0].imag, 2.5)
        self.assertAlmostEqual(weight[0], 4.0)

    def test_parallel_hand_average_ignores_flagged_samples(self):
        visibility, weight = weighted_parallel_hand_average(
            np.array([[2.0, 100.0]]),
            np.array([[1.0, 100.0]]),
            np.array([[2.0, 0.0]]),
        )
        self.assertEqual(visibility[0], 2.0 + 1.0j)
        self.assertEqual(weight[0], 2.0)


class VisibilityLikelihoodTests(unittest.TestCase):
    def test_subchannel_centres_cover_averaged_if(self):
        frequencies = subchannel_frequencies_hz(1.65e9)
        self.assertEqual(frequencies.shape, (32,))
        self.assertAlmostEqual(float(np.mean(frequencies)), 1.65e9)
        self.assertAlmostEqual(float(np.diff(frequencies).mean()), 0.25e6)
        self.assertAlmostEqual(float(frequencies[0]), 1.646125e9)
        self.assertAlmostEqual(float(frequencies[-1]), 1.653875e9)

    def test_band_averaged_point_visibility_matches_direct_average(self):
        uu_seconds = 0.031
        vv_seconds = -0.014
        l_radians = 2.3e-6
        m_radians = -1.2e-6
        frequencies = subchannel_frequencies_hz(1.65e9)
        expected = np.mean(np.exp(
            -2j * np.pi * frequencies
            * (uu_seconds * l_radians + vv_seconds * m_radians)
        ))
        actual = band_averaged_point_visibility(
            uu_seconds, vv_seconds, l_radians, m_radians, 1.65e9
        )
        self.assertAlmostEqual(actual.real, expected.real)
        self.assertAlmostEqual(actual.imag, expected.imag)

    def test_complex_gaussian_likelihood_includes_two_components(self):
        observed = np.array([1.0 + 2.0j])
        model = np.array([0.0 + 0.0j])
        sigma = 2.0
        expected = -0.5 * 5.0 / 4.0 - np.log(2.0 * np.pi * 4.0)
        self.assertAlmostEqual(
            complex_gaussian_log_likelihood(observed, model, sigma), expected
        )

    def test_noise_classification_and_lookup(self):
        frame = __import__("pandas").DataFrame({
            "baseline_code": [258, 3350, 3350, 1000],
            "time_bin": [0, 0, 1, 0],
            "complex_difference_count": [128, 128, 128, 128],
            "adjacent_sigma_jy": [0.01, 0.02, 2.0, 0.10],
        })
        classified = classify_noise_bins(frame, robust_z_threshold=0.5)
        self.assertTrue(bool(classified.loc[0, "published_triangle_flag"]))
        self.assertTrue(bool(classified.loc[2, "data_driven_rfi_flag"]))
        model = VisibilityNoiseModel.from_frame(
            classified, start_jd=2450000.0, apply_classification=False
        )
        sigma, retained = model.lookup(
            np.array([258, 3350, 3350, 1000, 9999]),
            2450000.0 + np.array([0.0, 0.0, 30.0, 0.0, 0.0]) / 1440.0,
        )
        np.testing.assert_array_equal(
            retained, np.array([False, True, False, True, True])
        )
        self.assertAlmostEqual(sigma[1], 0.02)
        self.assertAlmostEqual(sigma[3], 0.10)
        self.assertAlmostEqual(sigma[4], model.global_fallback)


class HostTidalPriorTests(unittest.TestCase):
    def test_published_normalization_is_recovered(self):
        reference = PublishedTidalReference()
        tidal = tidal_radius_power_law(
            reference.projected_radius_pc, reference.total_mass_msun,
            reference=reference,
        )
        self.assertAlmostEqual(float(tidal), reference.tidal_radius_pc)

    def test_minimum_radius_inverts_tidal_scaling(self):
        reference = PublishedTidalReference()
        imaging = FreePseudoJaffeSummary()
        radius = minimum_current_radius_pc(
            imaging.tidal_radius_pc, imaging.total_mass_msun,
            reference=reference,
        )
        recovered = tidal_radius_power_law(
            radius, imaging.total_mass_msun, reference=reference
        )
        self.assertAlmostEqual(float(recovered), imaging.tidal_radius_pc)
        self.assertGreater(float(radius), 4.0 * reference.projected_radius_pc)

    def test_weighted_quantile_tracks_dominant_sample(self):
        result = weighted_quantile(
            np.array([1.0, 2.0, 10.0]), (0.5,), np.array([0.01, 0.98, 0.01])
        )
        self.assertAlmostEqual(float(result[0]), 2.0, places=1)


class LensForwardModelTests(unittest.TestCase):
    def test_pseudo_jaffe_mass_is_monotonic_and_finite(self):
        radii = np.array([0.0, 10.0, 100.0, 1.0e8])
        mass = pseudo_jaffe_projected_mass(radii, 2.0e6, 100.0)
        self.assertEqual(mass[0], 0.0)
        self.assertTrue(np.all(np.diff(mass) > 0))
        self.assertAlmostEqual(mass[-1] / 2.0e6, 1.0, places=5)

    def test_bilinear_source_operator_preserves_constant_brightness(self):
        source = SourceGrid(4, 4, 0.1)
        beta_x = np.array([[-0.05, 0.05], [-0.05, 0.05]])
        beta_y = np.array([[-0.05, -0.05], [0.05, 0.05]])
        operator = bilinear_source_operator(
            beta_x, beta_y, source, image_pixel_area_arcsec2=0.01
        )
        image_flux = operator @ np.ones(16)
        np.testing.assert_allclose(image_flux, 0.01)

    def test_nufft_forward_matches_direct_sum_and_adjoint(self):
        grid = ImageGrid(4, 4, 0.01, center_x_arcsec=0.07, center_y_arcsec=-0.03)
        image = np.arange(16, dtype=float).reshape(4, 4) / 100.0
        uu = np.array([0.001, -0.002, 0.003])
        vv = np.array([-0.0015, 0.0025, 0.0005])
        frequency = np.array([1.65e9, 1.66e9, 1.67e9])
        operator = BandIntegratedNufft(
            grid, bandwidth_hz=1.0, channel_count=1, tolerance=1.0e-12
        )
        actual = operator.forward(image, uu, vv, frequency)
        xx, yy = grid.coordinates()
        expected = np.array([
            np.sum(image * np.exp(
                -2j * np.pi * f
                * (u * xx * np.deg2rad(1.0 / 3600.0)
                   + v * yy * np.deg2rad(1.0 / 3600.0))
            ))
            for u, v, f in zip(uu, vv, frequency)
        ])
        np.testing.assert_allclose(actual, expected, rtol=1e-10, atol=1e-11)
        values = np.array([0.5 + 0.2j, -0.1 + 0.4j, 0.3 - 0.7j])
        forward_inner = np.vdot(values, actual)
        adjoint_inner = np.vdot(operator.adjoint(values, uu, vv, frequency), image)
        self.assertAlmostEqual(forward_inner.real, adjoint_inner.real, places=10)
        self.assertAlmostEqual(forward_inner.imag, adjoint_inner.imag, places=10)

    def test_explicit_marginalization_matches_scalar_formula(self):
        data = np.array([1.0 + 0.5j, 0.2 - 0.1j])
        design = np.array([[1.0 + 0.0j], [0.5 + 0.2j]])
        sigma = np.array([0.3, 0.4])
        prior = np.array([[2.0]])
        result = marginalized_source_log_likelihood_explicit(
            data, design, sigma, prior
        )
        inverse_variance = sigma**-2
        hessian = float(np.sum(np.abs(design[:, 0])**2 * inverse_variance) + 2.0)
        rhs = float(np.real(np.sum(
            design[:, 0].conj() * data * inverse_variance
        )))
        self.assertAlmostEqual(result.source_map[0], rhs / hessian)
        expected = -np.sum(np.log(2.0 * np.pi * sigma**2)) - 0.5 * (
            np.sum(np.abs(data) ** 2 * inverse_variance)
            - rhs**2 / hessian + np.log(hessian) - np.log(2.0)
        )
        self.assertAlmostEqual(result.log_likelihood, expected)


if __name__ == "__main__":
    unittest.main()
