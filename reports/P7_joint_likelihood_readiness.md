# GM068 imaging likelihood and host-orbit prior: readiness audit

Date: 2026-08-03

## Status

The public GM068 data layer, complex Gaussian noise model, published baseline
mask, bandwidth-integrated visibility operator, pixel-source lens operator,
and a host-orbit/tidal sensitivity calculation are implemented and tested.
An exact reproduction of the Powell et al. / Vegetti et al. PRONTO posterior
is not yet possible from public files. No formal joint evidence or exclusion
has been generated.

## Public data verified

- EVN experiment: GM068, observation DOI
  [10.48717/wch4-m437](https://doi.org/10.48717/wch4-m437).
- Pipeline UVFITS: 608,572,800 bytes; SHA-256
  `4d78a1d81e39d819318714faac7a1901fab907ec5a090a9175c9260967537982`.
- Structure: 2,716,193 groups, eight IFs, RR/LL, 43,459,088 stored complex
  correlation samples.
- Positive finite pipeline weights: 30,173,833 samples (69.4304%).
- Frequency centres: 1.630615--1.686365 GHz. Each public IF is an AIPS
  `SPLIT` average of the original 32 channels over 8 MHz.
- Adjacent-difference noise: 2,255 baseline/time bins at 30-minute cadence.
- Likelihood mask: 76 bins in the published EF-JB-WB triangle, 16 bins with
  insufficient adjacent pairs, and two LA-PT bins at 9.0--10.0 h identified
  by the declared robust RFI rule.

## Implemented likelihood components

`src/P7_lens_joint/visibility_likelihood.py` implements the normalized
complex Gaussian likelihood with one noise standard deviation per real or
imaginary component. It also integrates a sky model over the 32 original
subchannel centres rather than evaluating the averaged UVFITS IF only at its
central frequency.

`src/P7_lens_joint/lens_forward_model.py` implements:

- the published elliptical power-law macro lens;
- the published radial-slope-dependent `m=3` and `m=4` convergence terms;
- external shear;
- singular pseudo-Jaffe or arbitrary tabulated spherical perturbers;
- a bilinear pixelated source operator;
- forward and adjoint band-integrated two-dimensional NUFFTs;
- a positive-definite gradient prior; and
- exact Gaussian source marginalization for a finite design matrix,
  including determinant normalization.

The full-data dirty-image QA uses all 15,923,560 retained Stokes-I samples.
It verifies visibility phases and data integrity but is not deconvolved and is
not a substitute for the lens likelihood.

## Host-orbit/tidal result available now

Powell et al. report `r_t=53+/-1 pc` when the 3D host distance is set equal to
the projected distance. Their free pseudo-Jaffe model gives
`r_t=149+/-18 pc` and `m_tot=(2.82+/-0.26)e6 Msun`. Scaling their own tidal
formula through the power-law host yields, at the mean parameters:

- projected radius: approximately 1.52 kpc (digitized macro geometry);
- minimum current 3D radius: 6.277 kpc;
- minimum absolute line-of-sight offset: 6.090 kpc.

The current-radius posterior medians are 5.97 kpc for the circular upper
envelope, 8.57 kpc for the phase-mixed sensitivity prior, and 14.97 kpc for
the radial-orbit sensitivity prior. Measurement uncertainty allows posterior
draws below the mean-parameter bound. These distributions are prior
sensitivity diagnostics, not formal orbit constraints.

For the phase-mixed family, varying the NFW-like host tracer scale radius from
10 to 30 to 100 kpc changes the current-radius posterior median from 8.15 to
8.57 to 8.92 kpc. The corresponding pericentre medians are 5.86, 5.96, and
6.06 kpc. This particular host-scale choice is therefore subdominant to the
orbital-family uncertainty.

## Inputs required for a formal joint posterior

The following numerical inputs are not supplied by the paper or archive:

- the exact 1--2 h RFI intervals or final flag table used by PRONTO;
- registration of PRONTO model coordinates to the UVFITS phase centre;
- the adaptive source-plane grid, image mask, and source regularization
  matrix for each lens model;
- the numerical macro-model priors, initialization, and posterior chain;
- the correlated perturber profile posterior samples;
- the log-determinant preconditioner and fast-chi-square configuration;
- MultiNest live-point, stopping, seed, and evidence settings; and
- the calibrated unaveraged 32-channel modelling dataset, if it differs from
  the public channel-averaged pipeline UVFITS.

The papers state that PRONTO is private and direct readers to the
corresponding author. The numerical posterior samples are not provided in the
supplementary files; only corner plots and one-dimensional summaries are
public.

## Completion criterion

A result can be called a complete joint imaging/orbit analysis only after one
of these paths is completed:

1. Obtain the PRONTO configuration, flag table, source grid, and posterior
   products from the authors, then reproduce their baseline evidence before
   adding the orbit prior.
2. Complete an independent production inference with the public visibilities,
   sample the global coordinate shift as a nuisance parameter, validate
   synthetic injection recovery and profile/evidence convergence, and release
   the new posterior chain and all configuration files.

Until then, the existing two-radius likelihood and the P7 host result must be
described as compressed and sensitivity analyses, respectively.
