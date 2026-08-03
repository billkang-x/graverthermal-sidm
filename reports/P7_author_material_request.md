# Request for B1938+666 modelling products

Subject: Reproducibility request for the GM068 B1938+666 visibility modelling

Dear Dr Powell and Dr Vegetti,

We are carrying out an independent follow-up analysis of the B1938+666
perturber using the public EVN GM068 data. We have reproduced the UVFITS data
audit, implemented the visibility-plane complex Gaussian likelihood, and are
adding a host-orbit/tidal prior to the perturber profile inference. To validate
the baseline analysis before interpreting the new prior, could you please
share any of the following research products that are available?

1. The calibrated visibility dataset used by PRONTO, including whether the
   32 channels in each 8 MHz IF were retained or averaged.
2. The final flag table, especially the exact RFI intervals on GB-HN, GB-OV,
   GB-PT, and LA-PT.
3. The coordinate transform between PRONTO image-plane coordinates and the
   GM068 UVFITS phase centre.
4. The fiducial image mask, adaptive source grid, source regularization
   operator, and regularization convention.
5. Numerical priors and either posterior samples or a best-fit configuration
   for the EPL plus multipoles, shear, object A, and object V.
6. Posterior samples for the free pseudo-Jaffe and NFW/profile models reported
   in Powell et al. (2025) and Vegetti et al. (2026), including evidence logs.
7. The PRONTO likelihood settings needed to reproduce the baseline evidence:
   NUFFT conventions, fast-chi-square configuration, log-determinant
   approximation, MultiNest live points, tolerance, and seed.

We would use these files only for scientific reproducibility and would cite
the original papers and any software/data release you specify. A minimal
package containing the final flags, coordinate registration, macro-model
chain, source grid, and perturber posterior samples would already allow us to
perform the host-prior reweighting without rerunning the full discovery
analysis.

Best regards,

[Name and affiliation]
