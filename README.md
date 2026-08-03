# Gravothermal dSIDM follow-up

This repository studies velocity-dependent dissipative SIDM with a spherical
gravothermal fluid solver, following Schmidt, Fischer, and Garny (2026).

## Current scientific status

The cooling implementation was corrected to use

```text
C_vol = (8/sqrt(pi)) (sigma_T/m_chi) (r_diss - 1) rho^2 nu^3
```

and the thermal average now uses the relative-speed Maxwell distribution for
two identical populations. Production outputs generated before these changes
are legacy artifacts and must not be used as quantitative results. In
particular, the old collapse times, B1938 match counts, resimulation offsets,
and exclusion figures require regeneration.

The current code now also evaluates the long-range identical-fermion emission
kernels from Lankester-Broche & Pradler (2026), Eqs. (3.43) and (3.47), and
uses their energy-weighted integral in the fluid cooling term for the `chi-V`
and `chi-phi` channels. The kernel is validated against Eq. (5.14a) in the
regression suite. This does not remove the Born-validity restriction: the
nominal calibrated M1/M2 points remain outside the controlled Born region.
The explicit Born-valid scan shows that the largest cross sections at
100 km/s and at the lenient `eta_B < 1` boundary are approximately
`3.23e-4 cm^2/g` (M1) and `9.72e-3 cm^2/g` (M2), far below the legacy
`50 cm^2/g` calibration.

At fixed threshold velocity, lowering the DM mass opens an astrophysically
relevant Born-valid window. The `eta_B=0.1` scan gives approximately
`m_chi=0.10--0.32 GeV` for M1 and `0.22--0.89 GeV` for M2 when requiring
`0.1 <= sigma/m <= 10 cm^2/g` at 100 km/s. The corresponding machine-readable
result is `data/born_valid_mass_hierarchy_scan.csv`.

The B1938 rescaling code now evaluates the two absolute projected masses with
a joint chi-square after the ratio preselection. Results under
`data/P2_born_valid_regression` are short interface diagnostics only; they are
not production lensing constraints.

The GADGET/N-body work under `src/P5_gadget` is outside the current analysis
and is not used to validate the fluid solver.

## Environment

Create a clean Python environment and install the core dependencies:

```powershell
python -m venv .venv
.venv\Scripts\python -m pip install --upgrade pip
.venv\Scripts\python -m pip install -r requirements.txt
```

The upstream solver can fall back to pure Python when Numba is unavailable,
which is sufficient for correctness tests. Production halo evolution should
use the Numba-enabled environment from `requirements.txt` because the fallback
is substantially slower.

Run the core regression tests:

```powershell
.venv\Scripts\python -m unittest discover -s tests -v
```

Generate the Born-valid parameter map and the short M1/M2 fluid regression:

```powershell
.venv\Scripts\python src\cross_sections\born_valid_scan.py
.venv\Scripts\python src\fluid_runner\run_born_valid.py --eta 0.1 --steps 25 --n-shells 48
.venv\Scripts\python src\fluid_runner\run_born_valid.py --eta 0.1 --target-sigma-100 1 --steps 25 --n-shells 48 --output-root data\P2_born_controlled_sigma1_regression_new
```

Run the resumable long-time Born-valid experiment and its diagnostics:

```powershell
.venv\Scripts\python src\fluid_runner\run_born_long.py --eta 0.1 --target-sigma-100 1 --t-end-gyr 0.5 --checkpoint-gyr 0.05 --n-shells 48 --output-root data\P2_born_long_eta0p1_sigma1
.venv\Scripts\python src\P3_rescaling\compare_fluid_runs.py --run baseline=data\P2_born_long_eta0p1_sigma1 --output data\P2_fluid_convergence_summary.csv
```

Rerunning the long-run command resumes from the latest checkpoint. The
snapshot diagnostics use the exact projected volume of each piecewise-constant
fluid shell. For a direct B1938 check, select the best evolved snapshot and
resimulate its physical parameters without assuming the rescaling symmetry:

```powershell
.venv\Scripts\python src\P5_resim\resim_born_point.py data\P2_born_long_eta0p1_sigma1\M1_eta0p1_sigma1 data\P5_born_direct_shell_projection\M1_best_n192 --n-shells 192 --t-epsilon 0.005
```

To audit the dependence on the selected evolved time, combine direct
resimulations of several source snapshots at two shell resolutions:

```powershell
.venv\Scripts\python src\P5_resim\summarize_time_scan.py `
  --root n96=data\P5_born_direct_time_scan_n96 `
  --root n192=data\P5_born_direct_time_scan_n192 `
  --output data\P5_born_direct_time_scan_summary.csv
```

The current scan finds its best direct physical fit near the 0.1 Gyr source
snapshot: `chi2_direct ~= 0.024` for both M1 and M2, with the 96-to-192-shell
change below `4.2e-4` in chi-square. This is a candidate-level diagnostic, not
yet a likelihood scan over the full microscopic parameter space.

The first coarse physical scan is driven by:

```powershell
.venv\Scripts\python src\P6_parameter_scan\run_parameter_scan.py `
  --eta-values 0.1 --sigma-values 0.1,1,10 `
  --threshold-values 20,100,500 --t-end-gyr 0.02 `
  --n-shells 48 --direct --resume `
  --output-root data\P6_parameter_scan
```

Its 18 points are screening data only. The direct columns are physical
resimulations; the `best_posthoc_*` columns are retained separately and must
not be interpreted as direct likelihood values.

The first long-candidate refinement also establishes a resolution requirement:
the 48-shell source run can produce an apparently excellent
`chi2_direct ~= 0.0145`, but the same 0.5 Gyr source state gives about
`0.14--0.15` with 96/192-shell source runs. The high-resolution source runs
instead select an earlier `~0.05 Gyr` candidate with
`chi2_direct ~= 0.0248` at 192 shells. Therefore 48 shells are screening-only;
final parameter-space claims require at least 96 shells and a fixed-resolution
time scan.

The corresponding audit tables are `data/P6_parameter_scan_long_summary.csv`
and `data/P6_parameter_scan_convergence.csv`.

The eta extension uses the same scanner at fixed 96-shell resolution:

```powershell
.venv\Scripts\python src\P6_parameter_scan\run_parameter_scan.py `
  --eta-values 0.03,0.1,0.3 --sigma-values 0.1,1,10 `
  --threshold-values 20,100,500 --t-end-gyr 0.5 --n-shells 96 `
  --resume --output-root data\P6_eta_scan_n96
```

The dedicated `v_star=100 km/s` branch is stored in
`data/P6_eta_scan_n96_vstar100/parameter_scan_summary.csv` and contains 18
fluid-complete points (M1/M2, eta = 0.03/0.1/0.3, and sigma/m = 0.1/1/10).
Its direct audit is in `data/P6_eta_direct_best_n96_vstar100`.

For selected source points, direct fits at a common time grid are generated
with:

```powershell
.venv\Scripts\python src\P6_parameter_scan\run_direct_time_scan.py `
  --scan-summary data\P6_eta_scan_n96\direct_candidate_sources.csv `
  --output-root data\P6_eta_direct_time_n96 `
  --source-times-gyr 0.05,0.1,0.2,0.3,0.4,0.5 --n-shells 96 --resume
```

On the earlier common grid (0.05--0.5 Gyr), the eta scan finds its most stable
direct candidate at `sigma/m ~= 0.1 cm^2/g`, `v_star=20 km/s`, with
`chi2_direct ~= 0.0244` at 96 shells and `0.0228` at 192 shells. The same
candidate is obtained for eta values 0.03, 0.1, and 0.3 within the current
precision. These numbers are screening results because the grid does not
include the initial state.

The current final analysis fixes `eta_B = 0.1` and uses 192 shells throughout.
It covers `sigma/m(100 km/s) = 0.1,1,10 cm^2/g` at `v_star=1,5,20 km/s`, plus
`sigma/m=0.1 cm^2/g` controls at 100 and 500 km/s. The source scans are:

- `data/P6_deep_threshold_n192` (`v_star=1 km/s`)
- `data/P6_low_threshold_n192` (`v_star=5 km/s`)
- `data/P6_high_sigma_vstar20_n192` (`v_star=20 km/s`, sigma 1 and 10)
- `data/P6_final_n192_early` (`v_star=20,100,500 km/s`, sigma 0.1)

Before direct resimulation, regenerate continuous-refined mass diagnostics.
The initial logarithmic radius grid locates a branch, then a bounded optimizer
removes discrete `r2D/r_s` switching artifacts:

```powershell
.venv\Scripts\python src\P6_parameter_scan\refine_scan_diagnostics.py `
  --scan-summary data\P6_final_n192_early\parameter_scan_summary.csv `
  --diagnostic-filename b1938_snapshot_diagnostics_refined.csv `
  --output data\P6_final_n192_early\refined_diagnostics_summary.csv

.venv\Scripts\python src\P6_parameter_scan\run_direct_time_scan.py `
  --scan-summary data\P6_final_n192_early\parameter_scan_summary.csv `
  --output-root data\P6_final_direct_time_n192_refined `
  --source-times-gyr 0,0.01,0.02,0.03,0.04,0.05,0.1 `
  --n-shells 192 --t-epsilon 0.005 --include-initial --resume `
  --diagnostic-filename b1938_snapshot_diagnostics_refined.csv
```

The four refined direct summaries contain 154 complete rows and no failures.
At `t_source=0`, two continuous NFW scales interpolate the two measured masses,
giving `chi2_direct ~= 2.1e-16`; this is not an interaction detection. The
time-dependent residuals are smooth after continuous refinement and M1/M2 are
identical at saved precision.

The physical-frame threshold/cooling audit is reproducible with:

```powershell
.venv\Scripts\python src\P6_parameter_scan\audit_physical_activity.py `
  --scan-summary data\P6_low_threshold_direct_time_n192_refined\direct_time_scan_summary.csv `
  --output data\P6_low_threshold_direct_time_n192_refined\physical_activity_audit.csv `
  --n-shells 192
```

The refined direct halo has `v_max ~= 5.34 km/s`. Thus `v_star=1` and 5 km/s
are physically open, while 20, 100, and 500 km/s are below threshold. The
shortest exact local cooling time is nevertheless `8.69e22 Gyr`, at
`v_star=5 km/s`, `sigma/m=10 cm^2/g`, M1. At `v_star=1 km/s`, where the
threshold is exceeded by a factor 5.34, the shortest cooling time is
`5.02e23 Gyr`.

Cooling-disabled controls for all 12 open-threshold candidates at
`t_source=0.1 Gyr` have exactly the same mass ratio, chi-square, velocity, and
step count at saved precision. Their comparison tables are:

- `data/P6_low_threshold_elastic_n192_t01_refined/elastic_comparison.csv`
- `data/P6_deep_threshold_elastic_n192_t01_refined/elastic_comparison.csv`

The fitted initial profile has `r_s=10.00 pc` and
`rho_s=55.88 Msun/pc^3`. Its isolated-NFW extrapolation is
`c200=219`, `r200=2.19 kpc`, and `M200=3.09e6 Msun`; since the modeled domain
ends at `50 r_s`, this is an extrapolation diagnostic rather than a stripped-
subhalo concentration measurement. The machine-readable result is
`data/P6_final_direct_time_n192_refined/nfw_isolated_extrapolation.csv`.

The two extended paper figures are regenerated from all four refined scans:

```powershell
.venv\Scripts\python src\P6_parameter_scan\plot_extended_paper_results.py `
  --time-csv data\P6_final_direct_time_n192_refined\direct_time_scan_summary.csv data\P6_low_threshold_direct_time_n192_refined\direct_time_scan_summary.csv data\P6_deep_threshold_direct_time_n192_refined\direct_time_scan_summary.csv data\P6_high_sigma_vstar20_direct_time_n192_refined\direct_time_scan_summary.csv `
  --audit-csv data\P6_final_direct_time_n192_refined\physical_activity_audit.csv data\P6_low_threshold_direct_time_n192_refined\physical_activity_audit.csv data\P6_deep_threshold_direct_time_n192_refined\physical_activity_audit.csv data\P6_high_sigma_vstar20_direct_time_n192_refined\physical_activity_audit.csv `
  --time-output figures\P6_extended_time_scan.pdf `
  --activity-output figures\P6_extended_activity_audit.pdf
```

The same physical-state audit now evaluates the cooling and local dynamical
times at the fastest-cooling shell, using
`t_dyn=(4*pi*G*rho)^(-1/2)`.  The dense, initial-state failure map is generated
without additional halo evolution:

```powershell
.venv\Scripts\python src\P6_parameter_scan\build_failure_map.py `
  --time-csv data\P6_final_direct_time_n192_refined\direct_time_scan_summary.csv data\P6_low_threshold_direct_time_n192_refined\direct_time_scan_summary.csv data\P6_deep_threshold_direct_time_n192_refined\direct_time_scan_summary.csv data\P6_high_sigma_vstar20_direct_time_n192_refined\direct_time_scan_summary.csv `
  --audit-csv data\P6_final_direct_time_n192_refined\physical_activity_audit.csv data\P6_low_threshold_direct_time_n192_refined\physical_activity_audit.csv data\P6_deep_threshold_direct_time_n192_refined\physical_activity_audit.csv data\P6_high_sigma_vstar20_direct_time_n192_refined\physical_activity_audit.csv `
  --output data\P6_failure_map\cooling_failure_map.csv `
  --figure figures\P6_cooling_failure_map.pdf `
  --eta 0.1 --n-shells 192 --vmax-kms 5.343777930616317
```

This produces 4018 Born-controlled post-processing points over
`0.1 <= v_star/v_max <= 100` and
`0.1 <= sigma/m(100 km/s) <= 10 cm^2/g`.  The minimum local ratio is
`t_cool/t_dyn = 8.32e26` (M1); the M2 minimum is `1.18e27`.

The initial-halo prior audit uses the Dutton--Maccio Planck concentration
relation and a smooth BMO-truncated NFW profile.  The baseline profiles over
`1 <= r_t/r_s <= 100`; a deliberately aggressive sensitivity branch extends
the truncation to `r_t/r_s = 0.5`:

```powershell
.venv\Scripts\python src\P6_parameter_scan\audit_halo_priors.py `
  --isolated-csv data\P6_final_direct_time_n192_refined\nfw_isolated_extrapolation.csv `
  --profile-output data\P6_halo_prior\concentration_profile.csv `
  --summary-output data\P6_halo_prior\halo_prior_summary.csv `
  --figure figures\P6_halo_prior_audit.pdf `
  --tau-min 1 --tau-max 100
```

The isolated extrapolation lies `10.46 sigma_logc` above the median relation.
After profiling over halo mass and the baseline tidal envelope, reaching
`chi2_2M <= 2.30` requires an interpolated concentration offset of
`6.36 sigma_logc`; the aggressive `r_t/r_s >= 0.5` branch lowers this to
`5.69 sigma_logc`.  These are prior-sensitivity diagnostics, not measurements
of an infall concentration or tidal radius.

The old `*_initial` direct directories are retained as provenance for the
pre-refinement 80-point radius-grid selection and must not be used for final
paper numbers.

## P7: public GM068 likelihood and host-tidal prior

The public pipeline-calibrated GM068 visibility file is intentionally kept
under the ignored `data/external/GM068` directory. Its expected SHA-256 is
`4d78a1d81e39d819318714faac7a1901fab907ec5a090a9175c9260967537982`.
The data-quality and likelihood audit is regenerated with:

```powershell
.venv\Scripts\python src\P7_lens_joint\audit_uvfits.py `
  data\external\GM068\gm068_B1938+666.UVDATA.FITS `
  --output-dir data\P7_lens_imaging `
  --figure figures\P7_uvfits_audit.pdf
```

The audit finds 30,173,833 valid RR/LL complex samples before the published
baseline removals. It estimates per-component noise from adjacent
visibilities in 30-minute baseline bins, removes the published EF-JB-WB
triangle, and records data-driven RFI flags separately. Because AIPS `SPLIT`
averaged each 8 MHz IF from 32 channels, the forward operator integrates over
the 32 original channel centres.

A full-data, naturally weighted dirty image can be generated as a phase and
coordinate QA check. It is not deconvolved and must not be used as the lens
likelihood:

```powershell
.venv\Scripts\python src\P7_lens_joint\make_dirty_image.py `
  data\external\GM068\gm068_B1938+666.UVDATA.FITS `
  --noise-csv data\P7_lens_imaging\noise_bins.csv `
  --metadata-json data\P7_lens_imaging\metadata.json `
  --fits-output data\P7_lens_imaging\gm068_dirty_image.fits `
  --summary-output data\P7_lens_imaging\dirty_image_summary.json `
  --figure figures\P7_gm068_dirty_image.pdf
```

The independent lensing core in `src/P7_lens_joint/lens_forward_model.py`
implements the published elliptical power law, `m=3,4` convergence
multipoles, external shear, spherical perturbers, bilinear pixel sources,
band-integrated NUFFTs, and exact finite-dimensional Gaussian source
marginalization. Synthetic forward/adjoint and evidence-normalization tests
are part of `tests/test_uvfits_likelihood.py`.

The host-orbit/tidal sensitivity calculation is regenerated with:

```powershell
.venv\Scripts\python src\P7_lens_joint\run_host_tidal_prior.py `
  --output-dir data\P7_host_tidal_prior `
  --figure figures\P7_host_tidal_prior.pdf `
  --samples 500000 --posterior-draws 20000
```

At the published mean free-PJ mass and truncation radius, the current 3D
radius must be at least 6.28 kpc under the published power-law tidal scaling,
compared with a projected radius of about 1.52 kpc. This is a geometry bound
within the stated tidal model. The orbit-family posteriors are sensitivity
priors, not a cosmological orbit measurement, and currently use the public
one-dimensional PJ summaries because the correlated PRONTO posterior samples
are not available.

The exact published imaging posterior is not reproducible from the archive
alone. The missing private inputs are listed in
`reports/P7_joint_likelihood_readiness.md`; no formal joint Bayes factor or
dSIDM exclusion may be quoted from P7 until those inputs are supplied or an
independent production fit is completed and validated.

The production workflow should not be run until these tests pass. A constant
cross section and constant `r_diss` must preserve the Appendix G rescaling
invariant before velocity-dependent models are interpreted.

## Analysis sequence

1. Validate the elastic upstream solver without dissipation.
2. Validate the corrected cooling closure with constant `sigma/m` and
   constant `r_diss`.
3. Reproduce the exact rescaling invariance for the constant model.
4. Establish the Born-validity mask and select a mass hierarchy with an
   astrophysically relevant cross section.
5. Establish timestep and shell-resolution convergence for long Born-valid runs.
6. Directly resimulate B1938 candidates with physical parameters and both
   absolute projected masses.
7. Extend the differential emission kernel to any additional channels.
8. Rebuild exclusion plots only from regenerated outputs.

The manuscript source is `paper/dsidm_paper.tex`. Disabled legacy sections are
kept in the source for provenance but are not compiled into the current PDF.
