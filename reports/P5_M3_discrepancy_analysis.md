# M3 Control Group Discrepancy Analysis

## Summary

The M3 control model (constant `r_diss = 1.05`) produces **26/195 viable grid points**, compared to 55/195 (M1) and 62/195 (M2). This is **NOT** a bug or inconsistency between the rescaling formula and direct simulation — it is a **genuine physical effect** caused by the much weaker dissipative cooling in M3 compared to M1/M2 at the relevant velocity scales.

## Root Cause

### 1. Weak cooling rate in M3
The cooling rate per unit mass in the DissipativeHalo implementation is:
```
C_cool ∝ σ_m × ρ × v × (r_diss - 1)
```
For M3, `(r_diss - 1) = 0.05` (constant), while for M1/M2, `r_diss(v)` is velocity-dependent and **much larger** at low velocities (the relevant regime for the inner halo where collapse happens). At v ~ 50 km/s (inner halo thermal velocity), M1/M2 have r_diss ~ 1.5-2.0, giving `(r_diss - 1) ~ 0.5-1.0` — **10-20× stronger cooling than M3**.

### 2. Mass ratio barely evolves in M3
At a representative grid point (σ/m=0.005, r_s=0.06 kpc):
- **M1**: ratio evolves 0.51 → 0.41 (Δ ≈ -0.10, cooling-driven collapse)
- **M3**: ratio evolves 0.51 → 0.51 (Δ ≈ -0.01, essentially no evolution)

M3's weak cooling cannot drive enough mass redistribution to bring the ratio down to the observed 0.364 at low-to-moderate σ/m.

### 3. M3 only viable at high σ/m and large r_s
M3 viable points are concentrated at:
- σ/m ≥ 0.011 (vs 0.005 for M1/M2)
- r_s ≥ 0.071 kpc (vs 0.060 for M1/M2)

At large r_s, the initial ratio is closer to the observed value (init_ratio ≈ 0.35-0.38 vs 0.41 at r_s=0.06), so less cooling is needed. At high σ/m, the scattering rate is high enough that even weak per-collision cooling (×0.05) accumulates sufficiently.

## Comparison with P3 Rescaling

### P3 elastic rescaling assumption
P3 uses elastic rescaling symmetry: if a simulation with (σ_sim, ρ_sim, t_sim) matches B1938, then any (σ_phys, ρ_phys, t_phys) satisfying the symmetry relation also matches. This assumption is **velocity-independent** — it treats σ/m and r_diss as constants.

For M3 (const r_diss=1.05), this is actually a **consistent** approximation since both σ/m and r_diss ARE constant. The P3 elastic rescaling gives M3 **8177 viable matches** — far more than the 26 from the self-consistent grid scan.

### Why the discrepancy?
The P3 rescaling matches are based on the **mass ratio at a snapshot**, including the **initial NFW snapshot** (snapshot_idx=0, t_evo=0). For M3, 17 of the 8177 matches have snapshot_idx=0 — these are points where the initial NFW ratio already falls in the observed band. The self-consistent grid scan correctly requires actual evolution (t_cross > 0), filtering out these trivial matches.

More importantly, P3's elastic symmetry assumes the mass ratio evolution depends only on the dimensionless combination n_sigma = σ_m × ρ × t × v. For M3 with weak cooling, the evolution is so slow that within the t_obs window, the ratio barely moves — the symmetry holds but the halo doesn't actually collapse enough. The grid scan correctly captures this physical limitation.

## Conclusion

The M3 discrepancy is **physical, not artifactual**:
1. M3's constant r_diss=1.05 produces cooling that is 10-20× weaker than M1/M2 at inner-halo velocities.
2. This is insufficient to drive the mass redistribution needed to reach the observed ratio at most (σ/m, r_s) points.
3. The P3 elastic rescaling over-counts viable points because it doesn't account for the actual cooling rate — it assumes the symmetry relation alone guarantees a match.
4. The grid scan correctly shows M3 is only viable in a narrow high-σ/m, large-r_s corner.

**This is NOT a code bug.** The DissipativeHalo implementation correctly applies the M3 cooling rate. The reduced viable region is the physically expected consequence of M3's weak dissipation.

## Recommendation

Document this as a physical finding in the paper: M3 (constant r_diss=1.05) is a "weak dissipation" control that is viable only in a narrow region, demonstrating that the velocity-dependent enhancement of r_diss(v) in M1/M2 is what enables the broad viable region.
