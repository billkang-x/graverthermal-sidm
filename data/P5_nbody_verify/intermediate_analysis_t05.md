# N-body Verification: Intermediate Analysis at t_code=0.5

## Status
- Date: 2026-07-27 09:00 CST
- Jobs: P1, P2, P3 all running, currently at t_code≈0.61
- snapshot_000 (t=0) and snapshot_001 (t=0.5) downloaded and analyzed
- snapshot_002 (t=1.0) pending — ETA ~10:00 CST

## Results at t_code=0.5 (t_phys=0.013 Gyr)

| Point | sigma/m | r_diss | N-body ratio | N-body rel | Fluid rel |
|-------|---------|--------|--------------|------------|-----------|
| P1 (elastic) | 3.33 | 1.0 | 0.068563 | 0.8155 (-18.5%) | 0.9599 (-4.0%) |
| P2 (low sigma) | 0.167 | 1.05 | 0.081230 | 0.9662 (-3.4%) | 1.0097 (+1.0%) |
| P3 (high sigma) | 7.33 | 1.05 | 0.069582 | 0.8276 (-17.2%) | 0.9905 (-0.9%) |

Reference (no SIDM) at t=0.5 (est): rel = 1.0118 (+1.18%)

## Key observations

1. **N-body shows much larger changes than fluid model** at the same physical time
2. **Correct ordering**: P2 (low sigma) shows smallest change, as expected
3. **P1 and P3 show similar large decreases** (~17-19%), despite different sigma/m values
4. **The discrepancy is ~4-20x** between N-body and fluid model

## Possible explanations

1. **IC not in hydrostatic equilibrium**: The cored NFW IC may not be in perfect
   equilibrium, causing an initial transient mass redistribution
2. **Rescaled sigma/m too large**: The rescaled values (3.33, 0.167, 7.33) may be
   in a regime where the fluid model breaks down
3. **Cored vs cusped response**: SIDM affects cored profiles differently from cusps
4. **2-body relaxation**: N-body simulation has numerical relaxation effects
5. **Fluid model underestimates early-time SIDM**: The fluid model may not capture
   the rapid initial response of the inner regions

## Next steps

1. Wait for snapshot_002 (t=1.0) to complete (~1 hour)
2. Analyze final state
3. Check if the trend continues or stabilizes
4. Investigate the discrepancy sources
