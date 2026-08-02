# P1 Relaxation Rerun — Execution Workflow

This document describes how to run the two-phase P1 relaxation rerun on the
HPC cluster. The goal is to test whether the -18.8% drop in the original P1
projected-mass-ratio run is caused by initial-condition disequilibrium.

## Overview

```
   Phase A (relax)                  Phase B (production)
   ┌──────────────────┐             ┌──────────────────┐
   │  ic.dat           │             │  ic_equilibrium  │
   │  (raw NFW IC)     │             │  .dat            │
   │       v           │             │  (from Phase A)  │
   │  elastic SIDM     │  snapshot   │       v          │
   │  t = 0 -> 0.5     │ ─────────> │  convert_restart │
   │       v           │   to_ic.py  │       v          │
   │  snapshot_000     │             │  elastic SIDM    │
   │  (at t=0.5)       │             │  t = 0 -> 1.0    │
   └──────────────────┘             │       v          │
                                    │  snapshot_000    │
                                    │  snapshot_001    │
                                    │  snapshot_002    │
                                    └──────────────────┘
                                            v
                                    analyze_relaxation.py
```

## Step 0: Prerequisites

1. **HPC account** on the cluster with the same modules as the original run.
2. **Compiled Gadget4 binaries** in `~/dsidm_project/source/`. You need two
   builds:
   - `Gadget4_P1_relax` — compiled with `relax_Config.sh`
   - `Gadget4_P1_relax` (same binary works for Phase B — the Config differs
     only in the absence of SIDM_DISSIPATIVE, but Phase B uses the same
     elastic-only Config). If you want to reuse the original P1 binary
     (`Gadget4_P1_elastic_control`), that also works because the elastic
     Config is identical — just point the submit scripts at it.

   **Simplest option**: use the original `Gadget4_P1_elastic_control` binary
   for both phases. Both `relax_Config.sh` and `prod_Config.sh` are elastic-
   only and compatible with that binary.

3. **The IC file** `ic.dat` must be present in the run directory. Copy from
   `D:/graverthermal-sidm/data/P5_nbody_verify/ics_sim_space/ic.dat`.

## Step 1: Upload files to HPC

```bash
# From local machine
scp -r D:/graverthermal-sidm/data/P5_nbody_verify/P1_relaxation_rerun \
       <hpc-user>@<hpc-host>:~/dsidm_project/P1_relaxation_rerun/

# Also upload the IC if not already there
scp D:/graverthermal-sidm/data/P5_nbody_verify/ics_sim_space/ic.dat \
    <hpc-user>@<hpc-host>:~/dsidm_project/P1_relaxation_rerun/ic.dat
```

On the HPC, the directory should contain:

```
P1_relaxation_rerun/
  ic.dat                       <- raw NFW IC (sim-space)
  relax_params.txt
  relax_Config.sh
  relax_submit.sh
  prod_params.txt
  prod_Config.sh
  prod_submit.sh
  output_list_relax.txt
  output_list_prod.txt
  convert_restart_to_ic.py
  analyze_relaxation.py
```

## Step 2: Compile the binary (only if reusing is not possible)

If you need a fresh build:

```bash
cd ~/dsidm_project/source
cp ~/dsidm_project/P1_relaxation_rerun/relax_Config.sh Config.sh
make clean && make -j
cp Gadget4 ~/dsidm_project/source/Gadget4_P1_relax
```

For Phase B, the Config is identical (elastic only), so the same binary works.

## Step 3: Run Phase A (relaxation)

```bash
cd ~/dsidm_project/P1_relaxation_rerun
sbatch relax_submit.sh
```

Monitor:

```bash
squeue -u $USER
tail -f slurm_relax_*.out
```

Phase A should take ~30-60 minutes (it runs to t=0.5, half the original time).
When done, `output_relax/snapshot_000` should exist (at t=0.5).

## Step 4: Convert Phase-A snapshot to Phase-B IC

```bash
# On HPC (or locally after downloading the snapshot)
module load python/3.x  # or your preferred Python
python convert_restart_to_ic.py output_relax/snapshot_000 ic_equilibrium.dat
```

This produces `ic_equilibrium.dat` — a Gadget2-binary IC with `time=0` in the
header, ready to be used as `InitCondFile` for Phase B.

## Step 5: Run Phase B (production)

```bash
cd ~/dsidm_project/P1_relaxation_rerun
sbatch prod_submit.sh
```

Monitor:

```bash
squeue -u $USER
tail -f slurm_prod_*.out
```

Phase B should take ~1-2 hours (same as original P1). When done,
`output_prod/snapshot_000`, `snapshot_001`, `snapshot_002` should exist at
t=0, 0.5, 1.0.

## Step 6: Download results

```bash
# From local machine
mkdir -p D:/graverthermal-sidm/data/P5_nbody_verify/P1_relaxation_rerun/output_prod
scp <hpc-user>@<hpc-host>:~/dsidm_project/P1_relaxation_rerun/output_prod/snapshot_* \
    D:/graverthermal-sidm/data/P5_nbody_verify/P1_relaxation_rerun/output_prod/
```

## Step 7: Analyze

```bash
cd D:/graverthermal-sidm/data/P5_nbody_verify/P1_relaxation_rerun
python analyze_relaxation.py
```

This prints a comparison table and a conclusion:

- **Drop < 5%** → IC disequilibrium was the main cause; P3 verification stands.
- **Drop 5-10%** → Partial IC effect; some residual issue remains.
- **Drop > 10%** → IC is not the cause; revisit P3 "pass" verdict.

A CSV `relaxation_comparison.csv` is also written for record-keeping.

## Step 8: Update the final report

Based on the outcome, update
`D:/graverthermal-sidm/data/P5_nbody_verify/final_analysis_report.md`:

- If IC was the cause: note the relaxation test confirms the IC transient
  hypothesis; P3 verdict unchanged.
- If not: downgrade P3 to "inconclusive" and recommend a deeper code/physics
  audit.

## Timing summary

| Phase | Wall time (est.) | Output |
|-------|-----------------|--------|
| A (relax)   | 30-60 min  | 1 snapshot |
| Conversion  | <1 min     | 1 IC file |
| B (prod)    | 1-2 hr     | 3 snapshots |
| Analysis    | <1 min     | 1 CSV + console |
| **Total**   | **~2-3 hr** | |

## Key invariants

- The IC, softening, code units, and SIDM parameters are identical to the
  original P1 run. The ONLY change is that Phase B starts from a relaxed
  particle distribution instead of the raw NFW sample.
- The projected-mass ratio is invariant under the Schmidt 2026 Appendix-G
  rescaling, so simulation-space radii (R_inner=0.847 kpc, R_outer=3.812 kpc)
  are used directly in the analysis.
