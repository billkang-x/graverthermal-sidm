#!/usr/bin/env bash
#SBATCH --job-name=P1_prod
#SBATCH --partition=amd_512
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --mem=48G
#SBATCH --output=slurm_prod_%j.out
#SBATCH --error=slurm_prod_%j.err

# ============================================================
# Phase B: Production run from equilibrium IC
# ============================================================
# Input:  ic_equilibrium.dat      (converted from Phase-A snapshot by
#                                  convert_restart_to_ic.py)
# Output: output_prod/snapshot_000 (t=0, equilibrium)
#         output_prod/snapshot_001 (t=0.5)
#         output_prod/snapshot_002 (t=1.0)
#
# These snapshots are then analyzed by analyze_relaxation.py.
# ============================================================

cd "$SLURM_SUBMIT_DIR"
mkdir -p output_prod

source /public1/soft/modules/module.sh
module purge
module load gcc/10.2.0 gsl/2.0 fftw/3.3.8-fjy mpi/openmpi/4.1.1-gcc7.3.0 hdf5/1.8.13-gcc-zyq 2>&1 || true

hdf5_root=/public3/soft/hdf5/1.8.13-gcc-zyq
gsl_root=/public3/soft/gsl/gsl2.0
fftw_root=/public3/soft/fftw/3.3.8-fjy
export LD_LIBRARY_PATH="$hdf5_root/lib:$gsl_root/lib:$fftw_root/lib:${LD_LIBRARY_PATH:-}"
export OMPI_MCA_btl=self,vader,tcp

echo "=== Starting P1 Phase B (production) at $(date) ==="
echo "Job ID: $SLURM_JOB_ID on $(hostname)"
echo "MPI ranks: $SLURM_NTASKS"
echo "Purpose: production run from equilibrium IC, snapshots at t=0,0.5,1.0"

mpirun -np 16 ~/dsidm_project/source/Gadget4_P1_elastic_control prod_params.txt 2>&1 | tee prod_run.log
rc=$?
echo "=== Finished P1 Phase B at $(date) with exit code $rc ==="
ls -la output_prod/ 2>&1 || true
echo ""
echo "Next step: download output_prod/snapshot_* and run"
echo "           analyze_relaxation.py to compare with original P1 and fluid model."
exit 0
