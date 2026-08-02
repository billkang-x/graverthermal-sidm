#!/usr/bin/env bash
#SBATCH --job-name=P1_relax
#SBATCH --partition=amd_512
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --mem=48G
#SBATCH --output=slurm_relax_%j.out
#SBATCH --error=slurm_relax_%j.err

# ============================================================
# Phase A: Elastic relaxation run (IC -> equilibrium snapshot)
# ============================================================
# Input:  ic.dat                  (same IC as original P1, sim-space NFW)
# Output: output_relax/snapshot_000   (new equilibrium IC at t=0.5 code)
#
# This snapshot is then converted to a new IC by convert_restart_to_ic.py
# and used as the t=0 initial condition for Phase B (prod_submit.sh).
# ============================================================

cd "$SLURM_SUBMIT_DIR"
mkdir -p output_relax

source /public1/soft/modules/module.sh
module purge
module load gcc/10.2.0 gsl/2.0 fftw/3.3.8-fjy mpi/openmpi/4.1.1-gcc7.3.0 hdf5/1.8.13-gcc-zyq 2>&1 || true

hdf5_root=/public3/soft/hdf5/1.8.13-gcc-zyq
gsl_root=/public3/soft/gsl/gsl2.0
fftw_root=/public3/soft/fftw/3.3.8-fjy
export LD_LIBRARY_PATH="$hdf5_root/lib:$gsl_root/lib:$fftw_root/lib:${LD_LIBRARY_PATH:-}"
export OMPI_MCA_btl=self,vader,tcp

echo "=== Starting P1 Phase A (relaxation) at $(date) ==="
echo "Job ID: $SLURM_JOB_ID on $(hostname)"
echo "MPI ranks: $SLURM_NTASKS"
echo "Purpose: relax IC to equilibrium, output snapshot at t=0.5"

mpirun -np 16 ~/dsidm_project/source/Gadget4_P1_elastic_control relax_params.txt 2>&1 | tee relax_run.log
rc=$?
echo "=== Finished P1 Phase A at $(date) with exit code $rc ==="
ls -la output_relax/ 2>&1 || true
echo ""
echo "Next step: run convert_restart_to_ic.py to build Phase B IC, then"
echo "           submit prod_submit.sh for the production run."
exit 0
