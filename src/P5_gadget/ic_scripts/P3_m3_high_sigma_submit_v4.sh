#!/usr/bin/env bash
#SBATCH --job-name=P3_m3_high_sigma
#SBATCH --partition=amd_512
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=24G
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err

cd "$SLURM_SUBMIT_DIR"
mkdir -p output

source /public1/soft/modules/module.sh
module purge
module load gcc/10.2.0 gsl/2.0 fftw/3.3.8-fjy mpi/openmpi/4.1.1-gcc7.3.0 hdf5/1.8.13-gcc-zyq 2>&1 || true

hdf5_root=/public3/soft/hdf5/1.8.13-gcc-zyq
gsl_root=/public3/soft/gsl/gsl2.0
fftw_root=/public3/soft/fftw/3.3.8-fjy
export LD_LIBRARY_PATH="$hdf5_root/lib:$gsl_root/lib:$fftw_root/lib:${LD_LIBRARY_PATH:-}"
export OMPI_MCA_btl=self,vader,tcp

echo "=== Starting P3_m3_high_sigma at $(date) ==="
echo "Job ID: $SLURM_JOB_ID on $(hostname)"
echo "MPI ranks: $SLURM_NTASKS"
echo "Mem avail: $(free -h | head -2)"
echo "Working dir: $(pwd)"

mpirun -np 2 ~/dsidm_project/source/Gadget4_P3_m3_high_sigma params.txt 2>&1 | tee run.log
rc=$?
echo "=== Finished P3_m3_high_sigma at $(date) with exit code $rc ==="
ls -la output/ 2>&1 || true
exit 0
