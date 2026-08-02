#!/bin/bash
#SBATCH --job-name=P4_m1_low_sigma
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err

source ~/env.sh

cd $SLURM_SUBMIT_DIR
echo "=== Starting P4_m1_low_sigma at $(date) ==="
echo "Job ID: $SLURM_JOB_ID on $(hostname)"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"

# Run Gadget4 with the params file
mpirun -np 8 ~/dsidm_project/source/Gadget4_dsidm params.txt

echo "=== Finished P4_m1_low_sigma at $(date) ==="
ls -la output/
