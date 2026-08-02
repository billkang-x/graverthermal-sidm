#!/bin/bash
#SBATCH -J dsidm_test
#SBATCH -p all
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -o dsidm_test_%j.out
#SBATCH -e dsidm_test_%j.err
#SBATCH -t 02:00:00

# Load environment
source ~/env.sh

# Build the dissipative Gadget4
cd ~/dsidm_project/source
make -j 16 EXEC=Gadget4_dsidm 2>&1 | tee build_dsidm.log

# Run an isolated halo test
cd ~/dsidm_project/runs/test_elastic
mpirun -np 32 ../../../source/Gadget4_dsidm params.txt 2>&1 | tee run_elastic.log
