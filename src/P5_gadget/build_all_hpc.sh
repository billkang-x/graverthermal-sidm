#!/bin/bash
# Build all 3 Gadget4 executables without SPH/COOLING and submit jobs.
# This script runs on the HPC login node directly.

set -e

SOURCE=~/dsidm_project/source
BASE=~/dsidm_project/nbody_verify

# Cancel any running jobs
scancel 40889954 40889955 40889958 2>&1 || true

# Load build environment
source /public1/soft/modules/module.sh
module purge
module load gcc/10.2.0 gsl/2.0 fftw/3.3.8-fjy mpi/openmpi/4.1.1-gcc7.3.0 hdf5/1.8.13-gcc-zyq 2>&1 || true
echo "Modules loaded"

# Backup Config.sh
cp $SOURCE/Config.sh $SOURCE/Config.sh.bak_final 2>/dev/null || true

# Build each executable
build_point() {
    local NAME=$1
    local SIGMA=$2
    local RDISS=$3
    local TEND=$4

    echo ""
    echo "=========================================="
    echo "Building $NAME (sigma/m=$SIGMA, r_diss=$RDISS, t_end=$TEND)"
    echo "=========================================="

    # Write Config.sh
    cat > $SOURCE/Config.sh <<EOF
SELFGRAVITY
NTYPES=3
EVALPOTENTIAL
SIDM
SIDM_SIGMA_OVER_MASS=$SIGMA
SIDM_K_NEIGHBORS=32
SIDM_PMAX=0.1
DOUBLEPRECISION=1
GADGET2_HEADER

# Dissipative fSIDM extension
SIDM_DISSIPATIVE
SIDM_R_DISS=$RDISS
EOF

    echo "--- Config.sh ---"
    cat $SOURCE/Config.sh

    # Clean previous build
    cd $SOURCE
    rm -f Gadget4_$NAME build/*.o build/*.mod 2>/dev/null
    mkdir -p build

    # Build
    echo "--- Building (timeout 600s) ---"
    if timeout 600 make build -j 4 EXEC=Gadget4_$NAME > /tmp/build_$NAME.log 2>&1; then
        echo "Build OK"
        ls -la $SOURCE/Gadget4_$NAME
    else
        echo "Build FAILED (exit $?)"
        echo "--- last 30 lines of build log ---"
        tail -30 /tmp/build_$NAME.log
        return 1
    fi

    # Write params.txt
    cat > $BASE/$NAME/params.txt <<EOF
%paramfile for N-body verification
% Point: $NAME
% SIDM sigma/m = $SIGMA cm^2/g
% SIDM r_diss = $RDISS
% Target evolution time: $TEND Gyr

InitCondFile        ic.dat
OutputDir           output
SnapshotFileBase    snapshot
OutputListFilename  output_list.txt

ICFormat            1
SnapFormat          1

TimeLimitCPU              7200
CpuTimeBetRestartFile     3600
MaxMemSize                2000

TimeBegin           0.0
TimeMax             $TEND

ComovingIntegrationOn    0
Omega0              0
OmegaLambda         0
OmegaBaryon         0
HubbleParam         1.0
Hubble              0
BoxSize             0

OutputListOn              1
TimeBetSnapshot           $TEND
TimeOfFirstSnapshot       $TEND
TimeBetStatistics         0.01
NumFilesPerSnapshot       1
MaxFilesWithConcurrentIO  1

ErrTolIntAccuracy        0.025
CourantFac               0.15
MaxSizeTimestep          0.005
MinSizeTimestep          0.0

TypeOfOpeningCriterion                1
ErrTolTheta                           0.5
ErrTolThetaMax                        0.5
ErrTolForceAcc                        0.005
TopNodeFactor                         3.0
ActivePartFracForNewDomainDecomp      0.01

DesNumNgb                             32
MaxNumNgbDeviation                    2

UnitLength_in_cm         3.0856775814913673e21
UnitMass_in_g            1.989e43
UnitVelocity_in_cm_per_s 1.0e5
GravityConstantInternal  0

SofteningComovingClass0       0.0005
SofteningMaxPhysClass0        0.0005
SofteningComovingClass1       0.0005
SofteningMaxPhysClass1         0.0005
SofteningComovingClass2        0.0005
SofteningMaxPhysClass2         0.0005
SofteningClassOfPartType0     0
SofteningClassOfPartType1     1
SofteningClassOfPartType2     2

ArtBulkViscConst       1.0
MinEgySpec             0
InitGasTemp            0
EOF

    # Write output_list.txt
    echo "$TEND" > $BASE/$NAME/output_list.txt

    # Write submit.sh
    cat > $BASE/$NAME/submit.sh <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=$NAME
#SBATCH --partition=amd_512
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=24G
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err

cd "\$SLURM_SUBMIT_DIR"
mkdir -p output

source /public1/soft/modules/module.sh
module purge
module load gcc/10.2.0 gsl/2.0 fftw/3.3.8-fjy mpi/openmpi/4.1.1-gcc7.3.0 hdf5/1.8.13-gcc-zyq 2>&1 || true

hdf5_root=/public3/soft/hdf5/1.8.13-gcc-zyq
gsl_root=/public3/soft/gsl/gsl2.0
fftw_root=/public3/soft/fftw/3.3.8-fjy
export LD_LIBRARY_PATH="\$hdf5_root/lib:\$gsl_root/lib:\$fftw_root/lib:\${LD_LIBRARY_PATH:-}"
export OMPI_MCA_btl=self,vader,tcp

echo "=== Starting $NAME at \$(date) ==="
echo "Job ID: \$SLURM_JOB_ID on \$(hostname)"
echo "MPI ranks: \$SLURM_NTASKS"

mpirun -np 2 ~/dsidm_project/source/Gadget4_$NAME params.txt 2>&1 | tee run.log
rc=\$?
echo "=== Finished $NAME at \$(date) with exit code \$rc ==="
ls -la output/ 2>&1 || true
exit 0
EOF
    chmod +x $BASE/$NAME/submit.sh

    # Clean previous output
    rm -rf $BASE/$NAME/output $BASE/$NAME/restart $BASE/$NAME/slurm_*.out $BASE/$NAME/slurm_*.err $BASE/$NAME/run.log 2>/dev/null
    mkdir -p $BASE/$NAME/output

    echo "  $NAME setup complete"
}

# Build all 3 points
build_point P1_elastic_control 0.1   1.0  0.68
build_point P2_m3_low_sigma    0.005 1.05 0.07
build_point P3_m3_high_sigma   0.220 1.05 0.10

# Restore Config.sh
cp $SOURCE/Config.sh.bak_sph $SOURCE/Config.sh 2>/dev/null || true
echo "Config.sh restored"

# Submit all jobs
echo ""
echo "=== Submitting jobs ==="
for name in P1_elastic_control P2_m3_low_sigma P3_m3_high_sigma; do
    cd $BASE/$name
    out=$(sbatch submit.sh)
    echo "  $name: $out"
done

echo ""
echo "All done. Poll with: squeue -u scg7816"
