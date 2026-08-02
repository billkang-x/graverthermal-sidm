#!/bin/bash
# =====================================================================
# Resubmit with working reference IC (cored NFW profile)
#
# The working reference IC has a cored profile (r_c ~ 0.3 kpc) that allows
# reasonable timesteps. We use this IC with our 3 rescaled sigma values.
#
# TimeMax = 1.0 code units (matching working reference, ~4 hours)
# Snapshots at t=0, 0.5, 1.0
#
# The working reference used sigma/m=3.0 and completed to t=1.0.
# Our 3 points use sigma/m = 3.33, 0.167, 7.33.
# =====================================================================

set -e

BASE=~/dsidm_project/nbody_verify_sim
SOURCE=~/dsidm_project/source

# --- Cancel any running jobs ---
scancel --user=scg7816 --state=RUNNING,PENDING 2>/dev/null || true
echo "Cancelled previous jobs"

# --- Load build environment ---
source /public1/soft/modules/module.sh
module purge
module load gcc/10.2.0 gsl/2.0 fftw/3.3.8-fjy mpi/openmpi/4.1.1-gcc7.3.0 hdf5/1.8.13-gcc-zyq 2>&1 || true
echo "Modules loaded"

# --- Backup Config.sh ---
cp $SOURCE/Config.sh $SOURCE/Config.sh.bak_ref 2>/dev/null || true

# =====================================================================
# Build function (using working reference's Config.sh)
# =====================================================================
build_point() {
    local NAME=$1
    local SIGMA=$2
    local RDISS=$3

    echo ""
    echo "=========================================="
    echo "Building $NAME"
    echo "  sigma/m_sim = $SIGMA cm^2/g"
    echo "  r_diss      = $RDISS"
    echo "  t_end       = 1.0 code units (matching working reference)"
    echo "=========================================="

    # --- Config.sh (matching working reference's setup) ---
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

    # --- Clean and build ---
    cd $SOURCE
    rm -f Gadget4_$NAME build/*.o build/*.mod 2>/dev/null || true
    mkdir -p build

    echo "--- Building (timeout 600s) ---"
    if timeout 600 make build -j 4 EXEC=Gadget4_$NAME > /tmp/build_${NAME}.log 2>&1; then
        echo "Build OK"
        ls -la $SOURCE/Gadget4_$NAME
    else
        echo "Build FAILED (exit $?)"
        tail -30 /tmp/build_${NAME}.log
        return 1
    fi

    # --- params.txt (matching working reference EXACTLY, except sigma/m via Config.sh) ---
    cat > $BASE/$NAME/params.txt <<EOF
%paramfile for rescaling-symmetry N-body verification
% Point: $NAME (simulation space, using working reference IC)
% SIDM sigma/m_sim = $SIGMA cm^2/g (via Config.sh)
% SIDM r_diss = $RDISS (via Config.sh)
% TimeMax = 1.0 code units (matching working reference)
% IC: cored NFW from working reference (r_s=3.6, rho_0=7.09e-3)

InitCondFile        ic.dat
OutputDir           output
SnapshotFileBase    snapshot
OutputListFilename  output_list.txt

ICFormat            1
SnapFormat          1

TimeLimitCPU              21000
CpuTimeBetRestartFile     3600
MaxMemSize                4000

TimeBegin           0.0
TimeMax             1.0

ComovingIntegrationOn    0
Omega0              0
OmegaLambda         0
OmegaBaryon         0
HubbleParam         1.0
Hubble              0
BoxSize             0

OutputListOn              1
TimeBetSnapshot           0.5
TimeOfFirstSnapshot       0.0
TimeBetStatistics         0.05
NumFilesPerSnapshot       1
MaxFilesWithConcurrentIO  1

ErrTolIntAccuracy        0.02
CourantFac               0.15
MaxSizeTimestep          0.002
MinSizeTimestep          0.0

TypeOfOpeningCriterion                1
ErrTolTheta                           0.7
ErrTolThetaMax                        0.9
ErrTolForceAcc                        0.0025
TopNodeFactor                         2.0
ActivePartFracForNewDomainDecomp      0.02

DesNumNgb                             64
MaxNumNgbDeviation                    8

UnitLength_in_cm         3.0856775814913673e21
UnitMass_in_g            1.989e43
UnitVelocity_in_cm_per_s 1.0e5
GravityConstantInternal  0

SofteningComovingClass0       0.1
SofteningMaxPhysClass0        0.1
SofteningComovingClass1       0.1
SofteningMaxPhysClass1        0.1
SofteningComovingClass2       0.1
SofteningMaxPhysClass2        0.1
SofteningClassOfPartType0     0
SofteningClassOfPartType1     1
SofteningClassOfPartType2     2

ArtBulkViscConst       1.0
MinEgySpec             0
InitGasTemp            0
EOF

    # --- output_list.txt: snapshots at t=0, 0.5, 1.0 ---
    cat > $BASE/$NAME/output_list.txt <<EOF
0.0
0.5
1.0
EOF

    # --- submit.sh (16 ranks matching working reference) ---
    cat > $BASE/$NAME/submit.sh <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=${NAME}
#SBATCH --partition=amd_512
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --mem=48G
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

mpirun -np 16 ~/dsidm_project/source/Gadget4_$NAME params.txt 2>&1 | tee run.log
rc=\$?
echo "=== Finished $NAME at \$(date) with exit code \$rc ==="
ls -la output/ 2>&1 || true
exit 0
EOF
    chmod +x $BASE/$NAME/submit.sh

    # --- Clean previous output ---
    rm -rf $BASE/$NAME/output $BASE/$NAME/restart \
           $BASE/$NAME/slurm_*.out $BASE/$NAME/slurm_*.err \
           $BASE/$NAME/run.log 2>/dev/null || true
    mkdir -p $BASE/$NAME/output

    echo "  $NAME setup complete"
}

# =====================================================================
# Build all 3 points (same sigma values, just different builds)
# =====================================================================
# P1: sigma/m=3.33 (elastic, r_diss=1.0) - slightly higher than ref's 3.0
# P2: sigma/m=0.167 (low, r_diss=1.05) - much lower, minimal scattering
# P3: sigma/m=7.33 (high, r_diss=1.05) - much higher, strong scattering

build_point P1_elastic_control 3.33  1.0
build_point P2_m3_low_sigma    0.167 1.05
build_point P3_m3_high_sigma   7.33  1.05

# --- Restore Config.sh ---
cp $SOURCE/Config.sh.bak_ref $SOURCE/Config.sh 2>/dev/null || true
echo "Config.sh restored"

# =====================================================================
# Submit all jobs
# =====================================================================
echo ""
echo "=========================================="
echo "Submitting jobs"
echo "=========================================="
for name in P1_elastic_control P2_m3_low_sigma P3_m3_high_sigma; do
    cd $BASE/$name
    out=$(sbatch submit.sh)
    echo "  $name: $out"
done

echo ""
echo "All submitted. Monitor with: squeue -u scg7816"
echo "Expected runtime: ~4 hours each (matching working reference)"
