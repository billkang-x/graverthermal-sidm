#!/bin/bash
# =====================================================================
# Plan 1: Rescaling symmetry N-body verification
# Build & submit 3 Gadget4 jobs in SIMULATION SPACE
#
# Sim-space halo: r_s=3.6 kpc, rho_0=7.09e-3 Msun/pc^3
#   (matches working reference HPC run that completed successfully)
#
# Rescaling: lambda=0.023611, mu=0.018565
#   sigma_sim = sigma_phys / 0.030028 = sigma_phys * 33.30
#   t_sim     = t_phys     / 0.026627 = t_phys * 37.56
#
# Matching fluid model ACTUAL evolution times:
#   P1: sigma_sim=3.33,   t_sim=8.75  Gyr = 8.95 code units
#   P2: sigma_sim=0.167,  t_sim=1.50  Gyr = 1.54 code units
#   P3: sigma_sim=7.33,   t_sim=3.76  Gyr = 3.84 code units
#
# Projected mass ratio M_2D(R_inner)/M_2D(R_outer) is INVARIANT under rescaling.
#   R_inner_sim = 0.847 kpc  (= 20 pc / lambda)
#   R_outer_sim = 3.812 kpc  (= 90 pc / lambda)
# Compare directly with fluid model predictions.
# =====================================================================

set -e

SOURCE=~/dsidm_project/source
BASE=~/dsidm_project/nbody_verify_sim
IC_DIR=~/dsidm_project/nbody_verify_sim/ic

# --- Cancel any previous jobs ---
scancel --user=scg7816 --state=RUNNING,PENDING 2>/dev/null || true
echo "Cancelled previous jobs"

# --- Load build environment ---
source /public1/soft/modules/module.sh
module purge
module load gcc/10.2.0 gsl/2.0 fftw/3.3.8-fjy mpi/openmpi/4.1.1-gcc7.3.0 hdf5/1.8.13-gcc-zyq 2>&1 || true
echo "Modules loaded"

# --- Create directories ---
mkdir -p $BASE/P1_elastic_control/output
mkdir -p $BASE/P2_m3_low_sigma/output
mkdir -p $BASE/P3_m3_high_sigma/output
mkdir -p $IC_DIR

# --- Upload IC (ic.dat should already be in $IC_DIR from upload step) ---
if [ ! -f "$IC_DIR/ic.dat" ]; then
    echo "ERROR: $IC_DIR/ic.dat not found!"
    echo "Please upload the simulation-space IC first."
    exit 1
fi
echo "IC found: $(ls -la $IC_DIR/ic.dat)"

# --- Backup Config.sh ---
cp $SOURCE/Config.sh $SOURCE/Config.sh.bak_sim 2>/dev/null || true

# =====================================================================
# Build function
# =====================================================================
build_point() {
    local NAME=$1
    local SIGMA=$2
    local RDISS=$3
    local TEND_CODE=$4   # in code units (t_sim / 0.978 Gyr)

    echo ""
    echo "=========================================="
    echo "Building $NAME"
    echo "  sigma/m_sim = $SIGMA cm^2/g"
    echo "  r_diss      = $RDISS"
    echo "  t_end       = $TEND_CODE code units"
    echo "=========================================="

    # --- Config.sh (matching working reference, no SPH/COOLING) ---
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

    # --- Copy IC to point directory ---
    cp $IC_DIR/ic.dat $BASE/$NAME/ic.dat

    # --- params.txt (matching working reference) ---
    # Using Softening=0.1 kpc (matching working reference that succeeded)
    # TimeMax in code units (1 code unit = 0.978 Gyr)
    cat > $BASE/$NAME/params.txt <<EOF
%paramfile for rescaling-symmetry N-body verification
% Point: $NAME (simulation space)
% SIDM sigma/m_sim = $SIGMA cm^2/g
% SIDM r_diss = $RDISS
% TimeMax = $TEND_CODE code units (sim space)

InitCondFile        ic.dat
OutputDir           output
SnapshotFileBase    snapshot
OutputListFilename  output_list.txt

ICFormat            1
SnapFormat          1

TimeLimitCPU              14400
CpuTimeBetRestartFile     10800
MaxMemSize                2000

TimeBegin           0.0
TimeMax             $TEND_CODE

ComovingIntegrationOn    0
Omega0              0
OmegaLambda         0
OmegaBaryon         0
HubbleParam         1.0
Hubble              0
BoxSize             0

OutputListOn              1
TimeBetSnapshot           $TEND_CODE
TimeOfFirstSnapshot       $TEND_CODE
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

% Softening = 0.1 kpc (matches working reference)
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

    # --- output_list.txt: snapshots at t=0, t_end/2, t_end ---
    # This gives us intermediate states to track evolution
    echo "0.0" > $BASE/$NAME/output_list.txt
    echo "$TEND_CODE" >> $BASE/$NAME/output_list.txt

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

# If we have restart files, submit continuation
if [ -f restart/restart.idx ]; then
    echo "Restart files found, submitting continuation..."
    sbatch submit_restart.sh
fi

exit 0
EOF
    chmod +x $BASE/$NAME/submit.sh

    # --- submit_restart.sh (for continuation if needed) ---
    cat > $BASE/$NAME/submit_restart.sh <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=${NAME}_r
#SBATCH --partition=amd_512
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --mem=48G
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --dependency=afterany:\$SLURM_JOB_ID

cd "\$SLURM_SUBMIT_DIR"
source /public1/soft/modules/module.sh
module purge
module load gcc/10.2.0 gsl/2.0 fftw/3.3.8-fjy mpi/openmpi/4.1.1-gcc7.3.0 hdf5/1.8.13-gcc-zyq 2>&1 || true

hdf5_root=/public3/soft/hdf5/1.8.13-gcc-zyq
gsl_root=/public3/soft/gsl/gsl2.0
fftw_root=/public3/soft/fftw/3.3.8-fjy
export LD_LIBRARY_PATH="\$hdf5_root/lib:\$gsl_root/lib:\$fftw_root/lib:\${LD_LIBRARY_PATH:-}"
export OMPI_MCA_btl=self,vader,tcp

echo "=== Restarting $NAME at \$(date) ==="
mpirun -np 16 ~/dsidm_project/source/Gadget4_$NAME params.txt 2>&1 | tee -a run.log
rc=\$?
echo "=== Restart finished at \$(date) with exit code \$rc ==="
ls -la output/ 2>&1 || true

if [ -f restart/restart.idx ]; then
    sbatch submit_restart.sh
fi
exit 0
EOF
    chmod +x $BASE/$NAME/submit_restart.sh

    # --- Clean previous output ---
    rm -rf $BASE/$NAME/output $BASE/$NAME/restart \
           $BASE/$NAME/slurm_*.out $BASE/$NAME/slurm_*.err \
           $BASE/$NAME/run.log 2>/dev/null || true
    mkdir -p $BASE/$NAME/output

    echo "  $NAME setup complete"
}

# =====================================================================
# Build all 3 points
# =====================================================================
# Matching fluid model ACTUAL evolution times:
#   P1: sigma_phys=0.1,   sigma_sim=3.33,   t_phys=0.233 -> t_sim=8.95 code
#   P2: sigma_phys=0.005, sigma_sim=0.167,  t_phys=0.040 -> t_sim=1.54 code
#   P3: sigma_phys=0.22,  sigma_sim=7.33,   t_phys=0.100 -> t_sim=3.84 code

build_point P1_elastic_control 3.33   1.0  8.95
build_point P2_m3_low_sigma    0.167  1.05 1.54
build_point P3_m3_high_sigma   7.33   1.05 3.84

# --- Restore Config.sh ---
cp $SOURCE/Config.sh.bak_sim $SOURCE/Config.sh 2>/dev/null || true
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
echo "Check progress: tail -f ~/dsidm_project/nbody_verify_sim/*/run.log"
