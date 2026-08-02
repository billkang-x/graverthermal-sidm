#!/usr/bin/env python3
"""Rebuild Gadget4 without PRESSURE_ENTROPY_SPH and COOLING for pure-DM SIDM runs.

The crash was in sph::init_entropy() because we have 0 gas particles but
PRESSURE_ENTROPY_SPH is enabled, causing a NULL pointer deref. Since these
are pure dark matter runs with SIDM, we don't need SPH or gas cooling.

For each test point (different SIDM_SIGMA_OVER_MASS, SIDM_R_DISS), we need
to rebuild with the appropriate compile-time constants. The SIDM module
reads these constants at compile time.

Strategy:
1. Update Config.sh to remove PRESSURE_ENTROPY_SPH and COOLING
2. For each point: set SIDM_SIGMA_OVER_MASS and SIDM_R_DISS, run `make build`
3. Submit jobs with the new executables
"""
import paramiko
import os
import time

HOST = "ssh.cn-zhongwei-1.paracloud.com"
PORT = 22
USER = "scg7816@ZC-M6"
PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"

# (name, sigma_m, r_diss, t_end_gyr)
POINTS = [
    ("P1_elastic_control", 0.1,   1.0,  0.68),
    ("P2_m3_low_sigma",    0.005, 1.05, 0.07),
    ("P3_m3_high_sigma",   0.220, 1.05, 0.10),
]


def run(ssh, cmd, timeout=120):
    print(f"$ {cmd[:200]}", flush=True)
    stdin, stdout, stderr = ssh.exec_command(cmd, timeout=timeout)
    out = stdout.read().decode(errors='replace').rstrip()
    err = stderr.read().decode(errors='replace').rstrip()
    if out:
        print(out[:5000], flush=True)
    if err:
        print(f"  [stderr] {err[:2000]}", flush=True)
    return out, err


def gen_config_sh(sigma_m, r_diss):
    """Generate Config.sh for a specific test point."""
    return f"""SELFGRAVITY
NTYPES=3
EVALPOTENTIAL
SIDM
SIDM_SIGMA_OVER_MASS={sigma_m}
SIDM_K_NEIGHBORS=32
SIDM_PMAX=0.1
DOUBLEPRECISION=1
GADGET2_HEADER

# Dissipative fSIDM extension
SIDM_DISSIPATIVE
SIDM_R_DISS={r_diss}
"""


def gen_submit_sh(name):
    """Submit script: ntasks=2, mpirun -np 2."""
    return f"""#!/usr/bin/env bash
#SBATCH --job-name={name}
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
export LD_LIBRARY_PATH="$hdf5_root/lib:$gsl_root/lib:$fftw_root/lib:${{LD_LIBRARY_PATH:-}}"
export OMPI_MCA_btl=self,vader,tcp

echo "=== Starting {name} at $(date) ==="
echo "Job ID: $SLURM_JOB_ID on $(hostname)"
echo "MPI ranks: $SLURM_NTASKS"

mpirun -np 2 ~/dsidm_project/source/Gadget4_{name} params.txt 2>&1 | tee run.log
rc=$?
echo "=== Finished {name} at $(date) with exit code $rc ==="
ls -la output/ 2>&1 || true
exit 0
"""


def gen_params_content(name, sigma_m, r_diss, t_end, description):
    """Generate params.txt without TreecoolFile (no COOLING)."""
    return f"""%paramfile for N-body verification
% Point: {name}
% Description: {description}
% SIDM sigma/m = {sigma_m} cm^2/g (compile-time via Config.sh)
% SIDM r_diss = {r_diss} (compile-time via Config.sh)
% Target evolution time: {t_end} Gyr

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
TimeMax             {t_end}

ComovingIntegrationOn    0
Omega0              0
OmegaLambda         0
OmegaBaryon         0
HubbleParam         1.0
Hubble              0
BoxSize             0

OutputListOn              1
TimeBetSnapshot           {t_end}
TimeOfFirstSnapshot       {t_end}
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
"""


def main():
    print("=== Connecting to HPC ===")
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(HOST, port=PORT, username=USER, password=PASS,
                timeout=30, look_for_keys=False, allow_agent=False)
    print("Connected.\n")
    sftp = ssh.open_sftp()
    base_remote = "/public3/home/scg7816/dsidm_project/nbody_verify"
    source_dir = "/public3/home/scg7816/dsidm_project/source"

    # Cancel stale jobs
    print("=== Cancelling stale jobs ===")
    run(ssh, "squeue -u scg7816 --noheader")
    run(ssh, "scancel 40889625 40889626 40889627 2>&1 || true")
    print()

    # Step 1: Backup the current Config.sh
    print("=== Backing up Config.sh ===")
    run(ssh, f"cp {source_dir}/Config.sh {source_dir}/Config.sh.bak_sph")
    print()

    # Step 2: For each point, write Config.sh, build, and prepare submission
    print("=== Building executables (one per point) ===")
    for name, sigma_m, r_diss, t_end in POINTS:
        print(f"\n--- Building Gadget4_{name} (sigma/m={sigma_m}, r_diss={r_diss}) ---")
        # Write Config.sh
        config_content = gen_config_sh(sigma_m, r_diss)
        local_config = f"D:/graverthermal-sidm/src/P5_gadget/ic_scripts/{name}_Config.sh"
        with open(local_config, 'w', newline='\n') as f:
            f.write(config_content)
        sftp.put(local_config, f"{source_dir}/Config.sh")
        # Strip CR
        run(ssh, f"sed -i 's/\\r$//g' {source_dir}/Config.sh")
        run(ssh, f"cat {source_dir}/Config.sh")

        # Clean previous build artifacts
        run(ssh, f"cd {source_dir} && rm -f Gadget4 Exec Gadget4_{name} build/*.o build/*.mod 2>/dev/null; mkdir -p build", timeout=60)

        # Build with make build (bypasses check.py)
        print(f"  Building...")
        out, err = run(ssh, f"cd {source_dir} && source /public1/soft/modules/module.sh && module load gcc/10.2.0 gsl/2.0 fftw/3.3.8-fjy mpi/openmpi/4.1.1-gcc7.3.0 hdf5/1.8.13-gcc-zyq && make build -j 4 EXEC=Gadget4_{name} 2>&1 | tail -30", timeout=300)

        # Verify executable
        run(ssh, f"ls -la {source_dir}/Gadget4_{name} 2>&1")

        # Update params.txt and submit.sh for this point
        remote_dir = f"{base_remote}/{name}"
        params_content = gen_params_content(name, sigma_m, r_diss, t_end,
                                             f"sigma/m={sigma_m}, r_diss={r_diss}")
        local_params = f"D:/graverthermal-sidm/src/P5_gadget/ic_scripts/{name}_params_v5.txt"
        with open(local_params, 'w', newline='\n') as f:
            f.write(params_content)
        sftp.put(local_params, f"{remote_dir}/params.txt")
        run(ssh, f"sed -i 's/\\r$//g' {remote_dir}/params.txt")

        # Write output_list.txt
        with sftp.open(f"{remote_dir}/output_list.txt", 'w') as f:
            f.write(f"{t_end:.6f}\n")

        # Write submit.sh
        submit_content = gen_submit_sh(name)
        local_submit = f"D:/graverthermal-sidm/src/P5_gadget/ic_scripts/{name}_submit_v5.sh"
        with open(local_submit, 'w', newline='\n') as f:
            f.write(submit_content)
        sftp.put(local_submit, f"{remote_dir}/submit.sh")
        sftp.chmod(f"{remote_dir}/submit.sh", 0o755)
        run(ssh, f"sed -i 's/\\r$//g' {remote_dir}/submit.sh")

        # Clean output dir
        run(ssh, f"rm -rf {remote_dir}/output {remote_dir}/restart {remote_dir}/slurm_*.out {remote_dir}/slurm_*.err {remote_dir}/run.log 2>/dev/null; mkdir -p {remote_dir}/output")

    print()

    # Restore the original Config.sh backup (we modified it for each build)
    print("=== Restoring Config.sh backup ===")
    run(ssh, f"cp {source_dir}/Config.sh.bak_sph {source_dir}/Config.sh")

    sftp.close()

    # Step 3: Submit jobs
    print("\n=== Submitting jobs ===")
    job_ids = []
    for name, _, _, _ in POINTS:
        remote_dir = f"{base_remote}/{name}"
        out, _ = run(ssh, f"cd {remote_dir} && sbatch submit.sh")
        for line in out.split('\n'):
            if 'Submitted batch job' in line:
                jid = line.split()[-1]
                job_ids.append((name, jid))
                print(f"  {name}: job {jid}")
                break
    print()

    # Wait 3 minutes for jobs to start and reach integration phase
    print("=== Waiting 3 minutes ===")
    time.sleep(180)
    for name, jid in job_ids[:1]:
        remote_dir = f"{base_remote}/{name}"
        print(f"--- {name} ({jid}) early state ---")
        run(ssh, f"squeue -u scg7816 --noheader")
        run(ssh, f"tail -60 {remote_dir}/run.log 2>&1")
    print()

    # Poll
    print("=== Polling for completion ===")
    start = time.time()
    max_wait = 3600
    while True:
        elapsed = time.time() - start
        if elapsed > max_wait:
            print(f"  Timeout after {max_wait/60:.0f} min")
            break
        out, _ = run(ssh, "squeue -u scg7816 --noheader 2>&1")
        our_jobs_running = 0
        for line in out.split('\n'):
            line = line.strip()
            if not line:
                continue
            for name, jid in job_ids:
                if jid in line:
                    our_jobs_running += 1
                    break
        if job_ids:
            n, j = job_ids[0]
            out2, _ = run(ssh, f"wc -l {base_remote}/{n}/run.log 2>&1 || echo no-log")
            last_line = out2.splitlines()[-1] if out2 else ''
        else:
            last_line = ''
        print(f"  [{int(elapsed/60)}m {int(elapsed%60)}s] {our_jobs_running}/{len(job_ids)} running | {last_line}")
        if our_jobs_running == 0:
            print("  All jobs finished!")
            break
        time.sleep(60)

    # Check results
    print("\n=== Check results ===")
    for name, jid in job_ids:
        remote_dir = f"{base_remote}/{name}"
        print(f"\n--- {name} (job {jid}) ---")
        run(ssh, f"sacct -j {jid} --format=JobID,JobName%30,State,Elapsed,ExitCode,MaxRSS --noheader")
        run(ssh, f"ls -la {remote_dir}/output/ 2>&1 | head -30")
        run(ssh, f"tail -80 {remote_dir}/run.log 2>&1")

    # Download snapshots
    print("\n=== Downloading snapshots ===")
    sftp = ssh.open_sftp()
    for name, jid in job_ids:
        remote_dir = f"{base_remote}/{name}"
        local_dir = f"D:/graverthermal-sidm/data/P5_nbody_verify/{name}"
        os.makedirs(local_dir, exist_ok=True)
        try:
            files = sftp.listdir(f"{remote_dir}/output")
            print(f"  {name}: remote files in output/: {files}")
            for f in files:
                if f.startswith('snapshot') or f.startswith('snap'):
                    remote_path = f"{remote_dir}/output/{f}"
                    local_path = f"{local_dir}/{f}"
                    print(f"    Downloading {f}...")
                    sftp.get(remote_path, local_path)
                    print(f"      -> {local_path} ({os.path.getsize(local_path)} bytes)")
        except Exception as e:
            print(f"  {name}: error listing/downloading: {e}")
    sftp.close()

    ssh.close()
    print("\nDone.")


if __name__ == "__main__":
    main()
