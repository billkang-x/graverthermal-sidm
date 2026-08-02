#!/usr/bin/env python3
"""Comprehensive fix for the N-body verification jobs.

Bugs found from slurm logs:
1. IC Header has NumPart_ThisFile=[100000,0,...] but particles are in PartType1 group
   -> Header should be [0, 100000, 0, 0, 0, 0] and MassTable [0, 624.67, ...]
2. params.txt has invalid tags: RestartDir, ErrTolForce, MaxRMSDisplacementFac,
   CpuTreeDomainUpdate, UnitLuminosity_in_erg_s, UnitEnergy_in_ergs, MinGasTemp,
   PeriodicBoundariesOn
3. params.txt is missing required: TimeLimitCPU
4. OutputListOn=0 + empty output_list.txt: use OutputListOn=1 with proper list
   (or just rely on TimeBetSnapshot - but reference job uses OutputListOn=1)
5. Need SnapshotFileBase='snapshot' to match what analyze_nbody.py expects

Strategy: patch the local IC file header, regenerate params.txt from the working
reference template, upload everything, resubmit.
"""
import h5py
import numpy as np
import paramiko
import os
import time
import sys

HOST = "ssh.cn-zhongwei-1.paracloud.com"
PORT = 22
USER = "scg7816@ZC-M6"
PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"

LOCAL_IC = "D:/graverthermal-sidm/data/P5_nbody_verify/ics/ic.hdf5"
LOCAL_IC_FIXED = "D:/graverthermal-sidm/data/P5_nbody_verify/ics/ic_fixed.hdf5"

POINTS = [
    # (name, sigma_m_compile, r_diss_compile, t_end_gyr, description)
    ("P1_elastic_control", 0.1,   1.0,  0.68,  "Elastic control: no dissipation"),
    ("P2_m3_low_sigma",    0.005, 1.05, 0.07,  "M3 low sigma/m: weak dissipative cooling"),
    ("P3_m3_high_sigma",   0.220, 1.05, 0.10,  "M3 high sigma/m: strong cooling"),
]


def fix_ic_header():
    """Patch the IC file Header to correctly report particles in PartType1."""
    print("=== Fixing IC header ===")
    # Read all data from original
    with h5py.File(LOCAL_IC, 'r') as fin:
        header_attrs = dict(fin['Header'].attrs)
        pt1_data = {}
        for k in fin['PartType1'].keys():
            pt1_data[k] = fin['PartType1'][k][:]

    print(f"  Original NumPart_ThisFile: {header_attrs['NumPart_ThisFile']}")
    print(f"  Original MassTable: {header_attrs['MassTable']}")
    print(f"  PartType1 group has {pt1_data['Coordinates'].shape[0]} particles")

    # Fix header: particles are in PartType1, not PartType0
    n = pt1_data['Coordinates'].shape[0]
    header_attrs['NumPart_ThisFile'] = np.array([0, n, 0, 0, 0, 0], dtype=np.uint32)
    header_attrs['NumPart_Total'] = np.array([0, n, 0, 0, 0, 0], dtype=np.uint32)
    header_attrs['NumPart_Total_HighWord'] = np.array([0, 0, 0, 0, 0, 0], dtype=np.uint32)
    mt = np.zeros(6)
    mt[1] = float(pt1_data['Masses'][0]) if 'Masses' in pt1_data else 0.0
    header_attrs['MassTable'] = mt

    # Remove InternalEnergy from PartType1 (DM doesn't have internal energy)
    if 'InternalEnergy' in pt1_data:
        del pt1_data['InternalEnergy']
        print("  Removed InternalEnergy from PartType1 (not applicable to DM)")

    # Write fixed file
    with h5py.File(LOCAL_IC_FIXED, 'w') as fout:
        h_header = fout.create_group('Header')
        for k, v in header_attrs.items():
            h_header.attrs[k] = v
        h_pt1 = fout.create_group('PartType1')
        for k, v in pt1_data.items():
            h_pt1.create_dataset(k, data=v)

    print(f"  Fixed NumPart_ThisFile: {header_attrs['NumPart_ThisFile']}")
    print(f"  Fixed MassTable: {header_attrs['MassTable']}")
    print(f"  Wrote: {LOCAL_IC_FIXED}")

    # Verify
    with h5py.File(LOCAL_IC_FIXED, 'r') as f:
        print(f"  Verify: groups={list(f.keys())}")
        print(f"  Verify: NumPart_ThisFile={f['Header'].attrs['NumPart_ThisFile']}")
        print(f"  Verify: MassTable={f['Header'].attrs['MassTable']}")
        if 'PartType1' in f:
            print(f"  Verify: PartType1 has {list(f['PartType1'].keys())}")


def gen_params_content(name, sigma_m, r_diss, t_end, description):
    """Generate a Gadget4 params.txt using the working reference template.

    Based on /public3/home/scg7816/sidm-diskbar-gate1/recheck/.../param.txt
    Adapted for our SIDM dissipative halo: kpc/Msun/km-s units, isolated, N=1e5.
    """
    return f"""%paramfile for N-body verification
% Point: {name}
% Description: {description}
% SIDM sigma/m = {sigma_m} cm^2/g (compile-time via Config.sh)
% SIDM r_diss = {r_diss} (compile-time via Config.sh)
% Target evolution time: {t_end} Gyr

InitCondFile        ic
OutputDir           output
SnapshotFileBase    snapshot
OutputListFilename  output_list.txt

ICFormat            3
SnapFormat          3

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
SofteningMaxPhysClass1        0.0005
SofteningComovingClass2       0.0005
SofteningMaxPhysClass2        0.0005
SofteningClassOfPartType0     0
SofteningClassOfPartType1     1
SofteningClassOfPartType2     2

ArtBulkViscConst       1.0
MinEgySpec             0
InitGasTemp            0
"""


def gen_output_list(t_end):
    """Generate output_list.txt with snapshot times.

    We want snapshots at t=0 (IC), t_end/2, and t_end (final).
    Format: one floating-point time per line.
    """
    # Just use t_end as the single output (we set TimeBetSnapshot=t_end too)
    # Format is just one float per line
    return f"{t_end:.6f}\n"


def gen_submit_sh(name):
    """Generate submit.sh using the working reference template pattern."""
    return f"""#!/usr/bin/env bash
#SBATCH --job-name={name}
#SBATCH --partition=amd_512
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=24G
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err

set -euo pipefail

source /public3/soft/modules/module.sh
module purge
module load gcc/10.2.0 gsl/2.0 fftw/3.3.8-fjy mpi/openmpi/4.1.1-gcc7.3.0 hdf5/1.8.13-gcc-zyq

cd "$SLURM_SUBMIT_DIR"
mkdir -p output

hdf5_root=/public3/soft/hdf5/1.8.13-gcc-zyq
gsl_root=/public3/soft/gsl/gsl2.0
fftw_root=/public3/soft/fftw/3.3.8-fjy
export LD_LIBRARY_PATH="$hdf5_root/lib:$gsl_root/lib:$fftw_root/lib:${{LD_LIBRARY_PATH:-}}"
export OMPI_MCA_btl=self,vader,tcp

echo "=== Starting {name} at $(date) ==="
echo "Job ID: $SLURM_JOB_ID on $(hostname)"
echo "Mem avail: $(free -h | head -2)"

~/dsidm_project/source/Gadget4_{name} params.txt 2>&1 | tee run.log

echo "=== Finished {name} at $(date) ==="
ls -la output/ 2>&1 || true
"""


def run(ssh, cmd, timeout=60):
    print(f"$ {cmd[:200]}", flush=True)
    stdin, stdout, stderr = ssh.exec_command(cmd, timeout=timeout)
    out = stdout.read().decode(errors='replace').rstrip()
    err = stderr.read().decode(errors='replace').rstrip()
    if out:
        print(out[:3000], flush=True)
    if err:
        print(f"  [stderr] {err[:2000]}", flush=True)
    return out, err


def main():
    # Step 1: Fix the IC file locally
    fix_ic_header()
    print()

    # Step 2: Connect to HPC
    print("=== Connecting to HPC ===")
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(HOST, port=PORT, username=USER, password=PASS,
                timeout=30, look_for_keys=False, allow_agent=False)
    print("Connected.\n")
    sftp = ssh.open_sftp()
    base_remote = "/public3/home/scg7816/dsidm_project/nbody_verify"

    # Step 3: Upload fixed IC and new params/submit scripts for each point
    print("=== Uploading fixed IC, params, submit scripts ===")
    for name, sigma_m, r_diss, t_end, desc in POINTS:
        remote_dir = f"{base_remote}/{name}"

        # Upload fixed IC
        sftp.put(LOCAL_IC_FIXED, f"{remote_dir}/ic.hdf5")
        print(f"  {name}: ic.hdf5 uploaded (fixed header)")

        # Write and upload params.txt
        params_content = gen_params_content(name, sigma_m, r_diss, t_end, desc)
        local_params = f"D:/graverthermal-sidm/src/P5_gadget/ic_scripts/{name}_params_v3.txt"
        with open(local_params, 'w') as f:
            f.write(params_content)
        sftp.put(local_params, f"{remote_dir}/params.txt")
        print(f"  {name}: params.txt uploaded ({len(params_content)} chars)")

        # Write and upload output_list.txt
        output_list_content = gen_output_list(t_end)
        with sftp.open(f"{remote_dir}/output_list.txt", 'w') as f:
            f.write(output_list_content)
        print(f"  {name}: output_list.txt uploaded (t_end={t_end})")

        # Write and upload submit.sh
        submit_content = gen_submit_sh(name)
        local_submit = f"D:/graverthermal-sidm/src/P5_gadget/ic_scripts/{name}_submit_v2.sh"
        with open(local_submit, 'w') as f:
            f.write(submit_content)
        sftp.put(local_submit, f"{remote_dir}/submit.sh")
        sftp.chmod(f"{remote_dir}/submit.sh", 0o755)
        print(f"  {name}: submit.sh uploaded")

        # Clean up previous output
        run(ssh, f"rm -rf {remote_dir}/output {remote_dir}/restart {remote_dir}/slurm_*.out {remote_dir}/slurm_*.err {remote_dir}/run.log 2>&1; mkdir -p {remote_dir}/output")
        print()

    sftp.close()

    # Step 4: Verify executables still exist
    print("=== Verify executables ===")
    for name, _, _, _, _ in POINTS:
        run(ssh, f"ls -la ~/dsidm_project/source/Gadget4_{name}")
    print()

    # Step 5: Submit jobs
    print("=== Submitting jobs ===")
    job_ids = []
    for name, _, _, _, _ in POINTS:
        remote_dir = f"{base_remote}/{name}"
        out, _ = run(ssh, f"cd {remote_dir} && sbatch submit.sh")
        # Parse "Submitted batch job 40886445"
        for line in out.split('\n'):
            if 'Submitted batch job' in line:
                jid = line.split()[-1]
                job_ids.append((name, jid))
                print(f"  {name}: job {jid}")
                break
    print()

    # Step 6: Poll for completion
    print("=== Polling for completion ===")
    start = time.time()
    max_wait = 1800  # 30 minutes
    while True:
        elapsed = time.time() - start
        if elapsed > max_wait:
            print(f"  Timeout after {max_wait/60:.0f} min")
            break

        out, _ = run(ssh, "squeue -u scg7816 --noheader 2>&1")
        # Count only our jobs (filter out sidm-gat 40880935)
        our_jobs_running = 0
        for line in out.split('\n'):
            line = line.strip()
            if not line:
                continue
            for name, jid in job_ids:
                if jid in line:
                    our_jobs_running += 1
                    break
        print(f"  [{elapsed/60:.0f}m {elapsed%60:.0f}s] {our_jobs_running}/{len(job_ids)} of our jobs still in queue")

        if our_jobs_running == 0:
            print("  All jobs finished!")
            break

        time.sleep(30)

    # Step 7: Check results
    print("\n=== Check results ===")
    for name, jid in job_ids:
        remote_dir = f"{base_remote}/{name}"
        print(f"\n--- {name} (job {jid}) ---")
        run(ssh, f"sacct -j {jid} --format=JobID,State,Elapsed,ExitCode,MaxRSS --noheader")
        run(ssh, f"ls -la {remote_dir}/output/ 2>&1 | head -20")
        # Tail the slurm log
        run(ssh, f"tail -30 {remote_dir}/slurm_{jid}.out 2>&1")
        if os.path.exists(f"D:/graverthermal-sidm/data/P5_nbody_verify/{name}"):
            pass  # local dir exists

    # Step 8: Download snapshots if any
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
