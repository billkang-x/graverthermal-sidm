"""
Submit N-body verification jobs for the 3 velocity-INDEPENDENT test points.

The Gadget4_dsidm module currently uses compile-time constants SIDM_SIGMA_OVER_MASS
and SIDM_R_DISS from Config.sh. Velocity-dependent mode (table loading) is not
yet implemented. Therefore, we run only the 3 points that use constant sigma/m
and constant r_diss:

  - P1_elastic_control: sigma/m=0.1, r_diss=1.0 (elastic, no dissipation)
  - P2_m3_low_sigma:    sigma/m=0.005, r_diss=1.05
  - P3_m3_high_sigma:   sigma/m=0.220, r_diss=1.05

For each point, we need to:
  1. Modify Config.sh with the appropriate sigma/m and r_diss values
  2. Recompile Gadget4_dsidm
  3. Submit the Slurm job

Since the elastic control (r_diss=1.0) is just the standard SIDM module, we
could potentially use the existing non-dissipative binary. But for a clean test,
we recompile all three with SIDM_DISSIPATIVE enabled (so the same code path is
exercised).

For P2 and P3 (same r_diss=1.05 but different sigma/m), we can actually submit
both with the same binary if sigma/m is also read from params.txt at runtime.
But the current code uses SIDM_SIGMA_OVER_MASS compile-time constant. So we
need to either recompile twice or modify the code.

Approach: Use a small wrapper that creates a per-point Config.sh, rebuilds with
a unique EXEC name, and submits the job. Three rebuilds ~ 1-2 min each = ~5 min
total.

Actually, simpler: we can use the BEGrun() function in Gadget4 to read sigma/m
from the paramfile at runtime. But sidm_dissipative.cc initializes from compile-
time constants in its constructor. We need to add runtime override.

For minimal effort, let's just recompile 3 times with different Config.sh.
"""
import paramiko
import os
import time
import re

HOST = "ssh.cn-zhongwei-1.paracloud.com"
PORT = 22
USER = "scg7816@ZC-M6"
PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"

# Each point: (name, sigma_m, r_diss)
POINTS = [
    ("P1_elastic_control", 0.1,   1.0),
    ("P2_m3_low_sigma",    0.005, 1.05),
    ("P3_m3_high_sigma",   0.220, 1.05),
]

# Submit script template (per point)
SUBMIT_TEMPLATE = """#!/bin/bash
#SBATCH --job-name={name}
#SBATCH --partition=amd_512
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err

source ~/env.sh

cd $SLURM_SUBMIT_DIR
echo "=== Starting {name} at $(date) ==="
echo "Job ID: $SLURM_JOB_ID on $(hostname)"

mpirun -np 8 ~/dsidm_project/source/Gadget4_{name} params.txt

echo "=== Finished {name} at $(date) ==="
ls -la output/
"""

# Config.sh template (rebuilds with per-point constants)
CONFIG_TEMPLATE = """SELFGRAVITY
NTYPES=3
EVALPOTENTIAL
PRESSURE_ENTROPY_SPH
COOLING
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


def run(client, cmd, timeout=300):
    print(f"$ {cmd[:200]}", flush=True)
    stdin, stdout, stderr = client.exec_command(cmd, timeout=timeout)
    try:
        out = stdout.read().decode(errors='replace').rstrip()
    except Exception:
        out = ""
    try:
        err = stderr.read().decode(errors='replace').rstrip()
    except Exception:
        err = ""
    if out:
        print(out[:5000], flush=True)
    if err:
        print(f"  [stderr] {err[:3000]}", flush=True)
    return out, err


def write_remote_file(client, remote_path, content):
    """Write a file to the remote path. Expand ~ to home dir."""
    if remote_path.startswith('~'):
        # SFTP doesn't expand ~ — use absolute path
        remote_path = '/public3/home/scg7816' + remote_path[1:]
    sftp = client.open_sftp()
    with sftp.open(remote_path, 'w') as f:
        f.write(content)
    sftp.close()


def main():
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(HOST, port=PORT, username=USER, password=PASS,
                   timeout=30, look_for_keys=False, allow_agent=False)
    print("Connected.\n", flush=True)

    submitted = []

    for name, sigma_m, r_diss in POINTS:
        print(f"\n{'=' * 60}")
        print(f"=== Point: {name} (sigma/m={sigma_m}, r_diss={r_diss}) ===")
        print(f"{'=' * 60}")

        # 1. Write Config.sh
        print("\n--- 1. Write Config.sh ---", flush=True)
        config_content = CONFIG_TEMPLATE.format(sigma_m=sigma_m, r_diss=r_diss)
        write_remote_file(client, "~/dsidm_project/source/Config.sh", config_content)
        run(client, "cat ~/dsidm_project/source/Config.sh")

        # 2. Recompile (use nohup + log file approach)
        print("\n--- 2. Recompile ---", flush=True)
        exec_name = f"Gadget4_{name}"
        log_file = f"~/dsidm_project/build_{name}.log"
        build_cmd = (
            "cd ~/dsidm_project/source && "
            "source ~/env.sh > /dev/null 2>&1 && "
            "make clean > /dev/null 2>&1 && "
            f"make build -j 8 EXEC={exec_name}"
        )
        full_cmd = f"nohup bash -c '{build_cmd}' > {log_file} 2>&1 & echo $!"
        out, _ = run(client, full_cmd)
        pid = out.strip().split('\n')[-1].strip() if out else ""
        print(f"  Build PID: {pid}", flush=True)

        # Poll for completion
        max_wait_s = 300
        check_interval = 15
        waited = 0
        while waited < max_wait_s:
            out, _ = run(client, f"if kill -0 {pid} 2>/dev/null; then echo RUNNING; else echo DONE; fi")
            status = out.strip().split('\n')[-1].strip() if out else "UNKNOWN"
            print(f"  [{waited}s] status: {status}", flush=True)
            if status == "DONE":
                break
            time.sleep(check_interval)
            waited += check_interval

        # Check build log tail
        run(client, f"tail -10 {log_file}")

        # Check executable
        run(client, f"ls -la ~/dsidm_project/source/{exec_name} 2>&1")

        # 3. Write submit script
        print(f"\n--- 3. Write submit script ---", flush=True)
        submit_content = SUBMIT_TEMPLATE.format(name=name)
        point_dir = f"~/dsidm_project/nbody_verify/{name}"
        write_remote_file(client, f"{point_dir}/submit.sh", submit_content)
        run(client, f"chmod +x {point_dir}/submit.sh")

        # 4. Make output dir and submit job
        print(f"\n--- 4. Submit Slurm job ---", flush=True)
        run(client, f"mkdir -p {point_dir}/output")
        out, _ = run(client, f"cd {point_dir} && sbatch submit.sh")
        # Parse job ID
        for line in out.split('\n'):
            if 'Submitted batch job' in line:
                job_id = line.split()[-1]
                submitted.append((name, job_id))
                print(f"  {name}: job {job_id}", flush=True)
                break

    print(f"\n\nSubmitted {len(submitted)} jobs:", flush=True)
    for name, jid in submitted:
        print(f"  {name}: {jid}", flush=True)

    # 5. Wait for jobs to finish
    print("\n=== 5. Poll for job completion ===", flush=True)
    max_wait_min = 60
    poll_interval = 30
    waited = 0
    while waited < max_wait_min * 60:
        out, _ = run(client, "squeue -u scg7816 --noheader 2>&1")
        our_jobs_running = sum(1 for _, jid in submitted if jid in out)
        print(f"  [{waited//60}m {waited%60}s] {our_jobs_running}/{len(submitted)} jobs still running",
              flush=True)
        if our_jobs_running == 0:
            print("  All jobs finished!", flush=True)
            break
        time.sleep(poll_interval)
        waited += poll_interval

    # 6. Check job status and download snapshots
    print("\n=== 6. Check job results ===", flush=True)
    for name, jid in submitted:
        point_dir = f"~/dsidm_project/nbody_verify/{name}"
        run(client, f"ls -la {point_dir}/output/ 2>&1 | head -10")
        run(client, f"tail -15 {point_dir}/slurm_{jid}.out 2>&1")

    # 7. Download snapshots
    print("\n=== 7. Download snapshots ===", flush=True)
    sftp = client.open_sftp()
    base_remote = "/public3/home/scg7816/dsidm_project/nbody_verify"
    local_base = "D:/graverthermal-sidm/data/P5_nbody_verify"
    for name, jid in submitted:
        point_local_dir = os.path.join(local_base, name)
        os.makedirs(point_local_dir, exist_ok=True)
        point_remote_dir = f"{base_remote}/{name}/output"

        try:
            files = sftp.listdir(point_remote_dir)
        except Exception as e:
            print(f"  {name}: cannot list remote output dir: {e}", flush=True)
            continue

        snap_files = sorted([f for f in files if f.startswith('snap_') and f.endswith('.hdf5')])
        if not snap_files:
            print(f"  {name}: no snapshot files found", flush=True)
            continue

        # Download the last (final) snapshot
        last_snap = snap_files[-1]
        remote_path = f"{point_remote_dir}/{last_snap}"
        local_path = os.path.join(point_local_dir, last_snap)
        print(f"  {name}: downloading {last_snap}...", flush=True)
        try:
            sftp.get(remote_path, local_path)
            print(f"    saved to {local_path} ({os.path.getsize(local_path)/1e6:.1f} MB)",
                  flush=True)
        except Exception as e:
            print(f"    ERROR: {e}", flush=True)

        # Also download slurm log
        try:
            sftp.get(f"{base_remote}/{name}/slurm_{jid}.out",
                     os.path.join(point_local_dir, f"slurm_{jid}.out"))
        except Exception:
            pass

    sftp.close()
    client.close()
    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
