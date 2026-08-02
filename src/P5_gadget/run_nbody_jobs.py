"""
Task #41 runner: Generate ICs and submit N-body verification jobs to the HPC.

Steps:
  1. For each of the 5 test points:
     a. Generate IC HDF5 (gen_ic.py)
     b. Submit Slurm job
  2. Wait for all jobs to finish (poll squeue)
  3. Download final snapshots back to local machine for analysis.
"""
import paramiko
import os
import time
import json

HOST = "ssh.cn-zhongwei-1.paracloud.com"
PORT = 22
USER = "scg7816@ZC-M6"
PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"

POINTS = [
    "P1_elastic_control",
    "P2_m3_low_sigma",
    "P3_m3_high_sigma",
    "P4_m1_low_sigma",
    "P5_m1_high_sigma",
]

LOCAL_DATA_DIR = "D:/graverthermal-sidm/data/P5_nbody_verify"


def run(client, cmd, timeout=120):
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


def main():
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(HOST, port=PORT, username=USER, password=PASS,
                   timeout=30, look_for_keys=False, allow_agent=False)
    print("Connected.\n", flush=True)

    os.makedirs(LOCAL_DATA_DIR, exist_ok=True)

    # 1. Check if h5py is available on the HPC
    print("=== 1. Check Python h5py on HPC ===", flush=True)
    run(client, "python3 -c 'import h5py, numpy; print(h5py.__version__, numpy.__version__)'")

    # 2. Generate ICs for each point
    print("\n=== 2. Generate ICs ===", flush=True)
    for point_name in POINTS:
        point_dir = f"~/dsidm_project/nbody_verify/{point_name}"
        # Run gen_ic.py
        out, _ = run(client, f"cd {point_dir} && python3 gen_ic.py 2>&1 | tail -15")
        # Verify ic.hdf5 was created
        run(client, f"ls -la {point_dir}/ic.hdf5 2>&1")

    # 3. Create output directories and submit Slurm jobs
    print("\n=== 3. Submit Slurm jobs ===", flush=True)
    submitted = []
    for point_name in POINTS:
        point_dir = f"~/dsidm_project/nbody_verify/{point_name}"
        # Make output dir
        run(client, f"mkdir -p {point_dir}/output")
        # Submit
        out, err = run(client, f"cd {point_dir} && sbatch submit.sh")
        # Parse job ID
        for line in out.split('\n'):
            if 'Submitted batch job' in line:
                job_id = line.split()[-1]
                submitted.append((point_name, job_id))
                print(f"  {point_name}: job {job_id}", flush=True)
                break

    if not submitted:
        print("\nERROR: No jobs submitted. Check error messages above.", flush=True)
        client.close()
        return

    print(f"\nSubmitted {len(submitted)} jobs:", flush=True)
    for name, jid in submitted:
        print(f"  {name}: {jid}", flush=True)

    # 4. Wait for jobs to finish (poll squeue)
    print("\n=== 4. Poll for job completion ===", flush=True)
    max_wait_min = 60  # 1 hour
    poll_interval = 30  # seconds
    waited = 0
    while waited < max_wait_min * 60:
        out, _ = run(client, "squeue -u scg7816 --noheader 2>&1")
        # Count lines that mention our jobs
        our_jobs_running = 0
        for name, jid in submitted:
            if jid in out:
                our_jobs_running += 1
        print(f"  [{waited//60}m {waited%60}s] {our_jobs_running}/{len(submitted)} jobs still running",
              flush=True)
        if our_jobs_running == 0:
            print("  All jobs finished!", flush=True)
            break
        time.sleep(poll_interval)
        waited += poll_interval

    # 5. Check job status for each
    print("\n=== 5. Job completion status ===", flush=True)
    for point_name, jid in submitted:
        point_dir = f"~/dsidm_project/nbody_verify/{point_name}"
        run(client, f"ls -la {point_dir}/output/ 2>&1 | head -10")
        run(client, f"cat {point_dir}/slurm_{jid}.out 2>&1 | tail -20")

    # 6. Download final snapshots
    print("\n=== 6. Download snapshots ===", flush=True)
    sftp = client.open_sftp()
    base_remote = "/public3/home/scg7816/dsidm_project/nbody_verify"
    for point_name, jid in submitted:
        point_local_dir = os.path.join(LOCAL_DATA_DIR, point_name)
        os.makedirs(point_local_dir, exist_ok=True)
        point_remote_dir = f"{base_remote}/{point_name}/output"

        # List snapshot files
        try:
            files = sftp.listdir(point_remote_dir)
        except Exception as e:
            print(f"  {point_name}: cannot list remote output dir: {e}", flush=True)
            continue

        snap_files = sorted([f for f in files if f.startswith('snap_') and f.endswith('.hdf5')])
        if not snap_files:
            print(f"  {point_name}: no snapshot files found", flush=True)
            continue

        # Download the last (final) snapshot
        last_snap = snap_files[-1]
        remote_path = f"{point_remote_dir}/{last_snap}"
        local_path = os.path.join(point_local_dir, last_snap)
        print(f"  {point_name}: downloading {last_snap}...", flush=True)
        try:
            sftp.get(remote_path, local_path)
            print(f"    saved to {local_path} ({os.path.getsize(local_path)/1e6:.1f} MB)",
                  flush=True)
        except Exception as e:
            print(f"    ERROR: {e}", flush=True)

        # Also download slurm log
        try:
            sftp.get(f"{base_remote}/{point_name}/slurm_{jid}.out",
                     os.path.join(point_local_dir, f"slurm_{jid}.out"))
        except Exception:
            pass

    sftp.close()

    client.close()
    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
