#!/usr/bin/env python3
"""Fix DOS line breaks in submit.sh and resubmit the 3 jobs."""
import paramiko
import os
import time
import io

HOST = "ssh.cn-zhongwei-1.paracloud.com"
PORT = 22
USER = "scg7816@ZC-M6"
PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"

POINTS = [
    "P1_elastic_control",
    "P2_m3_low_sigma",
    "P3_m3_high_sigma",
]


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
    print("=== Connecting to HPC ===")
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(HOST, port=PORT, username=USER, password=PASS,
                timeout=30, look_for_keys=False, allow_agent=False)
    print("Connected.\n")
    sftp = ssh.open_sftp()
    base_remote = "/public3/home/scg7816/dsidm_project/nbody_verify"

    # First cancel any CG-state jobs left over
    print("=== Cancelling stale jobs ===")
    run(ssh, "squeue -u scg7816 --noheader")
    run(ssh, "scancel 40886445 2>&1 || true")  # CG-state P1 from before
    print()

    # Fix submit.sh on HPC: convert \r\n -> \n in place using sed
    print("=== Converting submit.sh line endings on HPC ===")
    for name in POINTS:
        remote_path = f"{base_remote}/{name}/submit.sh"
        # Use sed to strip \r
        run(ssh, f"sed -i 's/\\r$//g' {remote_path} && chmod +x {remote_path}")
        # Verify
        out, _ = run(ssh, f"file {remote_path}")
    print()

    # Verify params.txt doesn't have \r either
    print("=== Cleaning params.txt and output_list.txt line endings ===")
    for name in POINTS:
        for f in ['params.txt', 'output_list.txt']:
            remote_path = f"{base_remote}/{name}/{f}"
            run(ssh, f"sed -i 's/\\r$//g' {remote_path}")
    print()

    # Submit jobs
    print("=== Submitting jobs ===")
    job_ids = []
    for name in POINTS:
        remote_dir = f"{base_remote}/{name}"
        out, err = run(ssh, f"cd {remote_dir} && sbatch submit.sh")
        if err:
            print(f"  {name}: FAILED - {err}")
            continue
        for line in out.split('\n'):
            if 'Submitted batch job' in line:
                jid = line.split()[-1]
                job_ids.append((name, jid))
                print(f"  {name}: job {jid}")
                break
    print()

    if not job_ids:
        print("No jobs submitted. Exiting.")
        ssh.close()
        return

    # Poll for completion
    print("=== Polling for completion ===")
    start = time.time()
    max_wait = 1800
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
        print(f"  [{int(elapsed/60)}m {int(elapsed%60)}s] {our_jobs_running}/{len(job_ids)} of our jobs still in queue")
        if our_jobs_running == 0:
            print("  All jobs finished!")
            break
        time.sleep(30)

    # Check results
    print("\n=== Check results ===")
    for name, jid in job_ids:
        remote_dir = f"{base_remote}/{name}"
        print(f"\n--- {name} (job {jid}) ---")
        run(ssh, f"sacct -j {jid} --format=JobID,JobName%30,State,Elapsed,ExitCode,MaxRSS --noheader")
        run(ssh, f"ls -la {remote_dir}/output/ 2>&1 | head -30")
        run(ssh, f"tail -50 {remote_dir}/slurm_{jid}.out 2>&1")

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
