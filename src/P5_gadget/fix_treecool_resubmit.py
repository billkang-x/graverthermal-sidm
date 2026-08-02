#!/usr/bin/env python3
"""Replace treecool.txt with the proper TREECOOL file and resubmit.

The previous treecool.txt was just an empty comment which caused Gadget4's
cooling module (COOLING is enabled in this build) to hang trying to parse it.

The source tree ships a proper TREECOOL file at:
  ~/dsidm_project/source/data/TREECOOL
Format: redshift  e1 e2 e3 e4 e5 e6 (photoheating rates, 171 lines)
"""
import paramiko
import os
import time

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

    # Cancel stale jobs
    print("=== Cancelling stale jobs ===")
    run(ssh, "squeue -u scg7816 --noheader")
    run(ssh, "scancel 40887851 40887852 40887855 2>&1 || true")
    print()

    # Replace treecool.txt with the proper TREECOOL file
    print("=== Replacing treecool.txt with proper TREECOOL file ===")
    for name in POINTS:
        remote_dir = f"{base_remote}/{name}"
        # Copy the TREECOOL file from source/data/TREECOOL
        run(ssh, f"cp ~/dsidm_project/source/data/TREECOOL {remote_dir}/treecool.txt")
        run(ssh, f"head -3 {remote_dir}/treecool.txt; wc -l {remote_dir}/treecool.txt")

        # Clean output dir
        run(ssh, f"rm -rf {remote_dir}/output {remote_dir}/restart {remote_dir}/slurm_*.out {remote_dir}/slurm_*.err {remote_dir}/run.log 2>/dev/null; mkdir -p {remote_dir}/output")
    print()

    sftp.close()

    # Submit jobs
    print("=== Submitting jobs ===")
    job_ids = []
    for name in POINTS:
        remote_dir = f"{base_remote}/{name}"
        out, _ = run(ssh, f"cd {remote_dir} && sbatch submit.sh")
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

    # Wait 90s for jobs to start and reach integration phase
    print("=== Waiting 90s for jobs to start ===")
    time.sleep(90)
    for name, jid in job_ids[:1]:
        remote_dir = f"{base_remote}/{name}"
        print(f"--- {name} ({jid}) early state ---")
        run(ssh, f"squeue -u scg7816 --noheader")
        run(ssh, f"tail -30 {remote_dir}/run.log 2>&1")
    print()

    # Poll for completion
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
        # Also show log progress
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
        run(ssh, f"tail -100 {remote_dir}/run.log 2>&1")

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
