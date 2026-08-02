#!/usr/bin/env python3
"""Add TreecoolFile tag to params.txt and resubmit. Also rewrite with explicit Unix EOL."""
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

    # Append TreecoolFile to each params.txt if not present
    print("=== Adding TreecoolFile to params.txt ===")
    for name in POINTS:
        remote_params = f"{base_remote}/{name}/params.txt"
        # Read current content
        with sftp.open(remote_params, 'r') as f:
            content = f.read().decode('utf-8')
        if 'TreecoolFile' not in content:
            # Append after InitGasTemp line
            new_content = content.rstrip() + "\nTreecoolFile             treecool.txt\n"
            with sftp.open(remote_params, 'w') as f:
                f.write(new_content)
            print(f"  {name}: TreecoolFile added")
        else:
            print(f"  {name}: TreecoolFile already present")
        # Strip CR if any
        run(ssh, f"sed -i 's/\\r$//g' {remote_params}")

        # Clean output dir
        run(ssh, f"rm -rf {base_remote}/{name}/output {base_remote}/{name}/restart {base_remote}/{name}/slurm_*.out {base_remote}/{name}/slurm_*.err {base_remote}/{name}/run.log 2>/dev/null; mkdir -p {base_remote}/{name}/output")
    print()

    sftp.close()

    # Submit
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

    # Wait 30s then check early state
    print("=== Waiting 30s ===")
    time.sleep(30)
    for name, jid in job_ids[:1]:
        remote_dir = f"{base_remote}/{name}"
        print(f"--- {name} ({jid}) early state ---")
        run(ssh, f"squeue -u scg7816 --noheader")
        run(ssh, f"ls -la {remote_dir}/output/")
        run(ssh, f"tail -30 {remote_dir}/slurm_{jid}.out 2>&1")
        run(ssh, f"tail -20 {remote_dir}/run.log 2>&1")
    print()

    # Poll
    print("=== Polling for completion ===")
    start = time.time()
    max_wait = 3600
    last_log_size = 0
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

        # Check log file growth for the first job
        if job_ids:
            name, jid = job_ids[0]
            try:
                st = sftp.stat(f"{base_remote}/{name}/run.log") if False else None
            except Exception:
                pass
            # Get log size
            out2, _ = run(ssh, f"wc -l {base_remote}/{name}/run.log 2>&1 || echo 'no log'")
            print(f"    run.log: {out2.splitlines()[-1] if out2 else 'no log'}")

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
        run(ssh, f"tail -100 {remote_dir}/slurm_{jid}.out 2>&1")
        run(ssh, f"tail -30 {remote_dir}/run.log 2>&1")

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
