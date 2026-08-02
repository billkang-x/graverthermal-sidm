"""Submit Slurm jobs for the 3 N-body verification points.

Skips rebuilding if the executable already exists.
Uses the correct partition 'amd_512'.
"""
import paramiko
import os
import time

HOST = "ssh.cn-zhongwei-1.paracloud.com"
PORT = 22
USER = "scg7816@ZC-M6"
PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"

POINTS = [
    ("P1_elastic_control", 0.1,   1.0),
    ("P2_m3_low_sigma",    0.005, 1.05),
    ("P3_m3_high_sigma",   0.220, 1.05),
]

SUBMIT_TEMPLATE = """#!/bin/bash
#SBATCH --job-name={name}
#SBATCH --partition=amd_512
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err

source ~/env.sh

cd $SLURM_SUBMIT_DIR
echo "=== Starting {name} at $(date) ==="
echo "Job ID: $SLURM_JOB_ID on $(hostname)"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo "Available memory: $(free -h | head -2)"

# Use single MPI task to avoid memory issues with small particle count
# (1e5 particles is too small to benefit from parallelism)
mpirun -np 1 ~/dsidm_project/source/Gadget4_{name} params.txt

echo "=== Finished {name} at $(date) ==="
ls -la output/
"""


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


def write_remote_file(client, remote_path, content):
    if remote_path.startswith('~'):
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
        print(f"\n{'=' * 60}\n=== {name} ===\n{'=' * 60}", flush=True)

        # 1. Check if executable exists
        exec_name = f"Gadget4_{name}"
        out, _ = run(client, f"ls -la ~/dsidm_project/source/{exec_name} 2>&1")
        if "No such file" in out:
            print(f"  Executable missing - need to rebuild", flush=True)
            # Rebuild
            config_content = (
                "SELFGRAVITY\nNTYPES=3\nEVALPOTENTIAL\nPRESSURE_ENTROPY_SPH\n"
                "COOLING\nSIDM\n"
                f"SIDM_SIGMA_OVER_MASS={sigma_m}\n"
                "SIDM_K_NEIGHBORS=32\nSIDM_PMAX=0.1\nDOUBLEPRECISION=1\n"
                "GADGET2_HEADER\n\n"
                "# Dissipative fSIDM extension\nSIDM_DISSIPATIVE\n"
                f"SIDM_R_DISS={r_diss}\n"
            )
            write_remote_file(client, "~/dsidm_project/source/Config.sh", config_content)
            run(client, "cat ~/dsidm_project/source/Config.sh")
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
            max_wait_s = 300
            waited = 0
            while waited < max_wait_s:
                out, _ = run(client, f"if kill -0 {pid} 2>/dev/null; then echo RUNNING; else echo DONE; fi")
                status = out.strip().split('\n')[-1].strip() if out else "UNKNOWN"
                print(f"  [{waited}s] {status}", flush=True)
                if status == "DONE":
                    break
                time.sleep(15)
                waited += 15
            run(client, f"tail -5 {log_file}")
            run(client, f"ls -la ~/dsidm_project/source/{exec_name} 2>&1")

        # 2. Write submit script
        submit_content = SUBMIT_TEMPLATE.format(name=name)
        point_dir = f"~/dsidm_project/nbody_verify/{name}"
        write_remote_file(client, f"{point_dir}/submit.sh", submit_content)
        run(client, f"chmod +x {point_dir}/submit.sh")

        # 3. Submit
        run(client, f"mkdir -p {point_dir}/output")
        out, _ = run(client, f"cd {point_dir} && sbatch submit.sh")
        for line in out.split('\n'):
            if 'Submitted batch job' in line:
                job_id = line.split()[-1]
                submitted.append((name, job_id))
                print(f"  {name}: job {job_id}", flush=True)
                break

    if not submitted:
        print("\nERROR: No jobs submitted.", flush=True)
        client.close()
        return

    print(f"\nSubmitted {len(submitted)} jobs:", flush=True)
    for name, jid in submitted:
        print(f"  {name}: {jid}", flush=True)

    # 4. Poll for completion
    print("\n=== Poll for job completion ===", flush=True)
    max_wait_min = 30
    poll_interval = 30
    waited = 0
    while waited < max_wait_min * 60:
        out, _ = run(client, "squeue -u scg7816 --noheader 2>&1")
        our_jobs = sum(1 for _, jid in submitted if jid in out)
        print(f"  [{waited//60}m {waited%60}s] {our_jobs}/{len(submitted)} jobs running", flush=True)
        if our_jobs == 0:
            print("  All jobs finished!", flush=True)
            break
        time.sleep(poll_interval)
        waited += poll_interval

    # 5. Check job results
    print("\n=== Check results ===", flush=True)
    for name, jid in submitted:
        point_dir = f"~/dsidm_project/nbody_verify/{name}"
        run(client, f"ls -la {point_dir}/output/ 2>&1 | head -10")
        run(client, f"tail -15 {point_dir}/slurm_{jid}.out 2>&1")

    # 6. Download snapshots
    print("\n=== Download snapshots ===", flush=True)
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
            print(f"  {name}: cannot list remote output: {e}", flush=True)
            continue

        snap_files = sorted([f for f in files if f.startswith('snap_') and f.endswith('.hdf5')])
        if not snap_files:
            print(f"  {name}: no snapshots", flush=True)
            continue

        last_snap = snap_files[-1]
        remote_path = f"{point_remote_dir}/{last_snap}"
        local_path = os.path.join(point_local_dir, last_snap)
        print(f"  {name}: downloading {last_snap}...", flush=True)
        try:
            sftp.get(remote_path, local_path)
            print(f"    saved ({os.path.getsize(local_path)/1e6:.1f} MB)", flush=True)
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
