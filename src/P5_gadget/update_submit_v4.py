#!/usr/bin/env python3
"""Final attempt: use ntasks=2 + mpirun, matching the working reference."""
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


def gen_submit_sh(name):
    """Use the working reference pattern exactly: ntasks=2, mpirun -np 2."""
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
echo "Mem avail: $(free -h | head -2)"
echo "Working dir: $(pwd)"

mpirun -np 2 ~/dsidm_project/source/Gadget4_{name} params.txt 2>&1 | tee run.log
rc=$?
echo "=== Finished {name} at $(date) with exit code $rc ==="
ls -la output/ 2>&1 || true
exit 0
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
    run(ssh, "scancel 40889399 40889402 40889403 2>&1 || true")
    print()

    # Upload new submit.sh with ntasks=2 + mpirun
    print("=== Uploading new submit.sh (ntasks=2, mpirun -np 2) ===")
    for name in POINTS:
        remote_dir = f"{base_remote}/{name}"
        submit_content = gen_submit_sh(name)
        local_submit = f"D:/graverthermal-sidm/src/P5_gadget/ic_scripts/{name}_submit_v4.sh"
        with open(local_submit, 'w', newline='\n') as f:
            f.write(submit_content)
        sftp.put(local_submit, f"{remote_dir}/submit.sh")
        sftp.chmod(f"{remote_dir}/submit.sh", 0o755)
        # Strip any CR
        run(ssh, f"sed -i 's/\\r$//g' {remote_dir}/submit.sh")
        # Clean output
        run(ssh, f"rm -rf {remote_dir}/output {remote_dir}/restart {remote_dir}/slurm_*.out {remote_dir}/slurm_*.err {remote_dir}/run.log 2>/dev/null; mkdir -p {remote_dir}/output")
        print(f"  {name}: submit.sh uploaded")
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

    # Wait 2 minutes for jobs to start and reach integration
    print("=== Waiting 2 minutes ===")
    time.sleep(120)
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
        # Show log progress
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
