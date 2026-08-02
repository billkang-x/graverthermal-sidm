#!/usr/bin/env python3
"""Upload simulation-space IC and build/submit rescaled N-body jobs to HPC.

Strategy:
1. Upload the new sim-space IC (ic.dat) to HPC
2. Upload build_submit_rescaled.sh
3. Execute the build script remotely (rebuilds Gadget4 with new sigma values,
   creates params.txt with sim-space TimeMax, submits SLURM jobs)

The HPC already has the Gadget4 source tree at ~/dsidm_project/source/ from
prior sessions. We just rebuild with new Config.sh values.
"""
import os
import sys
import time
import paramiko

HPC_HOST = "ssh.cn-zhongwei-1.paracloud.com"
HPC_USER = "scg7816@ZC-M6"
HPC_PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"
HPC_PORT = 22

# Local files
LOCAL_IC = "D:/graverthermal-sidm/data/P5_nbody_verify/ics_sim_space/ic.dat"
LOCAL_BUILD_SCRIPT = "D:/graverthermal-sidm/src/P5_gadget/build_submit_rescaled.sh"

# Remote paths
REMOTE_BASE = "/public3/home/scg7816/dsidm_project"
REMOTE_IC_DIR = f"{REMOTE_BASE}/nbody_verify_sim/ic"
REMOTE_BUILD_SCRIPT = f"{REMOTE_BASE}/build_submit_rescaled.sh"


def connect():
    print(f"Connecting to {HPC_USER}@{HPC_HOST}...")
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(HPC_HOST, port=HPC_PORT, username=HPC_USER, password=HPC_PASS,
                timeout=30, allow_agent=False, look_for_keys=False)
    print("Connected.")
    return ssh


def run_cmd(ssh, cmd, timeout=60):
    stdin, stdout, stderr = ssh.exec_command(cmd, timeout=timeout)
    out = stdout.read().decode('utf-8', errors='replace')
    err = stderr.read().decode('utf-8', errors='replace')
    return out, err


def upload_file(sftp, local, remote, desc=""):
    size = os.path.getsize(local)
    print(f"Uploading {desc or local} ({size/1e6:.1f} MB) -> {remote}")
    sftp.put(local, remote)
    print(f"  Done.")


def main():
    # Check local files
    for f in [LOCAL_IC, LOCAL_BUILD_SCRIPT]:
        if not os.path.exists(f):
            print(f"ERROR: {f} not found")
            sys.exit(1)
    print(f"Local IC: {os.path.getsize(LOCAL_IC)/1e6:.1f} MB")

    ssh = connect()
    sftp = ssh.open_sftp()

    # Create remote directories
    print("\nCreating remote directories...")
    for d in [REMOTE_BASE, f"{REMOTE_BASE}/nbody_verify_sim",
              REMOTE_IC_DIR,
              f"{REMOTE_BASE}/nbody_verify_sim/P1_elastic_control",
              f"{REMOTE_BASE}/nbody_verify_sim/P2_m3_low_sigma",
              f"{REMOTE_BASE}/nbody_verify_sim/P3_m3_high_sigma"]:
        try:
            sftp.stat(d)
            print(f"  Exists: {d}")
        except FileNotFoundError:
            run_cmd(ssh, f"mkdir -p {d}")
            print(f"  Created: {d}")

    # Upload IC
    print("\n--- Uploading simulation-space IC ---")
    upload_file(sftp, LOCAL_IC, f"{REMOTE_IC_DIR}/ic.dat", "ic.dat (sim space, 14 MB)")

    # Upload build script
    print("\n--- Uploading build script ---")
    upload_file(sftp, LOCAL_BUILD_SCRIPT, REMOTE_BUILD_SCRIPT, "build_submit_rescaled.sh")
    run_cmd(ssh, f"chmod +x {REMOTE_BUILD_SCRIPT}")

    sftp.close()

    # Execute build script
    print("\n" + "=" * 60)
    print("Executing build script on HPC")
    print("This will: rebuild 3 executables + submit 3 SLURM jobs")
    print("Expected time: ~10-15 minutes for builds")
    print("=" * 60)

    cmd = f"cd {REMOTE_BASE} && bash build_submit_rescaled.sh 2>&1"
    stdin, stdout, stderr = ssh.exec_command(cmd, timeout=1800)  # 30 min timeout

    # Stream output line by line
    while True:
        line = stdout.readline()
        if not line:
            break
        print(line, end='', flush=True)

    err = stderr.read().decode('utf-8', errors='replace')
    if err:
        print(f"\nSTDERR:\n{err}")

    # Check job status
    print("\n" + "=" * 60)
    print("Job status")
    print("=" * 60)
    out, _ = run_cmd(ssh, "squeue -u scg7816 2>&1")
    print(out)

    # Also check if executables were built
    print("\n--- Built executables ---")
    out, _ = run_cmd(ssh, f"ls -la {REMOTE_BASE}/source/Gadget4_P* 2>&1")
    print(out)

    ssh.close()
    print("\nDone. Monitor jobs with: squeue -u scg7816")


if __name__ == '__main__':
    main()
