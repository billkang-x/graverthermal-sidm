#!/usr/bin/env python3
"""Upload simulation-space IC and build script to HPC, then execute.

Steps:
1. Upload ic.dat (simulation-space IC) to HPC
2. Upload build_submit_rescaled.sh
3. Execute the build script on the HPC login node
4. Monitor build + job submission

Uses paramiko for SFTP/SSH. The HPC hostname is ssh.cn-zhongwei-1.paracloud.com
with username scg7816.
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
LOCAL_ANALYZE_SCRIPT = "D:/graverthermal-sidm/src/P5_gadget/analyze_rescaled.py"
LOCAL_READ_SNAP = "D:/graverthermal-sidm/src/P5_gadget/read_binary_snap.py"

# Remote paths
REMOTE_BASE = "/public3/home/scg7816/dsidm_project"
REMOTE_IC_DIR = f"{REMOTE_BASE}/nbody_verify_sim/ic"
REMOTE_BUILD_SCRIPT = f"{REMOTE_BASE}/build_submit_rescaled.sh"


def connect():
    """Establish SSH connection."""
    print(f"Connecting to {HPC_USER}@{HPC_HOST}...")
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(HPC_HOST, port=HPC_PORT, username=HPC_USER, password=HPC_PASS,
                timeout=30, allow_agent=False, look_for_keys=False)
    print("Connected.")
    return ssh


def run_cmd(ssh, cmd, timeout=60):
    """Run a command and return stdout/stderr."""
    stdin, stdout, stderr = ssh.exec_command(cmd, timeout=timeout)
    out = stdout.read().decode('utf-8', errors='replace')
    err = stderr.read().decode('utf-8', errors='replace')
    return out, err


def upload_file(sftp, local, remote, desc=""):
    """Upload a file with progress."""
    size = os.path.getsize(local)
    print(f"Uploading {desc or local} ({size/1e6:.1f} MB) -> {remote}")
    sftp.put(local, remote)
    print(f"  Done.")


def main():
    # Check local files exist
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
              REMOTE_IC_DIR, f"{REMOTE_BASE}/source"]:
        try:
            sftp.stat(d)
        except FileNotFoundError:
            run_cmd(ssh, f"mkdir -p {d}")
            print(f"  Created {d}")

    # Upload IC
    print("\n--- Uploading IC ---")
    upload_file(sftp, LOCAL_IC, f"{REMOTE_IC_DIR}/ic.dat", "ic.dat (sim space)")

    # Upload build script
    print("\n--- Uploading build script ---")
    upload_file(sftp, LOCAL_BUILD_SCRIPT, REMOTE_BUILD_SCRIPT, "build_submit_rescaled.sh")
    run_cmd(ssh, f"chmod +x {REMOTE_BUILD_SCRIPT}")

    # Upload analysis scripts (for later use)
    print("\n--- Uploading analysis scripts ---")
    for f in [LOCAL_ANALYZE_SCRIPT, LOCAL_READ_SNAP]:
        if os.path.exists(f):
            remote_f = f"{REMOTE_BASE}/nbody_verify_sim/{os.path.basename(f)}"
            upload_file(sftp, f, remote_f, os.path.basename(f))

    sftp.close()

    # Execute build script
    print("\n--- Executing build script on HPC ---")
    print("This may take several minutes (3 builds x ~3 min each)...")
    cmd = f"cd {REMOTE_BASE} && bash build_submit_rescaled.sh 2>&1"
    stdin, stdout, stderr = ssh.exec_command(cmd, timeout=1200)  # 20 min timeout

    # Stream output
    while True:
        line = stdout.readline()
        if not line:
            break
        print(line, end='')

    err = stderr.read().decode('utf-8', errors='replace')
    if err:
        print(f"\nSTDERR:\n{err}")

    # Check job status
    print("\n--- Checking job status ---")
    out, _ = run_cmd(ssh, "squeue -u scg7816")
    print(out)

    ssh.close()
    print("\nDone. Monitor with: squeue -u scg7816")


if __name__ == '__main__':
    main()
