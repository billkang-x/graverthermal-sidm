#!/usr/bin/env python3
"""Upload the P1 relaxation rerun package to HPC and submit Phase A.

Steps:
  1. Create remote dir ~/dsidm_project/P1_relaxation_rerun/
  2. Upload all config/param/script files + ic.dat
  3. Verify the existing Gadget4_P1_elastic_control binary
  4. Submit Phase A (relax_submit.sh)
  5. Report job ID and monitoring instructions

Phase B is NOT submitted here — it requires the Phase-A snapshot to be
converted to an IC first (see relaxation_workflow.md step 4-5).
"""
import os
import sys
import time
import paramiko

HPC_HOST = "ssh.cn-zhongwei-1.paracloud.com"
HPC_USER = "scg7816@ZC-M6"
HPC_PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"
HPC_PORT = 22
REMOTE_BASE = "/public3/home/scg7816/dsidm_project"
REMOTE_RUN_DIR = f"{REMOTE_BASE}/P1_relaxation_rerun"

# Files to upload (local_path, remote_path, description)
LOCAL_PKG_DIR = "D:/graverthermal-sidm/data/P5_nbody_verify/P1_relaxation_rerun"
LOCAL_IC = "D:/graverthermal-sidm/data/P5_nbody_verify/ics_sim_space/ic.dat"

FILES = [
    (f"{LOCAL_PKG_DIR}/relax_params.txt",     f"{REMOTE_RUN_DIR}/relax_params.txt",     "Phase A params"),
    (f"{LOCAL_PKG_DIR}/relax_Config.sh",       f"{REMOTE_RUN_DIR}/relax_Config.sh",       "Phase A Config"),
    (f"{LOCAL_PKG_DIR}/relax_submit.sh",       f"{REMOTE_RUN_DIR}/relax_submit.sh",       "Phase A SLURM"),
    (f"{LOCAL_PKG_DIR}/prod_params.txt",       f"{REMOTE_RUN_DIR}/prod_params.txt",       "Phase B params"),
    (f"{LOCAL_PKG_DIR}/prod_Config.sh",         f"{REMOTE_RUN_DIR}/prod_Config.sh",         "Phase B Config"),
    (f"{LOCAL_PKG_DIR}/prod_submit.sh",        f"{REMOTE_RUN_DIR}/prod_submit.sh",        "Phase B SLURM"),
    (f"{LOCAL_PKG_DIR}/output_list_relax.txt", f"{REMOTE_RUN_DIR}/output_list_relax.txt", "Phase A output list"),
    (f"{LOCAL_PKG_DIR}/output_list_prod.txt",  f"{REMOTE_RUN_DIR}/output_list_prod.txt",  "Phase B output list"),
    (f"{LOCAL_PKG_DIR}/convert_restart_to_ic.py", f"{REMOTE_RUN_DIR}/convert_restart_to_ic.py", "IC converter"),
    (f"{LOCAL_PKG_DIR}/analyze_relaxation.py", f"{REMOTE_RUN_DIR}/analyze_relaxation.py", "Analyzer"),
    (f"{LOCAL_PKG_DIR}/relaxation_workflow.md", f"{REMOTE_RUN_DIR}/relaxation_workflow.md", "Workflow doc"),
    (LOCAL_IC,                                f"{REMOTE_RUN_DIR}/ic.dat",               "IC (sim-space NFW, 42 MB)"),
]


def run(ssh, cmd, timeout=120):
    _, stdout, stderr = ssh.exec_command(cmd, timeout=timeout)
    out = stdout.read().decode('utf-8', errors='replace')
    err = stderr.read().decode('utf-8', errors='replace')
    return out, err


def main():
    # Verify all local files exist
    print("=== Local file check ===")
    for local, remote, desc in FILES:
        if not os.path.exists(local):
            print(f"  MISSING: {local}")
            sys.exit(1)
        print(f"  OK  {desc:24s}  {os.path.getsize(local)/1e6:7.2f} MB")

    print("\n=== Connecting to HPC ===")
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(HPC_HOST, port=HPC_PORT, username=HPC_USER, password=HPC_PASS,
                timeout=30, allow_agent=False, look_for_keys=False)
    print("Connected.")
    sftp = ssh.open_sftp()

    # Create remote directory
    print(f"\n=== Creating {REMOTE_RUN_DIR} ===")
    run(ssh, f"mkdir -p {REMOTE_RUN_DIR}")
    out, _ = run(ssh, f"ls -la {REMOTE_RUN_DIR}")
    print(out)

    # Verify existing binary
    print("=== Verifying Gadget4_P1_elastic_control binary ===")
    out, _ = run(ssh, f"ls -la {REMOTE_BASE}/source/Gadget4_P1_elastic_control 2>&1")
    print(out)
    if "No such file" in out or not out.strip():
        print("ERROR: binary not found, aborting.")
        ssh.close()
        sys.exit(1)

    # Upload files
    print("\n=== Uploading files ===")
    for local, remote, desc in FILES:
        size = os.path.getsize(local)
        print(f"  {desc:24s}  {size/1e6:7.2f} MB  -> {remote}")
        sftp.put(local, remote)
    print("All files uploaded.")

    # Make scripts executable
    print("\n=== Setting permissions ===")
    run(ssh, f"chmod +x {REMOTE_RUN_DIR}/relax_submit.sh {REMOTE_RUN_DIR}/prod_submit.sh")
    run(ssh, f"chmod +x {REMOTE_RUN_DIR}/*.py")

    # Verify upload
    print("\n=== Remote directory listing ===")
    out, _ = run(ssh, f"ls -la {REMOTE_RUN_DIR}")
    print(out)

    # Submit Phase A
    print("\n" + "=" * 60)
    print("=== Submitting Phase A (relaxation) ===")
    print("=" * 60)
    out, err = run(ssh, f"cd {REMOTE_RUN_DIR} && sbatch relax_submit.sh 2>&1", timeout=60)
    print(f"stdout: {out}")
    if err:
        print(f"stderr: {err}")

    # Check job status
    time.sleep(3)
    print("\n=== Job queue ===")
    out, _ = run(ssh, "squeue -u scg7816 2>&1")
    print(out)

    sftp.close()
    ssh.close()

    print("\n" + "=" * 60)
    print("Phase A submitted. Next steps:")
    print("=" * 60)
    print(f"1. Monitor:  ssh {HPC_USER}@{HPC_HOST}")
    print(f"            squeue -u scg7816")
    print(f"            tail -f {REMOTE_RUN_DIR}/slurm_relax_*.out")
    print("2. When done, convert snapshot to IC:")
    print(f"   cd {REMOTE_RUN_DIR}")
    print(f"   python convert_restart_to_ic.py output_relax/snapshot_000 ic_equilibrium.dat")
    print("3. Submit Phase B:")
    print(f"   cd {REMOTE_RUN_DIR}")
    print(f"   sbatch prod_submit.sh")
    print("4. Download output_prod/snapshot_* and run analyze_relaxation.py locally")
    print()


if __name__ == '__main__':
    main()
