#!/usr/bin/env python3
"""Re-upload fixed param files and resubmit Phase A."""
import paramiko
import time
import os

HPC_HOST = "ssh.cn-zhongwei-1.paracloud.com"
HPC_USER = "scg7816@ZC-M6"
HPC_PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"
HPC_PORT = 22
REMOTE_RUN_DIR = "/public3/home/scg7816/dsidm_project/P1_relaxation_rerun"

LOCAL_PKG = "D:/graverthermal-sidm/data/P5_nbody_verify/P1_relaxation_rerun"


def run(ssh, cmd, timeout=120):
    _, stdout, stderr = ssh.exec_command(cmd, timeout=timeout)
    return stdout.read().decode('utf-8', errors='replace'), stderr.read().decode('utf-8', errors='replace')


def main():
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(HPC_HOST, port=HPC_PORT, username=HPC_USER, password=HPC_PASS,
                timeout=30, allow_agent=False, look_for_keys=False)
    sftp = ssh.open_sftp()

    print("=== Re-uploading fixed param files ===")
    for name in ["relax_params.txt", "prod_params.txt"]:
        local = f"{LOCAL_PKG}/{name}"
        remote = f"{REMOTE_RUN_DIR}/{name}"
        print(f"  {name}  ({os.path.getsize(local)} bytes)")
        sftp.put(local, remote)

    # Clean up old output
    print("\n=== Cleaning old output_relax/ ===")
    run(ssh, f"rm -rf {REMOTE_RUN_DIR}/output_relax {REMOTE_RUN_DIR}/slurm_relax_* {REMOTE_RUN_DIR}/relax_run.log 2>&1")
    out, _ = run(ssh, f"ls -la {REMOTE_RUN_DIR}/ 2>&1")
    print(out)

    sftp.close()

    # Resubmit
    print("\n=== Resubmitting Phase A ===")
    out, err = run(ssh, f"cd {REMOTE_RUN_DIR} && sbatch relax_submit.sh 2>&1", timeout=60)
    print(f"stdout: {out}")
    if err:
        print(f"stderr: {err}")

    time.sleep(5)
    print("\n=== squeue ===")
    out, _ = run(ssh, "squeue -u scg7816 2>&1")
    print(out)

    ssh.close()
    print("Done.")


if __name__ == '__main__':
    main()
