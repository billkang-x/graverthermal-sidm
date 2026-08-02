#!/usr/bin/env python3
"""Check Phase A job status and tail early output."""
import paramiko
import time

HPC_HOST = "ssh.cn-zhongwei-1.paracloud.com"
HPC_USER = "scg7816@ZC-M6"
HPC_PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"
HPC_PORT = 22
REMOTE_RUN_DIR = "/public3/home/scg7816/dsidm_project/P1_relaxation_rerun"


def run(ssh, cmd, timeout=60):
    _, stdout, stderr = ssh.exec_command(cmd, timeout=timeout)
    return stdout.read().decode('utf-8', errors='replace'), stderr.read().decode('utf-8', errors='replace')


def main():
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(HPC_HOST, port=HPC_PORT, username=HPC_USER, password=HPC_PASS,
                timeout=30, allow_agent=False, look_for_keys=False)

    print("=== squeue ===")
    out, _ = run(ssh, "squeue -u scg7816 2>&1")
    print(out)

    print("\n=== slurm log files ===")
    out, _ = run(ssh, f"ls -la {REMOTE_RUN_DIR}/slurm_relax_* 2>&1")
    print(out)

    print("\n=== output_relax/ ===")
    out, _ = run(ssh, f"ls -la {REMOTE_RUN_DIR}/output_relax/ 2>&1")
    print(out)

    print("\n=== run.log (tail) ===")
    out, _ = run(ssh, f"tail -40 {REMOTE_RUN_DIR}/relax_run.log 2>&1")
    print(out)

    print("\n=== slurm .out (tail) ===")
    out, _ = run(ssh, f"tail -40 {REMOTE_RUN_DIR}/slurm_relax_*.out 2>&1")
    print(out)

    print("\n=== slurm .err (tail) ===")
    out, _ = run(ssh, f"tail -20 {REMOTE_RUN_DIR}/slurm_relax_*.err 2>&1")
    print(out)

    ssh.close()


if __name__ == '__main__':
    main()
