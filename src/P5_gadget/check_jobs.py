#!/usr/bin/env python3
"""Check HPC state of N-body verification jobs - inspect output dirs, slurm logs, snapshots."""
import paramiko
import sys
import time

HOST = "ssh.cn-zhongwei-1.paracloud.com"
PORT = 22
USER = "scg7816@ZC-M6"
PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"

def run(ssh, cmd, timeout=60):
    stdin, stdout, stderr = ssh.exec_command(cmd, timeout=timeout)
    out = stdout.read().decode(errors='replace')
    err = stderr.read().decode(errors='replace')
    return out, err

def main():
    print(f"Connecting to {USER}@{HOST}:{PORT}...")
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(HOST, port=PORT, username=USER, password=PASS, timeout=30, allow_agent=False, look_for_keys=False)
    print("Connected.\n")

    # 1. Job state
    out, _ = run(ssh, "squeue -u scg7816 --noheader")
    print("=== squeue ===")
    print(out)

    # 2. sacct for the three jobs
    out, _ = run(ssh, "sacct -j 40886445,40886449,40886453 --format=JobID,JobName%20,State,Elapsed,ExitCode,MaxRSS --noheader")
    print("=== sacct ===")
    print(out)

    # 3. For each point: list output dir, slurm-*.out tail
    points = ["P1_elastic_control", "P2_m3_low_sigma", "P3_m3_high_sigma"]
    for p in points:
        base = f"~/dsidm_project/nbody_verify/{p}"
        print(f"\n=== {p} ===")
        out, _ = run(ssh, f"ls -la {base}/ 2>&1 | head -30")
        print("--- dir listing ---")
        print(out)
        out, _ = run(ssh, f"ls -la {base}/output/ 2>&1 | head -30")
        print("--- output/ ---")
        print(out)
        out, _ = run(ssh, f"ls -la {base}/restart/ 2>&1 | head -30")
        print("--- restart/ ---")
        print(out)
        out, _ = run(ssh, f"cat {base}/slurm-*.out 2>&1 | tail -60")
        print("--- slurm log tail ---")
        print(out)

    ssh.close()
    print("\nDone.")

if __name__ == "__main__":
    main()
