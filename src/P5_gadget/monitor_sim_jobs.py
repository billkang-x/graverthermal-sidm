#!/usr/bin/env python3
"""Monitor the rescaled N-body jobs on HPC.

Checks:
1. Job status (squeue)
2. Run log tails
3. Snapshot files produced
4. Estimated completion time
"""
import paramiko
import sys
import time

HPC_HOST = "ssh.cn-zhongwei-1.paracloud.com"
HPC_USER = "scg7816@ZC-M6"
HPC_PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"
HPC_PORT = 22

REMOTE_BASE = "/public3/home/scg7816/dsidm_project/nbody_verify_sim"
POINTS = ["P1_elastic_control", "P2_m3_low_sigma", "P3_m3_high_sigma"]


def run(ssh, cmd, timeout=30):
    stdin, stdout, stderr = ssh.exec_command(cmd, timeout=timeout)
    return (stdout.read().decode('utf-8', errors='replace'),
            stderr.read().decode('utf-8', errors='replace'))


def main():
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(HPC_HOST, port=HPC_PORT, username=HPC_USER, password=HPC_PASS,
                timeout=30, allow_agent=False, look_for_keys=False)
    print("Connected to HPC\n")

    # 1. squeue
    print("=" * 60)
    print("1. Job queue")
    print("=" * 60)
    out, _ = run(ssh, "squeue -u scg7816 2>&1")
    print(out)

    # 2. For each point: check run.log and snapshots
    for name in POINTS:
        d = f"{REMOTE_BASE}/{name}"
        print(f"\n{'=' * 60}")
        print(f"  {name}")
        print(f"{'=' * 60}")

        # run.log tail
        out, _ = run(ssh, f"tail -20 {d}/run.log 2>/dev/null")
        if out.strip():
            print("--- run.log (last 20 lines) ---")
            print(out)
        else:
            print("(no run.log yet)")

        # snapshots
        out, _ = run(ssh, f"ls -la {d}/output/ 2>/dev/null")
        if out.strip() and "total" in out:
            print("--- snapshots ---")
            print(out)
        else:
            print("(no output yet)")

        # slurm logs
        out, _ = run(ssh, f"ls {d}/slurm_*.out {d}/slurm_*.err 2>/dev/null")
        if out.strip():
            for f in out.strip().split('\n'):
                content, _ = run(ssh, f"tail -15 {f} 2>/dev/null")
                if content.strip():
                    print(f"--- {f} ---")
                    print(content)

    # 3. Check for errors in any run.log
    print(f"\n{'=' * 60}")
    print("Error check")
    print(f"{'=' * 60}")
    for name in POINTS:
        d = f"{REMOTE_BASE}/{name}"
        out, _ = run(ssh, f"grep -i 'error\\|fatal\\|crash\\|timestep.*zero\\|segfault' {d}/run.log 2>/dev/null | tail -5")
        if out.strip():
            print(f"  {name}: ERRORS FOUND")
            print(out)
        else:
            print(f"  {name}: no errors found")

    # 4. Check restart files (for continuation)
    print(f"\n{'=' * 60}")
    print("Restart files")
    print(f"{'=' * 60}")
    for name in POINTS:
        d = f"{REMOTE_BASE}/{name}"
        out, _ = run(ssh, f"ls -la {d}/restart/ 2>/dev/null")
        if out.strip():
            print(f"  {name}: restart files exist")
            print(out)
        else:
            print(f"  {name}: no restart files")

    ssh.close()
    print("\nDone.")


if __name__ == '__main__':
    main()
