#!/usr/bin/env python3
"""Check HPC status with the CORRECT hostname (ssh.cn-zhongwei-1.paracloud.com)."""
import paramiko, sys, os

HPC_HOST = "ssh.cn-zhongwei-1.paracloud.com"  # CORRECT hostname with ssh. prefix
HPC_PORT = 22
HPC_USER = "scg7816@ZC-M6"
HPC_PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"

REMOTE_BASE = "/public3/home/scg7816/dsidm_project"
REMOTE_SOURCE = f"{REMOTE_BASE}/source"
REMOTE_NBODY = f"{REMOTE_BASE}/nbody_verify"

POINTS = ["P1_elastic_control", "P2_m3_low_sigma", "P3_m3_high_sigma"]


def run(ssh, cmd, timeout=60):
    stdin, stdout, stderr = ssh.exec_command(cmd, timeout=timeout)
    return (stdout.read().decode('utf-8', errors='replace'),
            stderr.read().decode('utf-8', errors='replace'))


def main():
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(HPC_HOST, port=HPC_PORT, username=HPC_USER, password=HPC_PASS,
                timeout=30, allow_agent=False, look_for_keys=False)
    print("Connected to HPC!")

    # 1. Build log
    print("\n" + "=" * 60)
    print("1. /tmp/build_all.log (last 80 lines)")
    print("=" * 60)
    out, _ = run(ssh, "tail -80 /tmp/build_all.log 2>/dev/null")
    print(out if out.strip() else "(empty or missing)")

    # 2. Build process
    print("\n" + "=" * 60)
    print("2. build_all_hpc.sh process")
    print("=" * 60)
    out, _ = run(ssh, "ps aux | grep -E 'build_all_hpc|make build' | grep -v grep")
    print(out if out.strip() else "(not running)")

    # 3. Slurm queue
    print("\n" + "=" * 60)
    print("3. squeue")
    print("=" * 60)
    out, _ = run(ssh, "squeue -u scg7816 2>&1")
    print(out)

    # 4. Executables
    print("\n" + "=" * 60)
    print("4. Executables built")
    print("=" * 60)
    out, _ = run(ssh, f"ls -la {REMOTE_SOURCE}/Gadget4_P* 2>/dev/null")
    print(out if out.strip() else "(none)")

    # 5. Output dirs
    print("\n" + "=" * 60)
    print("5. Output directories and snapshots")
    print("=" * 60)
    for name in POINTS:
        d = f"{REMOTE_NBODY}/{name}/output"
        out, _ = run(ssh, f"ls -la {d}/ 2>/dev/null || echo '(no output dir)'")
        n_snaps = out.count("snapshot_")
        print(f"\n  {name}: {n_snaps} snapshot files")
        if n_snaps > 0:
            print(out.strip())

    # 6. run.log
    print("\n" + "=" * 60)
    print("6. run.log tails")
    print("=" * 60)
    for name in POINTS:
        d = f"{REMOTE_NBODY}/{name}"
        out, _ = run(ssh, f"tail -30 {d}/run.log 2>/dev/null")
        if out.strip():
            print(f"\n--- {name} run.log ---")
            print(out)

    # 7. slurm logs
    print("\n" + "=" * 60)
    print("7. Slurm logs")
    print("=" * 60)
    for name in POINTS:
        d = f"{REMOTE_NBODY}/{name}"
        out, _ = run(ssh, f"ls {d}/slurm_*.out {d}/slurm_*.err {d}/slurm-*.out {d}/slurm-*.err 2>/dev/null")
        if out.strip():
            print(f"\n--- {name} ---")
            for f in out.strip().split('\n'):
                print(f"  {f}")
                content, _ = run(ssh, f"tail -20 {f} 2>/dev/null")
                if content.strip():
                    print(content.strip())

    ssh.close()
    print("\nDone.")


if __name__ == "__main__":
    main()
