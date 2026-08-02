#!/usr/bin/env python3
"""Pre-flight check on HPC: verify existing source/binary/IC layout before
uploading the P1 relaxation rerun package.
"""
import paramiko

HPC_HOST = "ssh.cn-zhongwei-1.paracloud.com"
HPC_USER = "scg7816@ZC-M6"
HPC_PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"
HPC_PORT = 22
REMOTE_BASE = "/public3/home/scg7816/dsidm_project"


def run(ssh, cmd, timeout=60):
    _, stdout, stderr = ssh.exec_command(cmd, timeout=timeout)
    out = stdout.read().decode('utf-8', errors='replace')
    err = stderr.read().decode('utf-8', errors='replace')
    return out, err


def main():
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(HPC_HOST, port=HPC_PORT, username=HPC_USER, password=HPC_PASS,
                timeout=30, allow_agent=False, look_for_keys=False)
    print("Connected to HPC.\n")

    print("=== Source tree ===")
    out, _ = run(ssh, f"ls -la {REMOTE_BASE}/source/ 2>&1 | head -30")
    print(out)

    print("=== Built Gadget4 binaries ===")
    out, _ = run(ssh, f"ls -la {REMOTE_BASE}/source/Gadget4_P* 2>&1")
    print(out)

    print("=== Existing P1 run dir (for reference) ===")
    out, _ = run(ssh, f"ls -la {REMOTE_BASE}/nbody_verify_sim/P1_elastic_control/ 2>&1 | head -20")
    print(out)

    print("=== Module availability ===")
    out, _ = run(ssh, "source /public1/soft/modules/module.sh 2>/dev/null; module list 2>&1 | head -20")
    print(out)

    print("=== squeue (any running jobs?) ===")
    out, _ = run(ssh, "squeue -u scg7816 2>&1")
    print(out)

    print("=== Disk usage ===")
    out, _ = run(ssh, f"df -h {REMOTE_BASE} 2>&1; echo '---'; du -sh {REMOTE_BASE} 2>&1")
    print(out)

    ssh.close()
    print("Pre-flight done.")


if __name__ == '__main__':
    main()
