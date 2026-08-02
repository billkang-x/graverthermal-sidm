#!/usr/bin/env python3
"""Check original P1 submit script and compare resource allocation."""
import paramiko

HPC_HOST = "ssh.cn-zhongwei-1.paracloud.com"
HPC_USER = "scg7816@ZC-M6"
HPC_PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"
HPC_PORT = 22
P1_DIR = "/public3/home/scg7816/dsidm_project/nbody_verify_sim/P1_elastic_control"


def run(ssh, cmd, timeout=60):
    _, stdout, stderr = ssh.exec_command(cmd, timeout=timeout)
    return stdout.read().decode('utf-8', errors='replace'), stderr.read().decode('utf-8', errors='replace')


def main():
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(HPC_HOST, port=HPC_PORT, username=HPC_USER, password=HPC_PASS,
                timeout=30, allow_agent=False, look_for_keys=False)

    print("=== Original P1 submit.sh ===")
    out, _ = run(ssh, f"cat {P1_DIR}/submit.sh")
    print(out)

    print("\n=== Original P1 submit_restart.sh ===")
    out, _ = run(ssh, f"cat {P1_DIR}/submit_restart.sh")
    print(out)

    print("\n=== Original P1 run.log: MPI tasks line ===")
    out, _ = run(ssh, f"grep 'Running on' {P1_DIR}/run.log | head -3")
    print(out)

    print("\n=== Original P1 run.log: first few sync-points ===")
    out, _ = run(ssh, f"grep 'Sync-Point' {P1_DIR}/run.log | head -10")
    print(out)

    print("\n=== Original P1 run.log: total sync-points ===")
    out, _ = run(ssh, f"grep -c 'Sync-Point' {P1_DIR}/run.log")
    print(out)

    print("\n=== Original P1 run.log: final time ===")
    out, _ = run(ssh, f"grep 'Final time' {P1_DIR}/run.log")
    print(out)

    print("\n=== Original P1 timings.txt (tail) ===")
    out, _ = run(ssh, f"tail -20 {P1_DIR}/output/timings.txt")
    print(out)

    print("\n=== Current Phase A: sync-point count ===")
    out, _ = run(ssh, "grep -c 'Sync-Point' /public3/home/scg7816/dsidm_project/P1_relaxation_rerun/relax_run.log 2>&1")
    print(out)

    print("\n=== Current Phase A: latest time ===")
    out, _ = run(ssh, "grep 'Sync-Point' /public3/home/scg7816/dsidm_project/P1_relaxation_rerun/relax_run.log 2>&1 | tail -5")
    print(out)

    print("\n=== Current Phase A: timebins.txt (tail) ===")
    out, _ = run(ssh, "tail -10 /public3/home/scg7816/dsidm_project/P1_relaxation_rerun/output_relax/timebins.txt 2>&1")
    print(out)

    ssh.close()


if __name__ == '__main__':
    main()
