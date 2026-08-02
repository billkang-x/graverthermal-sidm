#!/usr/bin/env python3
"""Fetch the actual params.txt and usedvalues from the original P1 run on HPC."""
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

    print("=== Original P1 params.txt ===")
    out, _ = run(ssh, f"cat {P1_DIR}/params.txt")
    print(out)

    print("\n=== Original P1 params.txt-usedvalues ===")
    out, _ = run(ssh, f"cat {P1_DIR}/params.txt-usedvalues")
    print(out)

    print("\n=== Original P1 output_list.txt ===")
    out, _ = run(ssh, f"cat {P1_DIR}/output_list.txt")
    print(out)

    print("\n=== Original P1 output/ listing ===")
    out, _ = run(ssh, f"ls -la {P1_DIR}/output/")
    print(out)

    print("\n=== Original P1 run.log (head 60) ===")
    out, _ = run(ssh, f"head -60 {P1_DIR}/run.log")
    print(out)

    print("\n=== Original P1 run.log (tail 40) ===")
    out, _ = run(ssh, f"tail -40 {P1_DIR}/run.log")
    print(out)

    ssh.close()


if __name__ == '__main__':
    main()
