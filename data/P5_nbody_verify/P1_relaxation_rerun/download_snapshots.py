#!/usr/bin/env python3
"""Download Phase B snapshots from HPC and run analysis."""
import paramiko
import os

HPC_HOST = "ssh.cn-zhongwei-1.paracloud.com"
HPC_USER = "scg7816@ZC-M6"
HPC_PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"
HPC_PORT = 22
REMOTE_RUN_DIR = "/public3/home/scg7816/dsidm_project/P1_relaxation_rerun"
LOCAL_DIR = "D:/graverthermal-sidm/data/P5_nbody_verify/P1_relaxation_rerun"


def main():
    os.makedirs(f"{LOCAL_DIR}/output_prod", exist_ok=True)

    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(HPC_HOST, port=HPC_PORT, username=HPC_USER, password=HPC_PASS,
                timeout=30, allow_agent=False, look_for_keys=False)
    sftp = ssh.open_sftp()

    print("=== Downloading Phase B snapshots ===")
    for snap in ["snapshot_000", "snapshot_001", "snapshot_002"]:
        remote = f"{REMOTE_RUN_DIR}/output_prod/{snap}"
        local = f"{LOCAL_DIR}/output_prod/{snap}"
        size = os.path.getsize(local) if os.path.exists(local) else 0
        if size == 42000288:
            print(f"  {snap} already downloaded ({size/1e6:.1f} MB)")
            continue
        print(f"  Downloading {snap}...", end=" ", flush=True)
        sftp.get(remote, local)
        print(f"done ({os.path.getsize(local)/1e6:.1f} MB)")

    sftp.close()
    ssh.close()
    print("\nAll snapshots downloaded.")


if __name__ == '__main__':
    main()
