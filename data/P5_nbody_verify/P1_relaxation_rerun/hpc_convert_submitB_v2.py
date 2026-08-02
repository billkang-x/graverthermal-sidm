#!/usr/bin/env python3
"""Upload pure-python converter, run it, verify IC, submit Phase B."""
import paramiko
import time
import os

HPC_HOST = "ssh.cn-zhongwei-1.paracloud.com"
HPC_USER = "scg7816@ZC-M6"
HPC_PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"
HPC_PORT = 22
REMOTE_RUN_DIR = "/public3/home/scg7816/dsidm_project/P1_relaxation_rerun"
LOCAL_PKG = "D:/graverthermal-sidm/data/P5_nbody_verify/P1_relaxation_rerun"


def run(ssh, cmd, timeout=300):
    _, stdout, stderr = ssh.exec_command(cmd, timeout=timeout)
    return stdout.read().decode('utf-8', errors='replace'), stderr.read().decode('utf-8', errors='replace')


def main():
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(HPC_HOST, port=HPC_PORT, username=HPC_USER, password=HPC_PASS,
                timeout=30, allow_agent=False, look_for_keys=False)
    sftp = ssh.open_sftp()

    # 1. Upload pure-python converter
    print("=== Uploading pure-python converter ===")
    local = f"{LOCAL_PKG}/convert_restart_to_ic_pure.py"
    remote = f"{REMOTE_RUN_DIR}/convert_restart_to_ic_pure.py"
    print(f"  {os.path.getsize(local)} bytes")
    sftp.put(local, remote)
    run(ssh, f"chmod +x {remote}")

    # 2. Run converter
    print("\n=== Converting snapshot to IC ===")
    cmd = f"cd {REMOTE_RUN_DIR} && python3 convert_restart_to_ic_pure.py output_relax/snapshot_000 ic_equilibrium.dat 2>&1"
    out, err = run(ssh, cmd, timeout=120)
    print(out)
    if err:
        print(f"stderr: {err}")

    # 3. Verify the new IC
    print("\n=== Verify new IC ===")
    out, _ = run(ssh, f"ls -la {REMOTE_RUN_DIR}/ic_equilibrium.dat 2>&1")
    print(out)
    out, _ = run(ssh, f"python3 -c \"import struct; f=open('{REMOTE_RUN_DIR}/ic_equilibrium.dat','rb'); sz=struct.unpack('I',f.read(4))[0]; hdr=f.read(sz); vals=struct.unpack('<6I6d2d2i6I2i4d2i6Ii',hdr[:struct.calcsize('<6I6d2d2i6I2i4d2i6Ii')]); print('npart:',list(vals[0:6])); print('mass:',list(vals[6:12])); print('time:',vals[12])\" 2>&1")
    print(out)

    # 4. Submit Phase B
    print("\n=== Submitting Phase B ===")
    out, err = run(ssh, f"cd {REMOTE_RUN_DIR} && sbatch prod_submit.sh 2>&1", timeout=60)
    print(f"stdout: {out}")
    if err:
        print(f"stderr: {err}")

    time.sleep(5)
    print("\n=== squeue ===")
    out, _ = run(ssh, "squeue -u scg7816 2>&1")
    print(out)

    sftp.close()
    ssh.close()
    print("\nDone.")


if __name__ == '__main__':
    main()
