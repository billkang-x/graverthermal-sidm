#!/usr/bin/env python3
"""Cancel current Phase A, copy the correct IC from original P1, fix submit
script to use 16 MPI ranks, and resubmit."""
import paramiko
import time

HPC_HOST = "ssh.cn-zhongwei-1.paracloud.com"
HPC_USER = "scg7816@ZC-M6"
HPC_PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"
HPC_PORT = 22
REMOTE_BASE = "/public3/home/scg7816/dsidm_project"
REMOTE_RUN_DIR = f"{REMOTE_BASE}/P1_relaxation_rerun"
ORIG_IC = f"{REMOTE_BASE}/nbody_verify_sim/P1_elastic_control/ic.dat"


def run(ssh, cmd, timeout=120):
    _, stdout, stderr = ssh.exec_command(cmd, timeout=timeout)
    return stdout.read().decode('utf-8', errors='replace'), stderr.read().decode('utf-8', errors='replace')


def main():
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(HPC_HOST, port=HPC_PORT, username=HPC_USER, password=HPC_PASS,
                timeout=30, allow_agent=False, look_for_keys=False)

    # 1. Cancel current jobs
    print("=== Cancelling current Phase A jobs ===")
    out, _ = run(ssh, "scancel 40981665 40981695 2>&1")
    print(out if out.strip() else "  (no output)")
    time.sleep(3)

    # 2. Copy the correct IC from original P1
    print("\n=== Copying correct IC from original P1 ===")
    out, _ = run(ssh, f"cp {ORIG_IC} {REMOTE_RUN_DIR}/ic.dat 2>&1")
    print(out if out.strip() else "  Copied.")
    out, _ = run(ssh, f"ls -la {REMOTE_RUN_DIR}/ic.dat 2>&1")
    print(out)

    # 3. Verify the copied IC
    print("\n=== Verify copied IC header ===")
    out, _ = run(ssh, f"python3 -c \"import struct; f=open('{REMOTE_RUN_DIR}/ic.dat','rb'); sz=struct.unpack('I',f.read(4))[0]; hdr=f.read(sz); vals=struct.unpack('<6I6d2d2i6I2i4d2i6Ii',hdr[:struct.calcsize('<6I6d2d2i6I2i4d2i6Ii')]); print('npart:',list(vals[0:6])); print('mass:',list(vals[6:12])); print('time:',vals[12])\" 2>&1")
    print(out)

    # 4. Clean old output
    print("\n=== Cleaning old output ===")
    run(ssh, f"rm -rf {REMOTE_RUN_DIR}/output_relax {REMOTE_RUN_DIR}/slurm_relax_* {REMOTE_RUN_DIR}/relax_run.log 2>&1")

    ssh.close()
    print("\nDone. Now need to fix submit script to use 16 ranks and resubmit.")


if __name__ == '__main__':
    main()
