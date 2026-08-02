#!/usr/bin/env python3
"""Check original IC particle counts and compare with our IC."""
import paramiko

HPC_HOST = "ssh.cn-zhongwei-1.paracloud.com"
HPC_USER = "scg7816@ZC-M6"
HPC_PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"
HPC_PORT = 22


def run(ssh, cmd, timeout=60):
    _, stdout, stderr = ssh.exec_command(cmd, timeout=timeout)
    return stdout.read().decode('utf-8', errors='replace'), stderr.read().decode('utf-8', errors='replace')


def main():
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(HPC_HOST, port=HPC_PORT, username=HPC_USER, password=HPC_PASS,
                timeout=30, allow_agent=False, look_for_keys=False)

    # Check original P1 IC
    print("=== Original P1 IC (ic.dat) header ===")
    out, _ = run(ssh, "python3 -c \"import struct; f=open('/public3/home/scg7816/dsidm_project/nbody_verify_sim/P1_elastic_control/ic.dat','rb'); sz=struct.unpack('I',f.read(4))[0]; hdr=f.read(sz); vals=struct.unpack('<6I6d2d2i6I2i4d2i6Ii',hdr[:struct.calcsize('<6I6d2d2i6I2i4d2i6Ii')]); print('npart:',list(vals[0:6])); print('mass:',list(vals[6:12])); print('time:',vals[12])\" 2>&1")
    print(out)

    # Check our IC
    print("\n=== Our IC (P1_relaxation_rerun/ic.dat) header ===")
    out, _ = run(ssh, "python3 -c \"import struct; f=open('/public3/home/scg7816/dsidm_project/P1_relaxation_rerun/ic.dat','rb'); sz=struct.unpack('I',f.read(4))[0]; hdr=f.read(sz); vals=struct.unpack('<6I6d2d2i6I2i4d2i6Ii',hdr[:struct.calcsize('<6I6d2d2i6I2i4d2i6Ii')]); print('npart:',list(vals[0:6])); print('mass:',list(vals[6:12])); print('time:',vals[12])\" 2>&1")
    print(out)

    # File sizes
    print("=== File sizes ===")
    out, _ = run(ssh, "ls -la /public3/home/scg7816/dsidm_project/nbody_verify_sim/P1_elastic_control/ic.dat /public3/home/scg7816/dsidm_project/P1_relaxation_rerun/ic.dat 2>&1")
    print(out)

    # md5sum
    print("=== md5sum ===")
    out, _ = run(ssh, "md5sum /public3/home/scg7816/dsidm_project/nbody_verify_sim/P1_elastic_control/ic.dat /public3/home/scg7816/dsidm_project/P1_relaxation_rerun/ic.dat 2>&1")
    print(out)

    ssh.close()


if __name__ == '__main__':
    main()
