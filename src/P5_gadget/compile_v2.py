"""Compile via a remote shell script (avoids quoting issues)."""
import paramiko

HOST = "ssh.cn-zhongwei-1.paracloud.com"
PORT = 22
USER = "scg7816@ZC-M6"
PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"

BUILD_SCRIPT = """#!/bin/bash
# Compile Gadget4 with SIDM_DISSIPATIVE
source /public1/soft/modules/module.sh 2>/dev/null
module purge || true
module load gcc/12.2 mpi/openmpi/4.1.5-gcc12.2 gsl/2.0 hdf5/1.8.13-gcc-zyq fftw/3.3.8-mpi

export CC=mpicc
export CXX=mpicxx
export FC=mpif90

cd ~/dsidm_project/source
echo "=== Clean ==="
make clean > /dev/null 2>&1
echo "=== Build ==="
make build -j8 EXEC=Gadget4_dsidm 2>&1
echo ""
echo "=== Result ==="
if [ -f Gadget4_dsidm ]; then
    ls -lh Gadget4_dsidm
    echo "SUCCESS"
else
    echo "FAILED"
fi
"""


def main():
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(HOST, port=PORT, username=USER, password=PASS,
                   timeout=30, look_for_keys=False, allow_agent=False)
    print("Connected.\n")

    # Upload the build script
    sftp = client.open_sftp()
    import os
    real_remote = os.path.expanduser("~/dsidm_build.sh")
    # Get the actual home path
    stdin, stdout, stderr = client.exec_command("echo $HOME/dsidm_build.sh", timeout=10)
    real_remote = stdout.read().decode(errors='replace').strip()
    import io
    sftp.putfo(io.BytesIO(BUILD_SCRIPT.encode()), real_remote)
    sftp.close()
    print(f"Uploaded build script -> {real_remote}\n")

    # Run it
    stdin, stdout, stderr = client.exec_command(f"bash {real_remote} 2>&1",
                                                 timeout=900)
    out = stdout.read().decode(errors='replace')
    # Print last 100 lines
    lines = out.split('\n')
    for line in lines[-100:]:
        print(line)

    client.close()
    print("\nDone.")


if __name__ == '__main__':
    main()
