"""Quick HPC connectivity and environment check.
Connect via SSH, run basic commands to verify the cluster environment.
"""
import paramiko
import sys

HOST = "ssh.cn-zhongwei-1.paracloud.com"
PORT = 22
USER = "scg7816@ZC-M6"
PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"

COMMANDS = [
    "hostname",
    "whoami",
    "uname -a",
    "cat /etc/os-release 2>/dev/null | head -3",
    "which gcc g++ mpicc mpirun srun sbatch 2>&1",
    "gcc --version 2>&1 | head -1",
    "mpicc --version 2>&1 | head -2 || echo 'no mpicc'",
    "module avail 2>&1 | head -30",
    "ls -la ~ 2>&1 | head -20",
    "find / -name 'gadget*' -type d 2>/dev/null | head -10",
    "find ~ -name 'Gadget*' 2>/dev/null | head -10",
    "ls /opt 2>/dev/null; ls /usr/local 2>/dev/null | head -10",
    "df -h ~ 2>&1 | head -3",
    "nvidia-smi 2>&1 | head -5 || echo 'no GPU'",
]


def main():
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    print(f"Connecting to {USER}@{HOST}:{PORT} ...", flush=True)
    try:
        client.connect(HOST, port=PORT, username=USER, password=PASS,
                       timeout=30, banner_timeout=30, auth_timeout=30,
                       look_for_keys=False, allow_agent=False)
    except Exception as e:
        print(f"Connection failed: {type(e).__name__}: {e}")
        return 1

    print("Connected.\n")
    for cmd in COMMANDS:
        print(f"$ {cmd}")
        try:
            stdin, stdout, stderr = client.exec_command(cmd, timeout=30)
            out = stdout.read().decode(errors='replace').rstrip()
            err = stderr.read().decode(errors='replace').rstrip()
            if out:
                print(out)
            if err:
                print(f"  [stderr] {err}")
        except Exception as e:
            print(f"  [error] {type(e).__name__}: {e}")
        print()

    client.close()
    print("Connection closed.")
    return 0


if __name__ == '__main__':
    sys.exit(main())
