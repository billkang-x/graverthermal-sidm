"""Find available Slurm partitions on the HPC."""
import paramiko

HOST = "ssh.cn-zhongwei-1.paracloud.com"
PORT = 22
USER = "scg7816@ZC-M6"
PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"


def run(client, cmd, timeout=60):
    print(f"$ {cmd}", flush=True)
    stdin, stdout, stderr = client.exec_command(cmd, timeout=timeout)
    out = stdout.read().decode(errors='replace').rstrip()
    err = stderr.read().decode(errors='replace').rstrip()
    if out:
        print(out[:5000], flush=True)
    if err:
        print(f"  [stderr] {err[:3000]}", flush=True)
    return out, err


def main():
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(HOST, port=PORT, username=USER, password=PASS,
                   timeout=30, look_for_keys=False, allow_agent=False)
    print("Connected.\n", flush=True)

    # Check available partitions
    run(client, "sinfo 2>&1 | head -30")
    run(client, "sinfo -s 2>&1")
    # Check the user's queue info
    run(client, "squeue -u scg7816 2>&1 | head -10")
    # List partitions
    run(client, "sinfo -h -o '%P %a %l %D %t %g' 2>&1 | head -20")

    # Check past jobs to see what partition was used
    run(client, "sacct -u scg7816 --format=JobID,Partition,State,ExitCode -n 2>&1 | head -10")

    client.close()
    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
