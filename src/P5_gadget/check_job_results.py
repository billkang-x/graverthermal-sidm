"""Quickly check job results without going through the polling script."""
import paramiko

HOST = "ssh.cn-zhongwei-1.paracloud.com"
PORT = 22
USER = "scg7816@ZC-M6"
PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"


def run(client, cmd, timeout=60):
    print(f"$ {cmd[:200]}", flush=True)
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

    # Check status
    print("=== Job status ===", flush=True)
    run(client, "sacct -j 40886445,40886449,40886453 --format=JobID,JobName,State,ExitCode,Elapsed 2>&1")

    # Look at output dirs for each point
    for name in ["P1_elastic_control", "P2_m3_low_sigma", "P3_m3_high_sigma"]:
        print(f"\n=== {name} ===", flush=True)
        run(client, f"ls -la ~/dsidm_project/nbody_verify/{name}/output/ 2>&1")
        # Get the slurm log
        run(client, f"ls ~/dsidm_project/nbody_verify/{name}/slurm_*.out 2>&1")
        # Tail the most recent one
        run(client, f"ls -t ~/dsidm_project/nbody_verify/{name}/slurm_*.out | head -1 | xargs tail -40")

    client.close()
    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
