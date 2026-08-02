"""Check the status of early jobs that may have finished quickly or failed."""
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
        print(out[:6000], flush=True)
    if err:
        print(f"  [stderr] {err[:3000]}", flush=True)
    return out, err


def main():
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(HOST, port=PORT, username=USER, password=PASS,
                   timeout=30, look_for_keys=False, allow_agent=False)
    print("Connected.\n", flush=True)

    # Check status of all 3 jobs
    print("=== Job status ===", flush=True)
    run(client, "sacct -j 40886164,40886167,40886168 --format=JobID,JobName,State,ExitCode,Elapsed 2>&1")

    # Check P1 output dir
    print("\n=== P1_elastic_control output ===", flush=True)
    run(client, "ls -la ~/dsidm_project/nbody_verify/P1_elastic_control/output/ 2>&1")
    run(client, "ls ~/dsidm_project/nbody_verify/P1_elastic_control/slurm_*.out 2>&1")
    # Show the slurm log
    run(client, "cat ~/dsidm_project/nbody_verify/P1_elastic_control/slurm_40886164.out 2>&1 | tail -30")

    # Check P2 output dir
    print("\n=== P2_m3_low_sigma output ===", flush=True)
    run(client, "ls -la ~/dsidm_project/nbody_verify/P2_m3_low_sigma/output/ 2>&1")
    run(client, "cat ~/dsidm_project/nbody_verify/P2_m3_low_sigma/slurm_40886167.out 2>&1 | tail -30")

    # Check P3
    print("\n=== P3_m3_high_sigma output ===", flush=True)
    run(client, "ls -la ~/dsidm_project/nbody_verify/P3_m3_high_sigma/output/ 2>&1")
    run(client, "cat ~/dsidm_project/nbody_verify/P3_m3_high_sigma/slurm_40886168.out 2>&1 | tail -30")

    # Also check err files
    print("\n=== Slurm err files ===", flush=True)
    run(client, "cat ~/dsidm_project/nbody_verify/P1_elastic_control/slurm_40886164.err 2>&1 | tail -20")
    run(client, "cat ~/dsidm_project/nbody_verify/P2_m3_low_sigma/slurm_40886167.err 2>&1 | tail -20")

    client.close()
    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
