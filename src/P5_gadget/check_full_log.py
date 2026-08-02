"""Look at the full Gadget4 startup output to debug the memory error."""
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
        print(out[:10000], flush=True)
    if err:
        print(f"  [stderr] {err[:3000]}", flush=True)
    return out, err


def main():
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(HOST, port=PORT, username=USER, password=PASS,
                   timeout=30, look_for_keys=False, allow_agent=False)
    print("Connected.\n", flush=True)

    # Look at the FULL slurm log (not just tail)
    print("=== Full P1 slurm log ===", flush=True)
    run(client, "cat ~/dsidm_project/nbody_verify/P1_elastic_control/slurm_40886164.out 2>&1")

    # Look at the params.txt that's on the HPC
    print("\n=== Params.txt on HPC ===", flush=True)
    run(client, "cat ~/dsidm_project/nbody_verify/P1_elastic_control/params.txt 2>&1")

    # Check what the previous successful sidm-gat run looked like
    print("\n=== Previous sidm-gat job (running) ===", flush=True)
    run(client, "squeue -u scg7816 2>&1")
    # Find the working directory of the running job
    run(client, "scontrol show job 40880935 2>&1 | head -30")

    client.close()
    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
