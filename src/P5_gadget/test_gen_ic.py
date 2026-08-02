"""Test IC generation on HPC for one point before submitting all 5."""
import paramiko

HOST = "ssh.cn-zhongwei-1.paracloud.com"
PORT = 22
USER = "scg7816@ZC-M6"
PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"


def run(client, cmd, timeout=300):
    print(f"$ {cmd[:200]}", flush=True)
    stdin, stdout, stderr = client.exec_command(cmd, timeout=timeout)
    try:
        out = stdout.read().decode(errors='replace').rstrip()
    except Exception:
        out = ""
    try:
        err = stderr.read().decode(errors='replace').rstrip()
    except Exception:
        err = ""
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

    # 1. Check Python h5py on HPC
    print("=== 1. Check Python packages on HPC ===", flush=True)
    run(client, "python3 -c 'import h5py, numpy; print(\"h5py\", h5py.__version__, \"numpy\", numpy.__version__)' 2>&1")
    # Fallback if python3 isn't the right interpreter
    run(client, "which python3 python 2>&1")
    run(client, "python -c 'import h5py, numpy; print(\"h5py\", h5py.__version__)' 2>&1")

    # 2. Test gen_ic.py for P1_elastic_control (only 1e5 particles - should be fast)
    print("\n=== 2. Test IC generation for P1_elastic_control ===", flush=True)
    point_dir = "~/dsidm_project/nbody_verify/P1_elastic_control"
    run(client, f"cd {point_dir} && python3 gen_ic.py 2>&1 | tail -20")
    run(client, f"ls -la {point_dir}/ic.hdf5 2>&1")

    # 3. If h5py is missing, install via pip --user
    out, _ = run(client, "python3 -c 'import h5py' 2>&1")
    if "ModuleNotFoundError" in out or "ImportError" in out:
        print("\n  Installing h5py for user...", flush=True)
        run(client, "pip3 install --user h5py 2>&1 | tail -5")
        # Retry
        run(client, f"cd {point_dir} && python3 gen_ic.py 2>&1 | tail -20")
        run(client, f"ls -la {point_dir}/ic.hdf5 2>&1")

    client.close()
    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
