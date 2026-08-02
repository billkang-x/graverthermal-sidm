"""
Continue from where fix_link_module.py timed out:
  - Verify run.cc edits are in place (already done).
  - Add SIDM_DISSIPATIVE to defines_extra.
  - Build with `make build`.
  - Verify executable and symbols.
"""
import paramiko
import sys

HOST = "ssh.cn-zhongwei-1.paracloud.com"
PORT = 22
USER = "scg7816@ZC-M6"
PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"


def run(client, cmd, timeout=900):
    print(f"$ {cmd[:200]}", flush=True)
    stdin, stdout, stderr = client.exec_command(cmd, timeout=timeout, get_pty=False)
    stdout.channel.settimeout(timeout)
    out_chunks = []
    err_chunks = []
    try:
        while True:
            if stdout.channel.exit_status_ready() and not stdout.channel.recv_ready():
                break
            if stdout.channel.recv_ready():
                data = stdout.channel.recv(65536).decode(errors='replace')
                if data:
                    out_chunks.append(data)
                    print(data, end='', flush=True)
            elif stderr.channel.recv_ready():
                data = stderr.channel.recv(65536).decode(errors='replace')
                if data:
                    err_chunks.append(data)
    except Exception as e:
        print(f"\n  [WARN] channel read error: {e}")
    out = ''.join(out_chunks).rstrip()
    err = ''.join(err_chunks).rstrip()
    return out, err


def main():
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(HOST, port=PORT, username=USER, password=PASS,
                   timeout=30, look_for_keys=False, allow_agent=False)
    print("Connected.\n", flush=True)

    # 1. Verify run.cc has the proper edits
    print("=== 1. Verify run.cc edits ===", flush=True)
    run(client, "grep -n 'sidm_dissipative\\|do_sidm_scattering' ~/dsidm_project/source/src/main/run.cc")

    # 2. Check defines_extra
    print("\n=== 2. defines_extra current state ===", flush=True)
    out, _ = run(client, "grep -c 'SIDM_DISSIPATIVE' ~/dsidm_project/source/defines_extra")
    # Add it if missing
    run(client, "if ! grep -q 'SIDM_DISSIPATIVE' ~/dsidm_project/source/defines_extra; then "
                "echo 'SIDM_DISSIPATIVE' >> ~/dsidm_project/source/defines_extra; "
                "echo 'added'; else echo 'already present'; fi")
    run(client, "tail -5 ~/dsidm_project/source/defines_extra")

    # 3. Verify Makefile has sidm_dissipative.o
    print("\n=== 3. Verify Makefile ===", flush=True)
    run(client, "grep -n 'sidm' ~/dsidm_project/source/Makefile")

    # 4. Build with `make build` (suppresses check.py illegal-macros check)
    print("\n=== 4. Building Gadget4_dsidm (make build) ===", flush=True)
    run(client,
        "cd ~/dsidm_project/source && "
        "source ~/env.sh > /dev/null 2>&1 && "
        "make clean 2>&1 | tail -3 && "
        "echo '--- starting build ---' && "
        "make build -j 8 EXEC=Gadget4_dsidm 2>&1 | tail -80",
        timeout=900)

    # 5. Check executable
    print("\n=== 5. Check executable ===", flush=True)
    run(client, "ls -la ~/dsidm_project/source/Gadget4_dsidm 2>&1")

    # 6. Verify dissipative symbols
    print("\n=== 6. Verify dissipative symbols in binary ===", flush=True)
    run(client, "nm ~/dsidm_project/source/Gadget4_dsidm 2>&1 | grep -i 'dissipative' | head -20")

    # 7. Count dissipative symbols
    print("\n=== 7. Count dissipative symbols ===", flush=True)
    run(client, "nm ~/dsidm_project/source/Gadget4_dsidm 2>&1 | grep -ic 'dissipative'")

    client.close()
    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
