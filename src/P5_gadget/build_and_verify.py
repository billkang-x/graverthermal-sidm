"""
Robust build script: writes build output to a remote log file, polls until done.

This replaces the unreliable finish_build.py approach. Uses nohup to detach
the build from SSH session, then polls for completion.
"""
import paramiko
import time

HOST = "ssh.cn-zhongwei-1.paracloud.com"
PORT = 22
USER = "scg7816@ZC-M6"
PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"


def run(client, cmd, timeout=120):
    """Run a quick command (with timeout). Returns stdout text."""
    print(f"$ {cmd[:200]}", flush=True)
    stdin, stdout, stderr = client.exec_command(cmd, timeout=timeout)
    try:
        out = stdout.read().decode(errors='replace').rstrip()
    except Exception as e:
        print(f"  [warn read stdout] {e}")
        out = ""
    try:
        err = stderr.read().decode(errors='replace').rstrip()
    except Exception:
        err = ""
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

    # 1. Verify run.cc edits (already applied by previous successful run)
    print("=== 1. Verify run.cc edits ===", flush=True)
    run(client, "grep -n 'sidm_dissipative' ~/dsidm_project/source/src/main/run.cc")

    # 2. Check defines_extra - use a different approach to avoid the hang.
    # Use cat and grep -F (avoid -c which sometimes hangs in this SSH session).
    print("\n=== 2. Check defines_extra ===", flush=True)
    run(client, "grep SIDM_DISSIPATIVE ~/dsidm_project/source/defines_extra || echo NOT_FOUND")
    # Add if missing (use single command, no conditional logic)
    run(client,
        "grep -q '^SIDM_DISSIPATIVE$' ~/dsidm_project/source/defines_extra || "
        "echo 'SIDM_DISSIPATIVE' >> ~/dsidm_project/source/defines_extra")
    run(client, "tail -5 ~/dsidm_project/source/defines_extra")

    # 3. Verify Makefile
    print("\n=== 3. Verify Makefile ===", flush=True)
    run(client, "grep -n 'sidm' ~/dsidm_project/source/Makefile")

    # 4. Start build in background using nohup, writing to log file
    # This detaches it from the SSH session so it survives any channel issues.
    print("\n=== 4. Launch build in background ===", flush=True)
    log_file = "~/dsidm_project/build.log"
    build_cmd = (
        "cd ~/dsidm_project/source && "
        "source ~/env.sh > /dev/null 2>&1 && "
        "make clean > /dev/null 2>&1 && "
        "make build -j 8 EXEC=Gadget4_dsidm"
    )
    # Use nohup to detach, then echo the PID
    full_cmd = f"nohup bash -c '{build_cmd}' > {log_file} 2>&1 & echo $!"
    out, _ = run(client, full_cmd)
    pid = out.strip().split('\n')[-1].strip() if out else ""
    print(f"  Build PID: {pid}", flush=True)

    # 5. Poll for completion
    print("\n=== 5. Poll for build completion ===", flush=True)
    max_wait_s = 600  # 10 minutes
    check_interval = 20  # seconds
    waited = 0
    while waited < max_wait_s:
        # Check if process is still running
        out, _ = run(client, f"if kill -0 {pid} 2>/dev/null; then echo RUNNING; else echo DONE; fi")
        status = out.strip().split('\n')[-1].strip() if out else "UNKNOWN"
        print(f"  [{waited}s] status: {status}", flush=True)
        if status == "DONE":
            break
        time.sleep(check_interval)
        waited += check_interval

    # 6. Show build log tail
    print("\n=== 6. Build log tail ===", flush=True)
    run(client, f"tail -50 {log_file}")

    # 7. Check executable
    print("\n=== 7. Check executable ===", flush=True)
    run(client, "ls -la ~/dsidm_project/source/Gadget4_dsidm 2>&1")

    # 8. Verify dissipative symbols
    print("\n=== 8. Verify dissipative symbols ===", flush=True)
    run(client, "nm ~/dsidm_project/source/Gadget4_dsidm 2>&1 | grep -i 'dissipative' | head -20")

    # 9. Count dissipative symbols
    print("\n=== 9. Count dissipative symbols ===", flush=True)
    run(client, "nm ~/dsidm_project/source/Gadget4_dsidm 2>&1 | grep -ic 'dissipative'")

    client.close()
    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
