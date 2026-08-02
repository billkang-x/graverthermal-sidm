"""Check what Python environments/modules are available on the HPC."""
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

    # 1. Check available modules
    print("=== 1. Available modules ===", flush=True)
    run(client, "module avail 2>&1 | head -60")

    # 2. Check for python/anaconda modules
    print("\n=== 2. Python modules ===", flush=True)
    run(client, "module avail python 2>&1 | head -20")
    run(client, "module avail anaconda 2>&1 | head -20")
    run(client, "module avail conda 2>&1 | head -20")

    # 3. Check for any pre-installed python with numpy
    print("\n=== 3. Look for python with numpy ===", flush=True)
    run(client, "ls /public1/soft/ 2>&1 | head -20")
    run(client, "ls /public1/soft/python* 2>&1 | head -10")
    run(client, "find /public1/soft -maxdepth 3 -name 'python*' -type f 2>&1 | head -10")

    # 4. Check user's home for any conda installations
    print("\n=== 4. User's conda installations ===", flush=True)
    run(client, "ls ~/anaconda3/bin/python 2>&1 || ls ~/miniconda3/bin/python 2>&1 || echo 'no conda'")
    run(client, "ls ~/.local/bin/ 2>&1 | head -10")

    # 5. Check if there's a system python with numpy somewhere
    print("\n=== 5. Search for numpy installations ===", flush=True)
    run(client, "find /usr -name 'numpy' -type d 2>/dev/null | head -5")
    run(client, "find /public1 -name 'numpy' -type d 2>/dev/null | head -5")
    run(client, "find /opt -name 'numpy' -type d 2>/dev/null | head -5")

    client.close()
    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
