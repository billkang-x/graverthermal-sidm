"""Check full build log to diagnose compile failure."""
import paramiko

HOST = "ssh.cn-zhongwei-1.paracloud.com"
PORT = 22
USER = "scg7816@ZC-M6"
PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"

def run(client, cmd, timeout=600):
    print(f"$ {cmd[:140]}")
    stdin, stdout, stderr = client.exec_command(cmd, timeout=timeout)
    out = stdout.read().decode(errors='replace').rstrip()
    err = stderr.read().decode(errors='replace').rstrip()
    if out:
        print(out[:8000])
    if err:
        print(f"  [stderr] {err[:3000]}")
    return out, err


def main():
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(HOST, port=PORT, username=USER, password=PASS,
                   timeout=30, look_for_keys=False, allow_agent=False)
    print("Connected.\n")

    # Check the source config and try a verbose build
    run(client, "cat ~/dsidm_project/source/Config.sh")
    run(client, "cat ~/dsidm_project/source/Makefile.systype")
    run(client, "ls ~/dsidm_project/source/src/sidm/")

    # Try make and capture full output
    print("\n=== Full build attempt ===")
    run(client,
        "cd ~/dsidm_project/source && source ~/env.sh 2>&1 | tail -3; "
        "make -j 4 EXEC=Gadget4_dsidm 2>&1 | head -200",
        timeout=600)

    # Look at build.log
    run(client, "tail -50 ~/dsidm_project/source/build.log 2>&1")

    client.close()


if __name__ == '__main__':
    main()
