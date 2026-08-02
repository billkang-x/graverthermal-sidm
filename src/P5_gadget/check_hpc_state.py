"""
Check the current state of run.cc and Makefile on HPC after the failed link attempt.
Restores from backup if needed and shows the actual line content.
"""
import paramiko

HOST = "ssh.cn-zhongwei-1.paracloud.com"
PORT = 22
USER = "scg7816@ZC-M6"
PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"


def run(client, cmd, timeout=120):
    print(f"$ {cmd[:200]}")
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
    print("Connected to HPC.\n")

    # 1. Check whether run.cc.bak exists (the original pristine copy)
    print("=== 1. Check backup existence ===")
    run(client, "ls -la ~/dsidm_project/source/src/main/run.cc* 2>&1")
    run(client, "ls -la ~/dsidm_project/source/Makefile* 2>&1")

    # 2. Show current state of run.cc (look for sidm/sidm_dissipative references)
    print("\n=== 2. Current run.cc lines mentioning sidm ===")
    run(client, "grep -n -i 'sidm' ~/dsidm_project/source/src/main/run.cc")

    # 3. Look at line 20-35 (where the include was added)
    print("\n=== 3. run.cc lines 20-40 ===")
    run(client, "sed -n '20,40p' ~/dsidm_project/source/src/main/run.cc")

    # 4. Look at line 250-270 (where the call is)
    print("\n=== 4. run.cc lines 250-275 ===")
    run(client, "sed -n '250,275p' ~/dsidm_project/source/src/main/run.cc")

    # 5. Makefile sidm lines
    print("\n=== 5. Makefile sidm lines ===")
    run(client, "grep -n 'sidm' ~/dsidm_project/source/Makefile")

    # 6. Check what defines are listed in Template-Config.sh and defines_extra
    print("\n=== 6. Template-Config.sh SIDM mentions ===")
    run(client, "grep -n -i 'SIDM' ~/dsidm_project/source/Template-Config.sh 2>&1 | head -20")
    run(client, "ls ~/dsidm_project/source/defines_extra 2>&1")
    run(client, "cat ~/dsidm_project/source/defines_extra 2>&1 | head -30")

    # 7. Check Config.sh
    print("\n=== 7. Config.sh content ===")
    run(client, "cat ~/dsidm_project/source/Config.sh 2>&1")

    # 8. Look at config.py to understand the illegal-macros check
    print("\n=== 8. config.py: search for 'Illegal macros' ===")
    run(client, "grep -rn 'Illegal macros' ~/dsidm_project/source/ 2>&1 | head -5")
    run(client, "grep -rn 'illegal' ~/dsidm_project/source/build_system/ 2>&1 | head -10")

    # 9. Check the build directory for any leftover .o files
    print("\n=== 9. Existing build dir sidm ===")
    run(client, "ls -la ~/dsidm_project/source/build/sidm/ 2>&1 | head -10")

    client.close()
    print("\nDone.")


if __name__ == "__main__":
    main()
