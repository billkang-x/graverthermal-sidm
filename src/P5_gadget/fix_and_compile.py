"""Fix missing files and retry GADGET-4 compile."""
import paramiko

HOST = "ssh.cn-zhongwei-1.paracloud.com"
PORT = 22
USER = "scg7816@ZC-M6"
PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"


def run(client, cmd, timeout=300):
    print(f"$ {cmd[:120]}")
    stdin, stdout, stderr = client.exec_command(cmd, timeout=timeout)
    out = stdout.read().decode(errors='replace').rstrip()
    err = stderr.read().decode(errors='replace').rstrip()
    if out:
        print(out[:5000])
    if err:
        print(f"  [stderr] {err[:2000]}")
    return out, err


def main():
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(HOST, port=PORT, username=USER, password=PASS,
                   timeout=30, look_for_keys=False, allow_agent=False)
    print("Connected.\n")

    # Copy the documentation directory and any other missing deps
    run(client, "cp -r ~/gadget4-master.bak/documentation "
                 "~/dsidm_project/source/ && "
                 "cp -r ~/gadget4-master.bak/data ~/dsidm_project/source/ 2>/dev/null; "
                 "cp ~/gadget4-master.bak/Doxyfile ~/dsidm_project/source/ 2>/dev/null; "
                 "cp ~/gadget4-master.bak/DEVELOPERS ~/dsidm_project/source/ 2>/dev/null; "
                 "cp ~/gadget4-master.bak/.clang-format ~/dsidm_project/source/ 2>/dev/null; "
                 "cp ~/gadget4-master.bak/.gitignore ~/dsidm_project/source/ 2>/dev/null; "
                 "cp ~/gadget4-master.bak/make-examples.sh ~/dsidm_project/source/ 2>/dev/null; "
                 "echo 'done'", timeout=60)

    # Also check if there's a Config.sh issue - use the existing gas_sidm_bars
    # Config.sh as a base (it's known to compile) and just add SIDM_DISSIPATIVE
    run(client, "cat ~/gas_sidm_bars/src/gadget4/Config.sh")

    # Use a minimal Config.sh that's known to work + SIDM_DISSIPATIVE
    # First check what the working one looks like:
    run(client, "diff ~/gas_sidm_bars/src/gadget4/Config.sh ~/dsidm_project/source/Config.sh")

    # Try the existing gas_sidm_bars Config.sh + SIDM_DISSIPATIVE
    run(client,
        "cp ~/gas_sidm_bars/src/gadget4/Config.sh ~/dsidm_project/source/Config.sh && "
        "echo '' >> ~/dsidm_project/source/Config.sh && "
        "echo '# Dissipative fSIDM extension' >> ~/dsidm_project/source/Config.sh && "
        "echo 'SIDM_DISSIPATIVE' >> ~/dsidm_project/source/Config.sh && "
        "echo 'SIDM_R_DISS=1.0' >> ~/dsidm_project/source/Config.sh && "
        "cat ~/dsidm_project/source/Config.sh")

    # Check SYSTYPE
    run(client, "cat ~/dsidm_project/source/Makefile.systype 2>&1; "
                "cat ~/gas_sidm_bars/src/gadget4/Makefile.systype 2>&1")

    # Copy the working Makefile.systype
    run(client, "cp ~/gas_sidm_bars/src/gadget4/Makefile.systype ~/dsidm_project/source/Makefile.systype")

    # Also need to remove the SIDM_DISSIPATIVE defines from defines_extra (they're now in Config.sh)
    run(client, "cd ~/dsidm_project/source && "
                "head -n -7 defines_extra > defines_extra.tmp && "
                "mv defines_extra.tmp defines_extra && "
                "tail -10 defines_extra")

    # Retry build
    print("\n=== Retry compile ===")
    out, err = run(client,
        "cd ~/dsidm_project/source && source ~/env.sh && "
        "make clean 2>&1 | tail -5; "
        "make -j 8 EXEC=Gadget4_dsidm 2>&1 | tail -60", timeout=600)

    run(client, "ls -la ~/dsidm_project/source/Gadget4_dsidm 2>&1")

    client.close()
    print("\nDone.")


if __name__ == '__main__':
    main()
