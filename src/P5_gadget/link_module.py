"""
Link the sidm_dissipative module into GADGET-4's main integration loop.

Changes needed:
1. run.cc: Change `Sidm.do_sidm_scattering` to `SidmDissipative.do_sidm_scattering`
   under #ifdef SIDM_DISSIPATIVE, with fallback to Sidm for elastic-only.
2. Makefile: Add sidm_dissipative.o to OBJS.
3. Recompile.
"""
import paramiko
import sys

HOST = "ssh.cn-zhongwei-1.paracloud.com"
PORT = 22
USER = "scg7816@ZC-M6"
PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"


def run(client, cmd, timeout=600):
    print(f"$ {cmd[:150]}")
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
    print("Connected to HPC.\n")

    # Step 1: Backup run.cc and modify it
    print("=== Step 1: Modify run.cc to call SidmDissipative ===")
    run(client, "cp ~/dsidm_project/source/src/main/run.cc ~/dsidm_project/source/src/main/run.cc.bak")

    # Add the sidm_dissipative.h include
    run(client,
        "sed -i '/#include \"\\.\\.\\/sidm\\/sidm.h\"/a\\"
        "#ifdef SIDM_DISSIPATIVE\\n"
        "#include \"../sidm/sidm_dissipative.h\"\\n"
        "#endif' "
        "~/dsidm_project/source/src/main/run.cc")

    # Verify the include was added
    run(client, "grep -n 'sidm_dissipative' ~/dsidm_project/source/src/main/run.cc")

    # Replace the Sidm.do_sidm_scattering call with conditional
    # Old: #ifdef SIDM\n  Sidm.do_sidm_scattering(&Sp, All.TimeStep);\n#endif
    # New: #ifdef SIDM_DISSIPATIVE\n  SidmDissipative.do_sidm_scattering(&Sp, All.TimeStep);\n#elif defined(SIDM)\n  Sidm.do_sidm_scattering(&Sp, All.TimeStep);\n#endif
    run(client,
        "cd ~/dsidm_project/source/src/main && "
        "sed -i 's|#ifdef SIDM|\\n#ifdef SIDM_DISSIPATIVE\\n  SidmDissipative.do_sidm_scattering(\\&Sp, All.TimeStep);\\n#elif defined(SIDM)|' run.cc && "
        "grep -n 'do_sidm_scattering' run.cc")

    # Step 2: Modify Makefile to compile sidm_dissipative.o
    print("\n=== Step 2: Add sidm_dissipative.o to Makefile ===")
    run(client, "cp ~/dsidm_project/source/Makefile ~/dsidm_project/source/Makefile.bak")

    # Add sidm_dissipative.o after sidm.o
    run(client,
        "sed -i 's|OBJS    += sidm/sidm.o|OBJS    += sidm/sidm.o\\n"
        "OBJS    += sidm/sidm_dissipative.o|' "
        "~/dsidm_project/source/Makefile")

    run(client, "grep -n 'sidm' ~/dsidm_project/source/Makefile | head -10")

    # Step 3: Recompile
    print("\n=== Step 3: Recompile ===")
    run(client,
        "cd ~/dsidm_project/source && source ~/env.sh 2>&1 | tail -2 && "
        "make clean 2>&1 | tail -3 && "
        "make -j 8 EXEC=Gadget4_dsidm 2>&1 | tail -40", timeout=600)

    # Check if executable was produced
    run(client, "ls -la ~/dsidm_project/source/Gadget4_dsidm 2>&1")

    # Verify the new binary has the dissipative symbols
    run(client, "nm ~/dsidm_project/source/Gadget4_dsidm 2>&1 | grep -i dissipative | head -10")

    client.close()
    print("\nDone.")


if __name__ == "__main__":
    main()
