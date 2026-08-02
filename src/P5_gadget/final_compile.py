"""Add SIDM macros to defines_extra and use the proper compile.sh."""
import paramiko
import os

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
        print(out[:6000])
    if err:
        print(f"  [stderr] {err[:2000]}")
    return out, err


def main():
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(HOST, port=PORT, username=USER, password=PASS,
                   timeout=30, look_for_keys=False, allow_agent=False)
    print("Connected.\n")

    # Check the existing defines_extra from gas_sidm_bars (the working one)
    run(client, "cat ~/gas_sidm_bars/src/gadget4/defines_extra | head -50")

    # Copy the working defines_extra
    run(client, "cp ~/gas_sidm_bars/src/gadget4/defines_extra ~/dsidm_project/source/defines_extra && "
               "tail -20 ~/dsidm_project/source/defines_extra")

    # Now add SIDM_DISSIPATIVE and SIDM_R_DISS to defines_extra
    run(client, "echo '' >> ~/dsidm_project/source/defines_extra && "
               "echo 'SIDM_DISSIPATIVE' >> ~/dsidm_project/source/defines_extra && "
               "echo 'SIDM_R_DISS=1.0' >> ~/dsidm_project/source/defines_extra && "
               "tail -10 ~/dsidm_project/source/defines_extra")

    # Update compile.sh to point to our directory and use the right modules
    run(client, "cat ~/dsidm_project/source/compile.sh")

    # Use the proper compile.sh that loads modules
    run(client,
        "cd ~/dsidm_project/source && "
        "sed -i 's|~/gas_sidm_bars/src/gadget4|~/dsidm_project/source|g' compile.sh && "
        "sed -i 's|make -j8|make -j8 EXEC=Gadget4_dsidm|g' compile.sh && "
        "cat compile.sh")

    # Run compile.sh
    print("\n=== Run compile.sh ===")
    out, err = run(client,
        "cd ~/dsidm_project/source && bash compile.sh 2>&1 | tail -50", timeout=600)

    # Check result
    run(client, "ls -la ~/dsidm_project/source/Gadget4_dsidm 2>&1")

    client.close()
    print("\nDone.")


if __name__ == '__main__':
    main()
