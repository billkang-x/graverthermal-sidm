"""Use the working gas_sidm_bars Gadget4 as base (it already has SIDM patched in),
and just add the dissipative module to it."""
import paramiko
import os

HOST = "ssh.cn-zhongwei-1.paracloud.com"
PORT = 22
USER = "scg7816@ZC-M6"
PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"

LOCAL_SRC = r"D:\graverthermal-sidm\src\P5_gadget"


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


def upload_file(client, local_path, remote_path):
    sftp = client.open_sftp()
    stdin, stdout, stderr = client.exec_command(f"echo {remote_path}", timeout=10)
    real_remote = stdout.read().decode(errors='replace').strip()
    sftp.put(local_path, real_remote)
    sftp.close()
    print(f"  uploaded {os.path.basename(local_path)} -> {real_remote}")


def main():
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(HOST, port=PORT, username=USER, password=PASS,
                   timeout=30, look_for_keys=False, allow_agent=False)
    print("Connected.\n")

    # 1. Remove the failed attempt
    run(client, "rm -rf ~/dsidm_project/source")

    # 2. Copy the WORKING gas_sidm_bars/src/gadget4 as base
    run(client, "mkdir -p ~/dsidm_project && "
                 "cp -r ~/gas_sidm_bars/src/gadget4 ~/dsidm_project/source && "
                 "ls ~/dsidm_project/source/ | head -20", timeout=180)

    # 3. Check the source structure and the SIDM module
    run(client, "ls ~/dsidm_project/source/src/sidm/ 2>&1")
    run(client, "head -30 ~/dsidm_project/source/compile.sh 2>&1")
    run(client, "cat ~/dsidm_project/source/Config.sh 2>&1")

    # 4. Upload the dissipative module
    upload_file(client, os.path.join(LOCAL_SRC, "sidm_dissipative.h"),
                "~/dsidm_project/source/src/sidm/sidm_dissipative.h")
    upload_file(client, os.path.join(LOCAL_SRC, "sidm_dissipative.cc"),
                "~/dsidm_project/source/src/sidm/sidm_dissipative.cc")

    # 5. Add SIDM_DISSIPATIVE to Config.sh
    run(client, "echo '' >> ~/dsidm_project/source/Config.sh && "
               "echo '# Dissipative fSIDM extension' >> ~/dsidm_project/source/Config.sh && "
               "echo 'SIDM_DISSIPATIVE' >> ~/dsidm_project/source/Config.sh && "
               "echo 'SIDM_R_DISS=1.0' >> ~/dsidm_project/source/Config.sh && "
               "cat ~/dsidm_project/source/Config.sh")

    # 6. Add SIDM_DISSIPATIVE / SIDM_R_DISS macro definitions
    # We need to add these to the buildsystem or defines so the preprocessor
    # recognizes them. The existing SIDM_* macros must be defined somewhere.
    run(client, "grep -r 'SIDM_SIGMA_OVER_MASS' ~/dsidm_project/source/buildsystem ~/dsidm_project/source/defines_extra 2>&1 | head -10")
    run(client, "grep -r '#define SIDM' ~/dsidm_project/source/src 2>&1 | head -10")

    # 7. Add the macro definitions to defines_extra
    run(client, "echo '' >> ~/dsidm_project/source/defines_extra && "
               "echo '#ifndef SIDM_DISSIPATIVE' >> ~/dsidm_project/source/defines_extra && "
               "echo '#define SIDM_DISSIPATIVE' >> ~/dsidm_project/source/defines_extra && "
               "echo '#endif' >> ~/dsidm_project/source/defines_extra && "
               "echo '#ifndef SIDM_R_DISS' >> ~/dsidm_project/source/defines_extra && "
               "echo '#define SIDM_R_DISS 1.0' >> ~/dsidm_project/source/defines_extra && "
               "echo '#endif' >> ~/dsidm_project/source/defines_extra")

    # 8. Try to compile using the existing compile.sh
    print("\n=== Attempt compile via compile.sh ===")
    out, err = run(client,
        "cd ~/dsidm_project/source && source ~/env.sh 2>&1 | tail -2; "
        "cat compile.sh", timeout=60)

    # Look at compile.sh to understand how they build
    run(client, "cat ~/dsidm_project/source/compile.sh")

    # Now try make directly with the env loaded
    print("\n=== make attempt ===")
    out, err = run(client,
        "cd ~/dsidm_project/source && source ~/env.sh && "
        "make clean 2>&1 | tail -3 && "
        "make -j 4 EXEC=Gadget4_dsidm 2>&1 | tail -80", timeout=600)

    run(client, "ls -la ~/dsidm_project/source/Gadget4_dsidm 2>&1")

    client.close()
    print("\nDone.")


if __name__ == '__main__':
    main()
