"""
Upload the dissipative fSIDM module to HPC and set up project structure.
1. Create ~/dsidm_project/{source,runs}
2. Copy gadget4-master.bak into source (so we don't disturb the original)
3. Add sidm_dissipative.{h,cc} to source/src/sidm/
4. Create Config.sh with SIDM_DISSIPATIVE option
5. Attempt to compile
"""
import paramiko
import os

HOST = "ssh.cn-zhongwei-1.paracloud.com"
PORT = 22
USER = "scg7816@ZC-M6"
PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"

LOCAL_SRC = r"D:\graverthermal-sidm\src\P5_gadget"


def run(client, cmd, timeout=120):
    print(f"$ {cmd[:120]}")
    stdin, stdout, stderr = client.exec_command(cmd, timeout=timeout)
    out = stdout.read().decode(errors='replace').rstrip()
    err = stderr.read().decode(errors='replace').rstrip()
    if out:
        print(out[:3000])
    if err:
        print(f"  [stderr] {err[:1000]}")
    return out, err


def upload_file(client, local_path, remote_path):
    sftp = client.open_sftp()
    # Expand ~ in remote path via shell
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

    # 1. Create project structure
    run(client, "mkdir -p ~/dsidm_project/source ~/dsidm_project/runs/test_elastic")
    run(client, "mkdir -p ~/dsidm_project/runs/test_rdiss105")

    # 2. Copy gadget4-master.bak as base (only the source tree + buildsystem)
    # To save time, we just symlink the buildsystem and copy src
    run(client, "cd ~/dsidm_project/source && "
               "cp -r ~/gadget4-master.bak/Makefile . && "
               "cp -r ~/gadget4-master.bak/buildsystem . && "
               "cp -r ~/gadget4-master.bak/src . && "
               "cp ~/gadget4-master.bak/Template-Config.sh . && "
               "cp ~/gadget4-master.bak/defines_extra . && "
               "cp ~/gadget4-master.bak/Makefile.systype . 2>/dev/null; "
               "echo 'done'", timeout=180)

    # 3. Upload sidm_dissipative.h and .cc to src/sidm/
    run(client, "mkdir -p ~/dsidm_project/source/src/sidm")
    upload_file(client, os.path.join(LOCAL_SRC, "sidm_dissipative.h"),
                "~/dsidm_project/source/src/sidm/sidm_dissipative.h")
    upload_file(client, os.path.join(LOCAL_SRC, "sidm_dissipative.cc"),
                "~/dsidm_project/source/src/sidm/sidm_dissipative.cc")

    # Also copy the existing sidm.h and sidm.cc from gas_sidm_bars
    # (so the base SIDM class is available)
    run(client, "cp ~/gas_sidm_bars/src/gadget4/src/sidm/sidm.h "
                "~/dsidm_project/source/src/sidm/sidm.h && "
                "cp ~/gas_sidm_bars/src/gadget4/src/sidm/sidm.cc "
                "~/dsidm_project/source/src/sidm/sidm.cc")

    # 4. Create Config.sh with SIDM + SIDM_DISSIPATIVE options
    # Use the existing gas_sidm_bars Config.sh as base, add SIDM_DISSIPATIVE
    config_content = """# Config.sh for dissipative fSIDM test
# Based on gas_sidm_bars Config.sh

SIDM
SIDM_K_NEIGHBORS=32
SIDM_SIGMA_OVER_MASS=3.0
SIDM_PMAX=0.1
SIDM_DISSIPATIVE
SIDM_R_DISS=1.0

EVALBACKWARDACCOUNTING
TVISCSPLINE_HYDRO
NONBARYONIC_BARYON_KICK

HYDROGEN_MASS=1.0
GAMMA=1.6666666667

SFR
WSGPS
CR_DIFFUSION
CR_ZOOM

POWERSPEC_GRID_2D
SAVE_HSML_IN_SNAPSHOT
SUBFIND
FOF

PERIODIC
DOUBLEPRECISION=1
OUTPUT_COORDINATES_DOUBLE
OUTPUT_VELOCITIES_DOUBLE
OUTPUT_ACCELERATIONS_DOUBLE
OUTPUT_POTENTIAL_DOUBLE
OUTPUT_TIMESTEP_DOUBLE
"""
    # Upload Config.sh
    config_local = os.path.join(LOCAL_SRC, "Config.sh.dsidm_test")
    with open(config_local, "w", newline="\n") as f:
        f.write(config_content)
    upload_file(client, config_local, "~/dsidm_project/source/Config.sh")

    # Upload Makefile.systype
    systype_local = os.path.join(LOCAL_SRC, "hpc_source", "Makefile.systype")
    if os.path.exists(systype_local):
        upload_file(client, systype_local, "~/dsidm_project/source/Makefile.systype")

    # 5. Add SIDM_DISSIPATIVE and SIDM_R_DISS macros to defines_extra
    # (so they can be used in Config.sh)
    run(client, "echo '' >> ~/dsidm_project/source/defines_extra && "
               "echo '#ifndef SIDM_DISSIPATIVE' >> ~/dsidm_project/source/defines_extra && "
               "echo '#define SIDM_DISSIPATIVE 1' >> ~/dsidm_project/source/defines_extra && "
               "echo '#endif' >> ~/dsidm_project/source/defines_extra && "
               "echo '#ifndef SIDM_R_DISS' >> ~/dsidm_project/source/defines_extra && "
               "echo '#define SIDM_R_DISS 1.0' >> ~/dsidm_project/source/defines_extra && "
               "echo '#endif' >> ~/dsidm_project/source/defines_extra && "
               "cat ~/dsidm_project/source/defines_extra | tail -20")

    # 6. Compile attempt
    print("\n=== Attempting to compile Gadget4 with SIDM_DISSIPATIVE ===")
    out, err = run(client,
        "cd ~/dsidm_project/source && source ~/env.sh && "
        "make -j 8 EXEC=Gadget4_dsidm 2>&1 | tail -40", timeout=300)

    # 7. Check if the executable was created
    run(client, "ls -la ~/dsidm_project/source/Gadget4_dsidm 2>&1 || "
                "echo 'BUILD FAILED — see errors above'")

    client.close()
    print("\nDone.")


if __name__ == '__main__':
    main()
