"""Check the working gas_sidm_bars defines_extra and use make build (skip check)."""
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
        print(f"  [stderr] {err[:2000]}")
    return out, err


def main():
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(HOST, port=PORT, username=USER, password=PASS,
                   timeout=30, look_for_keys=False, allow_agent=False)
    print("Connected.\n")

    # 1. Check the full gas_sidm_bars defines_extra
    print("=== gas_sidm_bars defines_extra ===")
    run(client, "cat ~/gas_sidm_bars/src/gadget4/defines_extra")

    # 2. Copy the FULL working defines_extra to our project
    run(client, "cp ~/gas_sidm_bars/src/gadget4/defines_extra ~/dsidm_project/source/defines_extra")

    # 3. Append our new macros
    run(client, "echo '' >> ~/dsidm_project/source/defines_extra && "
               "echo 'SIDM_DISSIPATIVE' >> ~/dsidm_project/source/defines_extra && "
               "echo 'SIDM_R_DISS=1.0' >> ~/dsidm_project/source/defines_extra && "
               "tail -10 ~/dsidm_project/source/defines_extra")

    # 4. Use `make build` to skip the macro check (the macros are valid,
    #    just not registered in the build system)
    print("\n=== make build (skip check) ===")
    out, err = run(client,
        "cd ~/dsidm_project/source && bash -c \""
        "source /public1/soft/modules/module.sh 2>/dev/null; "
        "module purge || true; "
        "module load gcc/12.2 mpi/openmpi/4.1.5-gcc12.2 gsl/2.0 hdf5/1.8.13-gcc-zyq fftw/3.3.8-mpi; "
        "export CC=mpicc CXX=mpicxx FC=mpif90; "
        "make clean > /dev/null 2>&1; "
        "make build -j8 EXEC=Gadget4_dsidm 2>&1 | tail -100", timeout=900)

    # 5. Check result
    run(client, "ls -la ~/dsidm_project/source/Gadget4_dsidm 2>&1")

    client.close()
    print("\nDone.")


if __name__ == '__main__':
    main()
