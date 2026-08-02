"""Fetch the SIDM source code from HPC for inspection."""
import paramiko
import os

HOST = "ssh.cn-zhongwei-1.paracloud.com"
PORT = 22
USER = "scg7816@ZC-M6"
PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"

OUT_DIR = r"D:\graverthermal-sidm\src\P5_gadget\hpc_source"
os.makedirs(OUT_DIR, exist_ok=True)

# Files to fetch from the existing SIDM Gadget4 setup
FILES = [
    ("~/gas_sidm_bars/src/gadget4/src/sidm/sidm.h", "sidm.h"),
    ("~/gas_sidm_bars/src/gadget4/src/sidm/sidm.cc", "sidm.cc"),
    ("~/gas_sidm_bars/src/gadget4/src/sidm/sidm_check.c", "sidm_check.c"),
    ("~/gas_sidm_bars/src/gadget4/Config.sh.sidm", "Config.sh.sidm"),
    ("~/gas_sidm_bars/src/gadget4/Config.sh", "Config.sh"),
    ("~/gas_sidm_bars/src/gadget4/Makefile.systype", "Makefile.systype"),
    ("~/gas_sidm_bars/src/gadget4/compile.sh", "compile.sh"),
    ("~/gadget4-master.bak/Template-Config.sh", "Template-Config.sh"),
]


def main():
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(HOST, port=PORT, username=USER, password=PASS,
                   timeout=30, look_for_keys=False, allow_agent=False)
    print("Connected.\n")
    sftp = client.open_sftp()
    for remote, local_name in FILES:
        try:
            # Use shell expansion via exec_command + cat
            stdin, stdout, stderr = client.exec_command(f"cat {remote}",
                                                          timeout=15)
            content = stdout.read().decode(errors='replace')
            err = stderr.read().decode(errors='replace').strip()
            if err:
                print(f"  [{local_name}] stderr: {err[:200]}")
            local_path = os.path.join(OUT_DIR, local_name)
            with open(local_path, "w", encoding="utf-8") as f:
                f.write(content)
            nlines = content.count("\n")
            print(f"  ✓ {local_name}: {len(content)} bytes, {nlines} lines")
        except Exception as e:
            print(f"  ✗ {local_name}: {type(e).__name__}: {e}")
    sftp.close()
    client.close()
    print("\nDone.")


if __name__ == '__main__':
    main()
