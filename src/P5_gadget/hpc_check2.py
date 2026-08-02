"""Check GADGET-4 source and SIDM-related code on HPC."""
import paramiko

HOST = "ssh.cn-zhongwei-1.paracloud.com"
PORT = 22
USER = "scg7816@ZC-M6"
PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"

COMMANDS = [
    "source ~/env.sh 2>&1; which mpicc mpirun; mpicc --version 2>&1 | head -2",
    "ls -la ~/gadget4-master.bak/ 2>&1 | head -20",
    "ls ~/gadget4-master.bak/src 2>&1 | head -40",
    "find ~/gadget4-master.bak -name 'Makefile' | head -5",
    "cat ~/gadget4-master.bak/Makefile 2>&1 | head -50 || cat ~/gadget4-master.bak/src/Makefile 2>&1 | head -50",
    "grep -l -r 'sidm\\|SIDM\\|self_interact\\|SelfInteract' ~/gadget4-master.bak/src 2>&1 | head -10",
    "grep -l -r 'sidm\\|SIDM\\|self_interact' ~/sidm-diskbar-gadget4-gate0-v2/source 2>&1 | head -10",
    "ls ~/sidm-diskbar-gadget4-gate0-v2/source 2>&1 | head -20",
    "ls ~/gas_sidm_bars/src 2>&1 | head -20",
    "ls ~/gas_sidm_bars/src/gadget4 2>&1 | head -20",
    "find ~/gas_sidm_bars -name '*.c' -o -name '*.h' 2>/dev/null | xargs grep -l 'sidm\\|dissipat\\|r_diss\\|rdiss' 2>/dev/null | head -10",
    "sinfo 2>&1 | head -10",
    "sacctmgr show user scg7816 2>&1 | head -5",
    "ls ~/bh-sidm_stage3 2>&1 | head -20",
]


def main():
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(HOST, port=PORT, username=USER, password=PASS,
                   timeout=30, look_for_keys=False, allow_agent=False)
    print("Connected.\n")
    for cmd in COMMANDS:
        # Wrap with bash -lc to source profile
        full = f"bash -lc \"{cmd}\""
        print(f"$ {cmd[:120]}")
        try:
            stdin, stdout, stderr = client.exec_command(full, timeout=30)
            out = stdout.read().decode(errors='replace').rstrip()
            err = stderr.read().decode(errors='replace').rstrip()
            if out:
                print(out[:2000])
            if err:
                print(f"  [stderr] {err[:500]}")
        except Exception as e:
            print(f"  [error] {type(e).__name__}: {e}")
        print()
    client.close()


if __name__ == '__main__':
    main()
