#!/usr/bin/env python3
"""Check HPC build_all_hpc.sh status, job queue, and outputs."""
import paramiko, sys, time

HPC = dict(host="cn-zhongwei-1.paracloud.com", port=22,
           username="scg7816@ZC-M6", password="ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390")

def run(ssh, cmd, timeout=60):
    stdin, stdout, stderr = ssh.exec_command(cmd, timeout=timeout)
    return (stdout.read().decode('utf-8', errors='replace'),
            stderr.read().decode('utf-8', errors='replace'))

def main():
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(HPC["host"], port=HPC["port"], username=HPC["username"],
                password=HPC["password"], timeout=30)

    print("=" * 60)
    print("1. build_all_hpc.sh process status")
    print("=" * 60)
    out, err = run(ssh, "ps aux | grep -E 'build_all_hpc|make build' | grep -v grep")
    print(out if out.strip() else "(no build process running)")

    print("\n" + "=" * 60)
    print("2. /tmp/build_all.log (last 120 lines)")
    print("=" * 60)
    out, err = run(ssh, "tail -120 /tmp/build_all.log 2>/dev/null")
    print(out)
    if err.strip():
        print("STDERR:", err)

    print("\n" + "=" * 60)
    print("3. Slurm queue")
    print("=" * 60)
    out, err = run(ssh, "squeue -u scg7816 2>&1")
    print(out)

    print("\n" + "=" * 60)
    print("4. Executables built")
    print("=" * 60)
    out, err = run(ssh, "ls -la /public3/home/scg7816/dsidm_project/source/Gadget4_P* 2>/dev/null")
    print(out if out.strip() else "(no Gadget4_P* executables found)")

    print("\n" + "=" * 60)
    print("5. nbody_verify directories")
    print("=" * 60)
    for p in ["P1_elastic_control", "P2_m3_low_sigma", "P3_m3_high_sigma"]:
        out, _ = run(ssh, f"ls -la /public3/home/scg7816/dsidm_project/nbody_verify/{p}/ 2>&1")
        print(f"--- {p} ---")
        print(out)

    print("\n" + "=" * 60)
    print("6. Slurm output logs (if any)")
    print("=" * 60)
    for p in ["P1_elastic_control", "P2_m3_low_sigma", "P3_m3_high_sigma"]:
        d = f"/public3/home/scg7816/dsidm_project/nbody_verify/{p}"
        out, _ = run(ssh, f"ls {d}/slurm-*.out {d}/run.log 2>/dev/null")
        if out.strip():
            print(f"--- {p} logs ---")
            for f in out.split():
                print(f"  {f}")
        # tail run.log if exists
        out, _ = run(ssh, f"tail -30 {d}/run.log 2>/dev/null")
        if out.strip():
            print(f"--- {p} run.log tail ---")
            print(out)

    ssh.close()

if __name__ == "__main__":
    main()
