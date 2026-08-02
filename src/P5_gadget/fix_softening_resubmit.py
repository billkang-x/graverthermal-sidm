#!/usr/bin/env python3
"""Fix the timestep crash by increasing softening and adjusting params.

Root cause: SofteningClass1 = 0.0005 kpc (0.5 pc) is too small for this
dense NFW halo (rho_0 = 10 Msun/pc^3). The gravitational accelerations
are enormous (ac ~ 1e13), causing dt -> 0.

Fix: Increase softening to 0.01 kpc (10 pc) - still much smaller than
the 90 pc scale we care about, but large enough to avoid the timestep
crash. Also reduce MaxSizeTimestep and ErrTolForceAcc to be more
conservative.

Based on working reference:
  SofteningClass1 = 0.1 kpc (100 pc)
  ErrTolForceAcc = 0.0025
  MaxSizeTimestep = 0.002
"""
import paramiko, sys, os

HPC_HOST = "ssh.cn-zhongwei-1.paracloud.com"
HPC_PORT = 22
HPC_USER = "scg7816@ZC-M6"
HPC_PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"

REMOTE_BASE = "/public3/home/scg7816/dsidm_project/nbody_verify"
POINTS = [
    ("P1_elastic_control", 0.68),
    ("P2_m3_low_sigma", 0.07),
    ("P3_m3_high_sigma", 0.10),
]

# New softening values
SOFTENING = 0.01  # kpc = 10 pc (was 0.0005 = 0.5 pc)
MAX_TIMESTEP = 0.001  # was 0.005
ERRTOL_FORCEACC = 0.0025  # was 0.005


def run(ssh, cmd, timeout=60):
    stdin, stdout, stderr = ssh.exec_command(cmd, timeout=timeout)
    return (stdout.read().decode('utf-8', errors='replace'),
            stderr.read().decode('utf-8', errors='replace'))


def main():
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(HPC_HOST, port=HPC_PORT, username=HPC_USER, password=HPC_PASS,
                timeout=30, allow_agent=False, look_for_keys=False)
    print("Connected to HPC")

    for name, t_end in POINTS:
        print(f"\n{'='*60}")
        print(f"Fixing {name}")
        print(f"{'='*60}")

        # Read current params.txt
        out, _ = run(ssh, f"cat {REMOTE_BASE}/{name}/params.txt")
        if not out.strip():
            print(f"  ERROR: params.txt not found")
            continue

        # Replace softening values
        import re
        new_params = out
        # Replace all SofteningClass values
        new_params = re.sub(
            r'(Softening(?:Comoving|MaxPhys)Class\d+)\s+[\d.eE+-]+',
            rf'\1 {SOFTENING}',
            new_params
        )
        # Replace MaxSizeTimestep
        new_params = re.sub(
            r'MaxSizeTimestep\s+[\d.eE+-]+',
            f'MaxSizeTimestep          {MAX_TIMESTEP}',
            new_params
        )
        # Replace ErrTolForceAcc
        new_params = re.sub(
            r'ErrTolForceAcc\s+[\d.eE+-]+',
            f'ErrTolForceAcc                        {ERRTOL_FORCEACC}',
            new_params
        )

        # Write new params.txt
        sftp = ssh.open_sftp()
        with sftp.file(f"{REMOTE_BASE}/{name}/params.txt", 'w') as f:
            f.write(new_params)
        sftp.close()

        # Verify
        out, _ = run(ssh, f"grep -E 'Softening|MaxSizeTimestep|ErrTolForceAcc' {REMOTE_BASE}/{name}/params.txt")
        print(f"  New params (key lines):")
        for line in out.strip().split('\n'):
            print(f"    {line}")

        # Clean previous output
        run(ssh, f"rm -rf {REMOTE_BASE}/{name}/output {REMOTE_BASE}/{name}/restart {REMOTE_BASE}/{name}/slurm_*.out {REMOTE_BASE}/{name}/slurm_*.err {REMOTE_BASE}/{name}/run.log 2>/dev/null")
        run(ssh, f"mkdir -p {REMOTE_BASE}/{name}/output")

        # Resubmit
        out, _ = run(ssh, f"cd {REMOTE_BASE}/{name} && sbatch submit.sh")
        print(f"  Submitted: {out.strip()}")

    ssh.close()
    print("\nDone. Poll with: squeue -u scg7816")


if __name__ == "__main__":
    main()
