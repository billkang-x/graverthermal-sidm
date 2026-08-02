#!/usr/bin/env python3
"""
Poll HPC connectivity, then check build_all_hpc.sh status, job queue, and
download snapshots once jobs complete.

Usage:
    python poll_and_download.py            # one-shot check
    python poll_and_download.py --poll      # retry until HPC reachable
"""
import os, sys, time, socket, argparse, posixpath, datetime
import paramiko

HPC_HOST = "cn-zhongwei-1.paracloud.com"
HPC_PORT = 22
HPC_USER = "scg7816@ZC-M6"
HPC_PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"

REMOTE_BASE = "/public3/home/scg7816/dsidm_project"
REMOTE_SOURCE = f"{REMOTE_BASE}/source"
REMOTE_NBODY = f"{REMOTE_BASE}/nbody_verify"

LOCAL_BASE = "D:/graverthermal-sidm/data/P5_nbody_verify"
POINTS = [
    ("P1_elastic_control", 0.1,   1.0,  0.68),
    ("P2_m3_low_sigma",    0.005, 1.05, 0.07),
    ("P3_m3_high_sigma",   0.22,  1.05, 0.10),
]


def hpc_reachable(timeout=8):
    try:
        s = socket.create_connection((HPC_HOST, HPC_PORT), timeout=timeout)
        s.close()
        return True
    except Exception:
        return False


def connect():
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(HPC_HOST, port=HPC_PORT, username=HPC_USER, password=HPC_PASS,
                timeout=30, allow_agent=False, look_for_keys=False)
    return ssh


def run(ssh, cmd, timeout=60):
    stdin, stdout, stderr = ssh.exec_command(cmd, timeout=timeout)
    return (stdout.read().decode('utf-8', errors='replace'),
            stderr.read().decode('utf-8', errors='replace'))


def ts():
    return datetime.datetime.now().strftime("%H:%M:%S")


def check_status(ssh):
    """Full status report: build log, jobs, executables, snapshots."""
    print(f"\n[{ts()}] === HPC STATUS REPORT ===")

    # 1. Build process still running?
    out, _ = run(ssh, "ps aux | grep -E 'build_all_hpc|make build' | grep -v grep")
    print(f"\n[1] build process: {'RUNNING' if out.strip() else 'not running'}")
    if out.strip():
        print(out.strip())

    # 2. Build log tail
    out, _ = run(ssh, "tail -60 /tmp/build_all.log 2>/dev/null")
    print(f"\n[2] /tmp/build_all.log (last 60 lines):")
    print(out if out.strip() else "(empty or missing)")

    # 3. Slurm queue
    out, _ = run(ssh, "squeue -u scg7816 2>&1")
    print(f"\n[3] squeue:")
    print(out)

    # 4. Executables
    out, _ = run(ssh, f"ls -la {REMOTE_SOURCE}/Gadget4_P* 2>/dev/null")
    print(f"\n[4] Executables:")
    print(out if out.strip() else "(none built yet)")

    # 5. Output dirs + snapshots
    print(f"\n[5] Snapshot directories:")
    for name, _, _, _ in POINTS:
        d = f"{REMOTE_NBODY}/{name}/output"
        out, _ = run(ssh, f"ls -la {d}/ 2>/dev/null || echo '(no output dir)'")
        n_snaps = out.count("snapshot_")
        print(f"  {name}: {n_snaps} snapshot files")
        if n_snaps > 0:
            print(out.strip())

    # 6. slurm logs / run.log
    print(f"\n[6] Slurm/run logs:")
    for name, _, _, _ in POINTS:
        d = f"{REMOTE_NBODY}/{name}"
        out, _ = run(ssh, f"ls {d}/slurm-*.out {d}/run.log 2>/dev/null")
        if out.strip():
            print(f"  {name}:")
            print("    " + "\n    ".join(out.strip().split("\n")))
            # tail run.log
            out, _ = run(ssh, f"tail -15 {d}/run.log 2>/dev/null")
            if out.strip():
                print(f"    run.log tail:")
                for line in out.strip().split("\n"):
                    print(f"      {line}")


def download_snapshots(ssh):
    """Download snapshots from each point dir via SFTP.

    Notes on snapshot naming:
      SnapshotFileBase = "snapshot", SnapFormat = 1 (binary Gadget2-style).
      Output is typically: output/snapshot_000  (no extension), possibly
      snapdir_000/snapshot_000.0.<nfiles> for parallel IO.
      Also check for hdf5 (SnapFormat=3) just in case.
    """
    sftp = ssh.open_sftp()
    n_downloaded = 0
    snapshot_base = "snapshot"  # matches SnapshotFileBase in params.txt
    for name, _, _, _ in POINTS:
        local_dir = os.path.join(LOCAL_BASE, name)
        os.makedirs(local_dir, exist_ok=True)
        remote_dir = f"{REMOTE_NBODY}/{name}/output"
        try:
            entries = sftp.listdir(remote_dir)
        except FileNotFoundError:
            print(f"  {name}: no output dir on HPC")
            continue

        # Filter snapshot files: snapshot_NNN or snapshot_NNN.hdf5
        snaps = sorted([f for f in entries
                        if f.startswith(snapshot_base + "_") and
                        (f.endswith(".hdf5") or
                         f.split("_")[-1].isdigit())])
        # snapdirs (parallel IO)
        snapdirs = sorted([f for f in entries if f.startswith("snapdir_")])

        print(f"  {name}: {len(snaps)} snapshot files, {len(snapdirs)} snapdirs")

        for snap in snaps:
            remote_path = f"{remote_dir}/{snap}"
            local_path = os.path.join(local_dir, snap)
            print(f"    downloading {snap} ...", end=" ", flush=True)
            sftp.get(remote_path, local_path)
            print(f"done ({os.path.getsize(local_path)/1e6:.1f} MB)")
            n_downloaded += 1

        for sd in snapdirs:
            local_sd = os.path.join(local_dir, sd)
            os.makedirs(local_sd, exist_ok=True)
            try:
                sub = sftp.listdir(f"{remote_dir}/{sd}")
            except Exception:
                sub = []
            for sf in sorted(sub):
                if sf.startswith(snapshot_base) and not sf.endswith(".txt"):
                    remote_path = f"{remote_dir}/{sd}/{sf}"
                    local_path = os.path.join(local_sd, sf)
                    print(f"    downloading {sd}/{sf} ...", end=" ", flush=True)
                    sftp.get(remote_path, local_path)
                    print(f"done ({os.path.getsize(local_path)/1e6:.1f} MB)")
                    n_downloaded += 1

        # Also download run.log, slurm logs, parameters for diagnostics
        for aux in ["run.log", "parameters", "output_list.txt", "params.txt"]:
            rp = f"{REMOTE_NBODY}/{name}/{aux}"
            lp = os.path.join(local_dir, aux)
            try:
                sftp.get(rp, lp)
                print(f"    + {aux} ({os.path.getsize(lp)} bytes)")
            except FileNotFoundError:
                pass
        # slurm logs
        try:
            top = sftp.listdir(f"{REMOTE_NBODY}/{name}")
        except Exception:
            top = []
        for f in top:
            if f.startswith("slurm-") and (f.endswith(".out") or f.endswith(".err")):
                rp = f"{REMOTE_NBODY}/{name}/{f}"
                lp = os.path.join(local_dir, f)
                try:
                    sftp.get(rp, lp)
                    print(f"    + {f} ({os.path.getsize(lp)} bytes)")
                except Exception:
                    pass
    sftp.close()
    return n_downloaded


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--poll", action="store_true",
                        help="retry until HPC is reachable")
    parser.add_argument("--interval", type=int, default=60,
                        help="poll interval seconds (default 60)")
    parser.add_argument("--max-wait", type=int, default=3600,
                        help="max total wait seconds (default 3600)")
    parser.add_argument("--download", action="store_true",
                        help="download snapshots if available")
    args = parser.parse_args()

    waited = 0
    while True:
        if not hpc_reachable():
            print(f"[{ts()}] HPC unreachable (DNS/conn failed).", end="")
            if args.poll and waited < args.max_wait:
                print(f" Waiting {args.interval}s before retry...")
                time.sleep(args.interval)
                waited += args.interval
                continue
            else:
                print(" Use --poll to keep retrying.")
                return 1
        print(f"[{ts()}] HPC reachable, connecting...")
        try:
            ssh = connect()
        except Exception as e:
            print(f"  connection failed: {e}")
            if args.poll and waited < args.max_wait:
                time.sleep(args.interval)
                waited += args.interval
                continue
            return 1

        check_status(ssh)

        if args.download:
            print(f"\n[{ts()}] === DOWNLOAD SNAPSHOTS ===")
            n = download_snapshots(ssh)
            print(f"\nDownloaded {n} files total.")

        ssh.close()
        return 0


if __name__ == "__main__":
    sys.exit(main())
