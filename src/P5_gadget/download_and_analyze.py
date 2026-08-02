#!/usr/bin/env python3
"""Download snapshots from HPC and run analysis.

Downloads all snapshots (snapshot_000, 001, 002) for each of P1, P2, P3
from the HPC, saves them locally, and runs the projected mass ratio analysis.

Usage:
    python download_and_analyze.py
    python download_and_analyze.py --skip-download  # use existing local snapshots
"""
import paramiko
import os
import sys
import time
import numpy as np

# Add parent dir for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

HPC_HOST = "ssh.cn-zhongwei-1.paracloud.com"
HPC_USER = "scg7816@ZC-M6"
HPC_PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"
REMOTE_BASE = "/public3/home/scg7816/dsidm_project/nbody_verify_sim"
POINTS = ["P1_elastic_control", "P2_m3_low_sigma", "P3_m3_high_sigma"]

LOCAL_BASE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "data", "P5_nbody_verify", "sim_snapshots"
)


def run_remote(ssh, cmd, timeout=30):
    stdin, stdout, stderr = ssh.exec_command(cmd, timeout=timeout)
    return (stdout.read().decode('utf-8', errors='replace'),
            stderr.read().decode('utf-8', errors='replace'))


def check_jobs_done(ssh):
    """Check if all jobs have completed (no longer in squeue)."""
    out, _ = run_remote(ssh, "squeue -u scg7816 -h | wc -l")
    return int(out.strip()) == 0


def download_snapshot(sftp, remote_path, local_path):
    """Download a single snapshot file via SFTP."""
    # Check remote file exists and get size
    try:
        st = sftp.stat(remote_path)
        remote_size = st.st_size
    except IOError:
        return False, "file not found"
    
    os.makedirs(os.path.dirname(local_path), exist_ok=True)
    
    # Check if local file already exists with same size
    if os.path.exists(local_path) and os.path.getsize(local_path) == remote_size:
        return True, f"already exists ({remote_size/1e6:.1f} MB)"
    
    # Download
    t0 = time.time()
    sftp.get(remote_path, local_path)
    dt = time.time() - t0
    local_size = os.path.getsize(local_path)
    speed = local_size / dt / 1e6 if dt > 0 else 0
    return True, f"downloaded {local_size/1e6:.1f} MB in {dt:.1f}s ({speed:.1f} MB/s)"


def main():
    skip_download = '--skip-download' in sys.argv
    
    if not skip_download:
        print("=" * 70)
        print("Connecting to HPC...")
        print("=" * 70)
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect(HPC_HOST, port=22, username=HPC_USER, password=HPC_PASS,
                    timeout=30, allow_agent=False, look_for_keys=False)
        
        # Check job status
        out, _ = run_remote(ssh, "squeue -u scg7816 2>&1")
        print("Current queue:")
        print(out)
        
        jobs_done = check_jobs_done(ssh)
        if not jobs_done:
            print("\nJobs still running. Downloading available snapshots...")
        else:
            print("\nAll jobs completed!")
        
        # List available snapshots
        sftp = ssh.open_sftp()
        
        for point in POINTS:
            remote_dir = f"{REMOTE_BASE}/{point}/output"
            out, _ = run_remote(ssh, f"ls -la {remote_dir}/snapshot_* 2>/dev/null")
            print(f"\n  {point} remote snapshots:")
            print(out)
        
        # Download all snapshots
        print("\n" + "=" * 70)
        print("Downloading snapshots...")
        print("=" * 70)
        
        downloaded = {}
        for point in POINTS:
            remote_dir = f"{REMOTE_BASE}/{point}/output"
            local_dir = os.path.join(LOCAL_BASE, point)
            downloaded[point] = []
            
            for snap_name in ['snapshot_000', 'snapshot_001', 'snapshot_002',
                              'snapshot_003', 'snapshot_004']:
                remote_path = f"{remote_dir}/{snap_name}"
                local_path = os.path.join(local_dir, snap_name)
                
                ok, msg = download_snapshot(sftp, remote_path, local_path)
                if ok:
                    print(f"  {point}/{snap_name}: {msg}")
                    downloaded[point].append(local_path)
                else:
                    pass  # snapshot doesn't exist yet
            
            if not downloaded[point]:
                print(f"  {point}: no snapshots downloaded")
        
        sftp.close()
        ssh.close()
    else:
        # Use existing local snapshots
        print("Skipping download, using local snapshots")
        for point in POINTS:
            local_dir = os.path.join(LOCAL_BASE, point)
            downloaded = {point: []}
        downloaded = {}
        for point in POINTS:
            local_dir = os.path.join(LOCAL_BASE, point)
            snaps = []
            if os.path.isdir(local_dir):
                for snap_name in sorted(os.listdir(local_dir)):
                    if snap_name.startswith('snapshot_'):
                        snaps.append(os.path.join(local_dir, snap_name))
            downloaded[point] = snaps
    
    # Run analysis
    print("\n" + "=" * 70)
    print("Running analysis...")
    print("=" * 70)
    
    from analyze_rescaled import analyze_snapshot, FLUID_PREDICTIONS, R_INNER_SIM, R_OUTER_SIM
    
    print(f"R_inner_sim = {R_INNER_SIM:.4f} kpc")
    print(f"R_outer_sim = {R_OUTER_SIM:.4f} kpc")
    
    all_results = {}
    for point in POINTS:
        snaps = downloaded.get(point, [])
        if not snaps:
            print(f"\n  {point}: no snapshots to analyze")
            continue
        
        print(f"\n{'='*70}")
        print(f"  Analyzing {point} ({len(snaps)} snapshots)")
        print(f"{'='*70}")
        
        results = []
        for snap_path in snaps:
            r = analyze_snapshot(snap_path, point)
            if r:
                results.append(r)
        all_results[point] = results
    
    # Combined comparison table
    print("\n" + "=" * 70)
    print("COMBINED COMPARISON TABLE")
    print("=" * 70)
    
    header = (f"{'Point':<25} {'snap':>12} {'t_code':>8} {'t_phys(Gyr)':>11} "
              f"{'ratio_2d':>10} {'rel_change':>11}")
    print(header)
    print("-" * len(header))
    
    for point in POINTS:
        results = all_results.get(point, [])
        if not results:
            continue
        
        # Find init (t=0)
        init_ratio = None
        for r in results:
            if r['t_code'] < 1e-8:
                init_ratio = r['ratio_2d']
                break
        
        for r in results:
            rel = r['ratio_2d'] / init_ratio if init_ratio else float('nan')
            snap_name = os.path.basename(r['file'])
            print(f"{point:<25} {snap_name:>12} {r['t_code']:>8.4f} "
                  f"{r['t_phys_gyr']:>11.4f} {r['ratio_2d']:>10.6f} "
                  f"{rel:>11.4f}")
    
    # Relative change comparison with fluid model
    print("\n" + "=" * 70)
    print("RELATIVE CHANGE vs FLUID MODEL")
    print("=" * 70)
    
    ref_rel = 0.0861 / 0.0841  # working reference (no SIDM)
    print(f"Working reference (no SIDM): rel_change = {ref_rel:.4f} ({(ref_rel-1)*100:+.2f}%)")
    
    print(f"\n{'Point':<25} {'N-body rel':>11} {'Fluid rel':>11} "
          f"{'SIDM-only':>11} {'Fluid SIDM':>12}")
    print("-" * 75)
    
    for point in POINTS:
        results = all_results.get(point, [])
        if not results:
            continue
        
        init_ratio = None
        final_ratio = None
        for r in results:
            if r['t_code'] < 1e-8:
                init_ratio = r['ratio_2d']
            else:
                final_ratio = r['ratio_2d']  # last one wins (largest t)
        
        if init_ratio is None or final_ratio is None:
            print(f"{point:<25} (insufficient data)")
            continue
        
        nb_rel = final_ratio / init_ratio
        fp = FLUID_PREDICTIONS[point]
        fluid_rel = fp['rel_change']
        
        # SIDM-only effect: subtract gravity (reference) contribution
        nb_sidm = nb_rel / ref_rel
        fluid_sidm = fluid_rel / ref_rel  # fluid model has no gravity-only run,
                                            # but this normalizes similarly
        
        print(f"{point:<25} {nb_rel:>11.4f} {fluid_rel:>11.4f} "
              f"{nb_sidm:>11.4f} {fluid_sidm:>12.4f}")
    
    print(f"\nLegend:")
    print(f"  N-body rel  = ratio_final / ratio_init (N-body)")
    print(f"  Fluid rel   = ratio_final / ratio_init (fluid model)")
    print(f"  SIDM-only   = N-body rel / ref rel (removes gravity effect)")
    print(f"  Fluid SIDM  = Fluid rel / ref rel")


if __name__ == '__main__':
    main()
