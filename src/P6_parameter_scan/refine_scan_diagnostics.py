"""Regenerate continuous-refined B1938 diagnostics for a scan summary."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(_ROOT, "src", "P3_rescaling"))

from snapshot_diagnostics import (  # noqa: E402
    best_mass_fit,
    diagnose_directory,
    write_diagnostic_csv,
)


def refine_scan(scan_summary: str, filename: str, n_scan: int) -> pd.DataFrame:
    frame = pd.read_csv(scan_summary)
    if "fluid_output_dir" not in frame:
        raise ValueError("scan summary lacks fluid_output_dir")
    rows = []
    for source in frame["fluid_output_dir"].dropna().astype(str).drop_duplicates():
        diagnostic_rows = diagnose_directory(source, n_scan=n_scan)
        output = Path(source) / filename
        write_diagnostic_csv(diagnostic_rows, str(output))
        best = best_mass_fit(diagnostic_rows)
        rows.append({
            "fluid_output_dir": source,
            "diagnostic_path": str(output),
            "n_rows": len(diagnostic_rows),
            "best_snapshot_idx": None if best is None else best["snapshot_idx"],
            "best_snapshot_time_gyr": None if best is None else best["snapshot_time_gyr"],
            "best_mass_chi2": None if best is None else best["mass_chi2"],
            "best_r2D_rs": None if best is None else best["r2D_rs"],
        })
        print(output, "complete")
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scan-summary", required=True)
    parser.add_argument("--diagnostic-filename", default="b1938_snapshot_diagnostics_refined.csv")
    parser.add_argument("--n-scan", type=int, default=80)
    parser.add_argument("--output", required=True, help="summary of regenerated diagnostics")
    args = parser.parse_args()
    if args.n_scan < 8:
        parser.error("n-scan must be at least 8")
    result = refine_scan(args.scan_summary, args.diagnostic_filename, args.n_scan)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, index=False)
    print(f"wrote {len(result)} rows to {output}")


if __name__ == "__main__":
    main()
