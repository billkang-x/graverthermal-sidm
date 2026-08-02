"""Combine direct B1938 time-scan summaries into one audit table."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def _read_summary(root: Path) -> pd.DataFrame:
    path = root / "summary.csv"
    if not path.is_file():
        raise FileNotFoundError(f"missing summary.csv under {root}")
    frame = pd.read_csv(path)
    if frame.empty:
        raise ValueError(f"empty direct summary: {path}")
    required = {"model", "source_snapshot_idx", "mass_chi2_direct", "n_shells"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"{path} lacks columns: {', '.join(sorted(missing))}")
    return frame


def combine_scan(roots: list[tuple[str, Path]]) -> pd.DataFrame:
    """Read labelled scan roots and return a stable, sorted audit table."""
    frames = []
    for label, root in roots:
        frame = _read_summary(root).copy()
        frame.insert(0, "scan_label", label)
        frame.insert(1, "scan_root", str(root.resolve()))
        frames.append(frame)
    result = pd.concat(frames, ignore_index=True, sort=False)
    return result.sort_values(
        ["model", "source_snapshot_idx", "n_shells", "scan_label"]
    ).reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root", action="append", required=True, metavar="LABEL=ROOT",
        help="direct-scan root containing summary.csv",
    )
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    roots = []
    for spec in args.root:
        if "=" not in spec:
            parser.error("--root must use LABEL=ROOT")
        label, path = spec.split("=", 1)
        root = Path(path)
        if not root.is_dir():
            parser.error(f"scan root does not exist: {root}")
        roots.append((label, root))

    frame = combine_scan(roots)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output, index=False)
    print(f"wrote {len(frame)} rows to {output}")
    print(
        frame[[
            "scan_label", "model", "source_snapshot_idx", "source_snapshot_time_gyr",
            "n_shells", "mass_chi2_direct", "mass_ratio_direct",
        ]].to_string(index=False)
    )


if __name__ == "__main__":
    main()
