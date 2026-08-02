"""Compare matched dissipative and cooling-disabled direct time scans."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


KEYS = ["scan_point_id", "requested_source_time_gyr"]
COMPARISON_COLUMNS = {
    "direct_mass_chi2": "chi2",
    "direct_mass_ratio": "mass_ratio",
    "direct_vmax_kms": "vmax_kms",
    "direct_max_cooling_code": "max_cooling_code",
    "direct_steps": "steps",
}


def _completed_rows(frame: pd.DataFrame, label: str) -> pd.DataFrame:
    required = set(KEYS) | {"status"} | set(COMPARISON_COLUMNS)
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"{label} summary lacks columns: {', '.join(sorted(missing))}")
    complete = frame.loc[frame["status"] == "complete", KEYS + list(COMPARISON_COLUMNS)].copy()
    if complete.duplicated(KEYS).any():
        raise ValueError(f"{label} summary contains duplicate point/time rows")
    return complete


def compare_frames(dissipative: pd.DataFrame, elastic: pd.DataFrame) -> pd.DataFrame:
    """Return matched differences, defined as dissipative minus elastic."""
    diss = _completed_rows(dissipative, "dissipative").rename(
        columns={column: f"diss_{name}" for column, name in COMPARISON_COLUMNS.items()}
    )
    ctrl = _completed_rows(elastic, "elastic").rename(
        columns={column: f"elastic_{name}" for column, name in COMPARISON_COLUMNS.items()}
    )
    # A control run commonly contains one selected time while the dissipative
    # scan retains a larger time grid.  Require every control row to have a
    # unique dissipative counterpart, but allow additional dissipative rows.
    result = ctrl.merge(
        diss, on=KEYS, how="left", validate="one_to_one", indicator=True
    )
    if not (result["_merge"] == "both").all():
        raise ValueError("at least one elastic point/time row has no dissipative match")
    result = result.drop(columns="_merge")

    for name in COMPARISON_COLUMNS.values():
        left = pd.to_numeric(result[f"diss_{name}"], errors="raise")
        right = pd.to_numeric(result[f"elastic_{name}"], errors="raise")
        result[f"delta_{name}"] = left - right
    result["identical_mass_observables"] = (
        np.isclose(result["delta_chi2"], 0.0, rtol=0.0, atol=1e-15)
        & np.isclose(result["delta_mass_ratio"], 0.0, rtol=0.0, atol=1e-15)
    )
    return result.sort_values(KEYS).reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dissipative", required=True, help="dissipative direct-time summary CSV")
    parser.add_argument("--elastic", required=True, help="cooling-disabled direct-time summary CSV")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    result = compare_frames(pd.read_csv(args.dissipative), pd.read_csv(args.elastic))
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, index=False)
    print(f"wrote {len(result)} matched rows to {output}")
    print(result["identical_mass_observables"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
