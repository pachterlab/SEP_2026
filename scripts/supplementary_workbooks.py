#!/usr/bin/env python3
"""Package saved supplementary results into uniformly formatted Excel workbooks.

Reads collated feature-weight and CHIRAL tables without rerunning analyses or
figures. Use --workbooks to select outputs; the default builds all four.
"""

import argparse
from pathlib import Path

import pandas as pd

from lib.workbooks import write_chiral_workbook, write_feature_workbook

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/supplementary_data"
CHIRAL = ROOT / "results/chiral_comparison"
WORKBOOKS = {
    1: ("Circadian", "Supplementary_Data_1_GSE54650.xlsx"),
    2: ("FUCCI", "Supplementary_Data_2_GSE121265.xlsx"),
    3: ("Endometrium", "Supplementary_Data_3_GSE234354_GSE111976.xlsx"),
    4: ("CHIRAL", "Supplementary_Data_4_CHIRAL_vs_clrcycle.xlsx"),
}
EXCEL_MAX_DATA_ROWS = 1_048_575


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workbooks", nargs="+", type=int, choices=WORKBOOKS,
                        default=list(WORKBOOKS), metavar="N", help="workbook numbers (default: 1 2 3 4)")
    parser.add_argument("--out-dir", type=Path, default=SOURCE, help="workbook output directory")
    args = parser.parse_args()
    selected = list(dict.fromkeys(args.workbooks))
    required = {}
    for number in selected:
        if number == 4:
            for name in ("chiral_vs_clrcycle.csv", "sample_phases.csv", "comparison_audit.json", "summary.csv"):
                required[CHIRAL / name] = "python scripts/supplementary_chiral_comparison.py"
        else:
            dataset, _ = WORKBOOKS[number]
            for name in ("ellipse_geometry.csv", f"{dataset.lower()}_feature_weights.csv.gz"):
                required[SOURCE / name] = "python scripts/supplementary_feature_weights.py"
    missing = [f"{path.relative_to(ROOT)} (generate with: {command})"
               for path, command in required.items() if not path.is_file()]
    if missing:
        parser.error("Missing saved results:\n" + "\n".join(missing))
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for number in selected:
        dataset, filename = WORKBOOKS[number]
        path = args.out_dir / filename
        if number == 4:
            write_chiral_workbook(path, CHIRAL)
        else:
            table = pd.read_csv(SOURCE / f"{dataset.lower()}_feature_weights.csv.gz",
                                dtype={"feature": str}, float_precision="round_trip")
            geometry = pd.read_csv(SOURCE / "ellipse_geometry.csv", float_precision="round_trip")
            geometry = geometry.loc[geometry.dataset.eq(dataset)].dropna(axis=1, how="all").copy()
            geometry["displacement_ratio"] = geometry.ellipse_center_distance / geometry.major_semiaxis
            row_counts = table.groupby("tissue").size() if "tissue" in table else [len(table)]
            if max(row_counts) > EXCEL_MAX_DATA_ROWS:
                raise ValueError(f"{dataset} exceeds the Excel row limit; use the saved CSV instead")
            write_feature_workbook(path, geometry, table, number)
        print(f"Wrote {path}", flush=True)


if __name__ == "__main__":
    main()
