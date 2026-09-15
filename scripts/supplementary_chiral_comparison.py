"""Reproduce the CHIRAL-vs-clrcycle apples-to-apples benchmark.

Both methods are evaluated blind (neither sees the true GSE54650 CT labels
during fitting; CT is used only for post-hoc circular alignment and error
scoring) on the exact same tissue-specific 96-gene panels that Figure 1
already selected and fit. clrcycle's numbers are read directly from the
validated ``results/fig1/tissue_summary.csv`` audit trail rather than
refitting, so this script only adds the CHIRAL side and the comparison
table. Tables are written to ``results/chiral_comparison/``.

Requires the CHIRAL R package (not on CRAN/conda-forge):

    Rscript -e 'devtools::install_github("naef-lab/CHIRAL/Pkg/CHIRAL")'
"""

from pathlib import Path
from io import StringIO
import json
import subprocess
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from lib.circadian import read_matrix, read_annotation, LABELS
from lib.common import write_versions

RAW = ROOT / "data/raw/GSE54650"
FIG1 = ROOT / "results/fig1"
OUT = ROOT / "results/chiral_comparison"
OUT.mkdir(parents=True, exist_ok=True)

TWO_PI = 2 * np.pi
PREFIX_OF = {label: prefix for prefix, label in LABELS.items()}

required = [
    RAW / "GSE54650_series_matrix.txt.gz",
    RAW / "GPL6246.annot.gz",
    FIG1 / "gene_orders.csv",
    FIG1 / "tissue_summary.csv",
]
missing = [str(p) for p in required if not p.exists()]
if missing:
    raise FileNotFoundError(
        "Missing source inputs (run scripts/download_data.py and "
        "scripts/fig1_circadian.py first):\n" + "\n".join(missing)
    )

# CHIRAL R code: read a genes-by-samples log2 matrix and CT labels, fit CHIRAL
# blind (clockgenes = every row, i.e. the full supplied panel), then align the
# inferred phases to true CT only for scoring and print one TSV result line.
CHIRAL_CODE = """
suppressMessages(library(CHIRAL))
a <- commandArgs(TRUE)
mat <- as.matrix(read.csv(a[1], row.names = 1, check.names = FALSE))
ct <- read.csv(a[2])$CT
set.seed(17)
out <- CHIRAL(mat, clockgenes = rownames(mat), pbar = FALSE)
true_phi <- 2 * pi * (ct %% 24) / 24
inf_phi <- delta.phi(true_phi, out$phi, mode = "say", median_scale = 12 / pi)
R <- abs(cor.c(true_phi, inf_phi))
adj_phi <- adjust.phases(true_phi, inf_phi)
circ_err <- pmin(abs(adj_phi - true_phi), 2 * pi - abs(adj_phi - true_phi))
# CHIRAL prints its own progress/diagnostic lines to stdout regardless of
# pbar, so the result row is marked and extracted by the caller rather than
# assumed to be the only output.
cat("===CHIRAL_RESULT===\\n")
cat("R\\tmean_abs_err_h\\titer\\n")
cat(R, "\\t", mean(circ_err) * 24 / (2 * pi), "\\t", out$iter, "\\n", sep = "")
"""
RESULT_MARKER = "===CHIRAL_RESULT==="


def build_tissue_panel(expr, annot, label, genes):
    """Return the genes-by-samples log2 matrix and CT labels for one tissue."""
    prefix = PREFIX_OF[label]
    cols = sorted(
        [c for c in expr if c.startswith(prefix + "_CT")],
        key=lambda c: int(c.split("_CT")[1]),
    )
    ct = np.array([int(c.split("_CT")[1]) for c in cols], float)
    tissue = expr[cols]
    tissue = tissue.loc[np.isfinite(tissue.to_numpy()).all(1) & (tissue.to_numpy().min(1) > 0)]

    symbol = annot.set_index("ID_REF").gene_symbol.reindex(tissue.index)
    keep = symbol.isin(genes)
    tissue, symbol = tissue.loc[keep], symbol.loc[keep]

    # One probe per gene symbol: the highest-variance probe, matching the
    # deduplication Figure 1's own selection already applied upstream.
    logvar = np.log2(tissue.to_numpy()).var(1, ddof=1)
    order = pd.DataFrame({"symbol": symbol.to_numpy(), "logvar": logvar}, index=tissue.index)
    keep_ids = order.sort_values("logvar", ascending=False).drop_duplicates("symbol").index
    tissue = np.log2(tissue.loc[keep_ids])
    tissue.index = symbol.loc[keep_ids].to_numpy()
    return tissue, ct


def run_chiral(matrix: pd.DataFrame, ct: np.ndarray, tissue_dir: Path) -> dict:
    """Write one tissue's panel to disk and fit CHIRAL on it via Rscript."""
    tissue_dir.mkdir(parents=True, exist_ok=True)
    mat_path = tissue_dir / "log2_matrix.csv"
    ct_path = tissue_dir / "ct.csv"
    matrix.to_csv(mat_path)
    pd.DataFrame({"CT": ct}).to_csv(ct_path, index=False)
    done = subprocess.run(
        ["Rscript", "-e", CHIRAL_CODE, str(mat_path), str(ct_path)],
        capture_output=True,
        text=True,
        check=True,
    )
    if RESULT_MARKER not in done.stdout:
        raise RuntimeError(f"CHIRAL produced no result marker for {tissue_dir}:\n{done.stdout}\n{done.stderr}")
    tsv = done.stdout.split(RESULT_MARKER, 1)[1].strip()
    row = pd.read_csv(StringIO(tsv), sep="\t").iloc[0]
    return {"R": float(row.R), "mean_abs_err_h": float(row.mean_abs_err_h), "iterations": int(row.iter)}


def main():
    expr = read_matrix(required[0])
    annot = read_annotation(required[1])
    gene_orders = pd.read_csv(FIG1 / "gene_orders.csv")
    clrcycle_summary = pd.read_csv(FIG1 / "tissue_summary.csv")
    clrcycle_summary["mean_abs_err_h"] = clrcycle_summary.circular_mae_rad * 24 / TWO_PI

    rows = []
    for label in LABELS.values():
        genes = gene_orders.loc[gene_orders.tissue == label, "feature"].tolist()
        matrix, ct = build_tissue_panel(expr, annot, label, genes)
        chiral = run_chiral(matrix, ct, OUT / "matrices" / label.replace(" ", "_"))
        clr = clrcycle_summary.set_index("tissue").loc[label]
        rows.append(
            {
                "tissue": label,
                "n_genes": matrix.shape[0],
                "chiral_R": chiral["R"],
                "chiral_mean_abs_err_h": chiral["mean_abs_err_h"],
                "chiral_iterations": chiral["iterations"],
                "clrcycle_R": float(clr.phase_resultant),
                "clrcycle_mean_abs_err_h": float(clr.mean_abs_err_h),
            }
        )
        print(
            f"{label:14s} n={matrix.shape[0]:3d}  "
            f"CHIRAL R={chiral['R']:.3f} err={chiral['mean_abs_err_h']:.2f}h "
            f"(iters={chiral['iterations']})  |  "
            f"clrcycle R={clr.phase_resultant:.3f} err={clr.mean_abs_err_h:.2f}h"
        )

    table = pd.DataFrame(rows)
    table.to_csv(OUT / "chiral_vs_clrcycle.csv", index=False)

    summary = pd.DataFrame(
        {
            "method": ["CHIRAL", "clrcycle"],
            "mean_R": [table.chiral_R.mean(), table.clrcycle_R.mean()],
            "mean_abs_err_h": [table.chiral_mean_abs_err_h.mean(), table.clrcycle_mean_abs_err_h.mean()],
            "tissues_with_lower_error": [
                int((table.chiral_mean_abs_err_h < table.clrcycle_mean_abs_err_h).sum()),
                int((table.clrcycle_mean_abs_err_h < table.chiral_mean_abs_err_h).sum()),
            ],
        }
    )
    summary.to_csv(OUT / "summary.csv", index=False)
    print("\n" + summary.to_string(index=False))

    r_version = subprocess.run(
        ["Rscript", "-e", "cat(as.character(getRversion()))"],
        capture_output=True, text=True, check=True,
    ).stdout.strip()
    chiral_version = subprocess.run(
        ["Rscript", "-e", 'cat(as.character(packageVersion("CHIRAL")))'],
        capture_output=True, text=True, check=True,
    ).stdout.strip()
    write_versions(OUT / "software_versions.json")
    versions = json.loads((OUT / "software_versions.json").read_text())
    versions["r"] = r_version
    versions["CHIRAL"] = chiral_version
    (OUT / "software_versions.json").write_text(json.dumps(versions, indent=2) + "\n")


if __name__ == "__main__":
    main()
