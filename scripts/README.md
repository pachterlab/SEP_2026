# Analysis, figure, and workbook scripts

This directory contains the complete analysis and rendering workflow for the
clrcycle paper. Run commands from the repository root after activating the
`sep_2026` Conda environment:

```bash
cd /path/to/SEP_2026
make system-deps  # Debian/Ubuntu; see the root README for other systems
make environment
conda activate sep_2026
make check-env
```

## Reproduction order

| Order | Command | Purpose | Main outputs |
|---:|---|---|---|
| 1 | `python scripts/download_data.py` | Download public source data and write provenance records | `data/raw/`, `data/processed/download_manifest.json` |
| 2 | `python scripts/fig1_architecture.py` | Render the clrcycle architecture with LaTeX math | `figures/Figure1.*` |
| 3 | `python scripts/fig2_circadian.py` | Analyze GSE54650 circadian tissues and render Figure 2 | `results/fig1/`, `figures/Figure2.*` |
| 4 | `python scripts/fig3_fucci.py --plot` | Analyze GSE121265 FUCCI cells, benchmark embeddings, and render Figure 3 | `results/fig2/`, `figures/Figure3.*` |
| 5 | `python scripts/fig4_endometrium.py` | Run endometrium discovery/transfer analysis and render Figure 4 | `results/fig3/`, `figures/Figure4.*` |
| 6 | `python scripts/supplementary_feature_weights.py` | Collate native feature weights and render the liver CT10 heatmap | `results/supp_fig1/`, `results/supplementary_data/`, `figures/supp_fig1_liver_CT10_feature_weights.*` |
| 7 | `python scripts/supplementary_figure2_benchmarks.py --plot` | Benchmark circular ordering and harmonic distortion | `results/supp_fig2/`, `figures/supp_fig2_ordering_benchmarks.*` |
| 8 | `python scripts/supplementary_chiral_comparison.py` | Benchmark CHIRAL on the exact current Figure 2 probe panels using shared R and error scoring | `results/chiral_comparison/` |
| 9 | `python scripts/supplementary_workbooks.py` | Package saved results into four Excel workbooks | `results/supplementary_data/Supplementary_Data_*.xlsx` |

`make all` downloads inputs and runs the main and supplementary figure workflow
(steps 1–7). `make figures` skips downloads.
`make supplements` runs steps 6 and 7 from existing main-analysis results.
These targets do not write Excel workbooks. Step 6 collates the feature-weight
and ellipse CSVs; `make chiral` runs step 8 separately when its saved results
are missing or need refreshing. `make workbooks` runs step 9 without rerunning
analyses or figures.

```bash
make workbooks                       # all four, using saved tables
make workbooks WORKBOOKS="1 2 3"       # no CHIRAL results required
python scripts/supplementary_workbooks.py --workbooks 4
```

The workbook command checks selected inputs before writing outputs. Workbooks
1–3 read `ellipse_geometry.csv` and the annotated `*_feature_weights.csv.gz`
files in `results/supplementary_data/`. Workbook 4 reads `chiral_vs_clrcycle.csv`,
`sample_phases.csv`, `summary.csv`, and `comparison_audit.json` from
`results/chiral_comparison/`. Use `--out-dir PATH` to export elsewhere.

The CHIRAL comparison reads the exact probe IDs and sample IDs from
`results/fig1/all_tissues/` and checks that their expression values reconstruct
the saved clrcycle coordinates. Both methods use `lib.circadian.align`, giving
the same phase-resultant R and mean absolute circular phase error in hours.
This supersedes the earlier 96-gene comparison with mismatched R definitions.
`sample_phases.csv` contains raw phases, aligned phases, and signed errors;
`comparison_audit.json` records the scoring definitions and input checks.
Use `python scripts/supplementary_chiral_comparison.py --rescore-only` to rescore
saved CHIRAL phases after validating their input fingerprints. The upstream
probe selection uses collection times; fitting CHIRAL does not. This single-seed
comparison does not assess robustness across initializations.


## Script structure

The top-level scripts are executable entry points. Analysis and export helpers
live in `scripts/lib/`:

- `lib/circadian.py` parses and analyzes the GSE54650 tissue atlas.
- `lib/fucci.py` performs the FUCCI clrcycle fit and embedding benchmark.
- `lib/endometrium.py` performs RNA-seq discovery and frozen cross-platform
  transfer.
- `lib/common.py` contains general parsing, projection, circular-statistics,
  styling, provenance, and portable R-execution helpers.
- `lib/workbooks.py` contains Excel writers, uniform table formatting, and
  linked contents summaries, used by `supplementary_workbooks.py`.

The sample panels use the core package's raw-coordinate ellipse fit. Ellipse
curves and vertices undergo only the existing plotting rotation/reflection.
Each analysis writes `results/fig1/ellipse_geometry.csv`,
`results/fig2/ellipse_geometry.csv`, or `results/fig3/training_ellipse_geometry.csv`.
Circadian also saves separate liver and atlas tables in its existing subfolders.
The first three supplementary workbooks group feature results by analysis: circadian
(GSE54650), FUCCI (GSE121265), and endometrium (GSE234354/GSE111976).
All four workbooks start with a brief `Summary` contents tab and use matching
formatting. The first three have an `Ellipse_features` tab with vertex coordinates and metadata.
Circadian feature weights have one tab per tissue; FUCCI and endometrium each
have a single `Feature_weights` tab. Combined ellipse CSV files are also saved.
Workbook 4 has `Overall_scores`, `Tissue_comparison`, `Sample_phases`, and `Methods` tabs.
Endometrium ellipses are computed only for the training projection in panels A and C.
Held-out and transfer panels have no ellipse fits. Anisotropy p-values use
1,000 radius-preserving random-angle draws with seed 17, without multiplicity
adjustment; leave-one-out columns report robustness, not significance.
Figure 4's gene-circle annotations are generated from training expression and
histology after the gene order is frozen; no cached annotation table is required.

Each analysis follows the same broad sequence:

1. Validate and parse source data.
2. Construct the analysis matrix and metadata tables.
3. Fit clrcycle or a comparison method using fixed seeds.
4. Write numerical results before rendering figures.
5. Render PDF, SVG, and high-resolution PNG outputs.
6. Record software versions alongside the numerical results.

## Inputs and generated files

Main figure scripts and exports use numbers 1–4 (architecture, circadian,
FUCCI, endometrium). Numerical results retain their historical directories:
`results/fig1/` is circadian, `results/fig2/` is FUCCI, and `results/fig3/` is
endometrium. Supplementary figure numbering is unchanged.

Downloaded source files belong under `data/raw/<accession>/`. Reusable prepared
tables belong under `data/processed/`. Numerical analysis artifacts are written
under `results/<figure>/`, while final figures are written under `figures/`.

The scripts resolve these paths relative to the repository root, so they do not
depend on the shell's current working directory. They never modify files under
`data/raw/` after download.

## Fast reruns

- `python scripts/fig2_circadian.py --plot-only` redraws Figure 2 from existing
  result tables without refitting.
- `python scripts/fig3_fucci.py` runs the numerical analysis and prints the
  benchmark table; add `--plot` to render Figure 3.
- `python scripts/supplementary_figure2_benchmarks.py` writes numerical
  benchmark outputs; add `--plot` to render the supplementary figure.
- `python scripts/supplementary_feature_weights.py` collates the native
  per-sample feature weights and renders Supplementary Figure 1, without writing Excel files.
- `python scripts/supplementary_workbooks.py` rebuilds the workbooks from saved
  tables only; `--workbooks 1 2 3` or `--workbooks 4` selects a subset.

These modes are useful when adjusting presentation while preserving validated
numerical artifacts.

## Reproducibility notes

- Non-math text requires genuine Helvetica through `lib.figure_style`. Install
  it on the system or supply font files in `fonts/Helvetica/` (see
  `fonts/README.md`). Font substitution is disabled; equation rendering is unchanged.
- Stochastic routines use fixed seeds defined in the relevant analysis module.
- Feature selection and any post-fit alignment are documented separately in
  the code because they answer different scientific questions.
- CSV and JSON outputs are the audit trail for each figure; plotting code reads
  those artifacts whenever a plot-only mode is available.
- Some FUCCI and endometrium preprocessing steps call `Rscript`, so the Conda
  environment and its R dependencies must be active.
- `make environment` installs R, Matrix, and Biobase through Conda, then installs
  the pinned CHIRAL revision with `scripts/install_chiral.R` into that R library.
  R 3.6.3 and reference BLAS/LAPACK match the original CHIRAL runtime; Python
  3.11.8 is pinned for compatibility with those libraries.
  `make check-env` checks imports, the environment-local R packages, a small
  CHIRAL fit, Excel export, and the actual methods projection panel rendered
  through LaTeX to PDF/SVG/PNG. It rejects a fallback to system R.
- LaTeX is required for the methods figure. `make system-deps` installs the
  Debian/Ubuntu packages listed in `system-requirements.txt`, including
  `dvipng`, CM-Super fonts, and the extra LaTeX packages. Other platforms need
  equivalent tools on `PATH`; see the root README. The Conda environment check
  reports missing tools/fonts before checking the remaining dependencies.
- The optional combined supplement-document builder and its outputs are ignored
  and excluded from Make targets.
- The public-data downloader records file sizes and SHA-256 checksums in
  `data/processed/download_manifest.json`.

For dataset citations and the high-level reproduction commands, see the
repository-level `README.md`.
