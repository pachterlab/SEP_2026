# Figure-analysis scripts

This directory contains the complete analysis and rendering workflow for the
clrcycle paper. Run commands from the repository root after activating the
`sep_2026` Conda environment:

```bash
cd /path/to/SEP_2026
conda activate sep_2026
```

## Reproduction order

| Order | Command | Purpose | Main outputs |
|---:|---|---|---|
| 1 | `python scripts/download_data.py` | Download public source data and write provenance records | `data/raw/`, `data/processed/download_manifest.json` |
| 2 | `python scripts/methods_figure.py` | Render the conceptual clrcycle workflow | `figures/Methods_Figure.*` |
| 3 | `python scripts/fig1_circadian.py` | Analyze GSE54650 circadian tissues and render Figure 1 | `results/fig1/`, `figures/Figure1.*` |
| 4 | `python scripts/fig2_fucci.py --plot` | Analyze GSE121265 FUCCI cells, benchmark embeddings, and render Figure 2 | `results/fig2/`, `figures/Figure2.*` |
| 5 | `python scripts/fig3_endometrium.py` | Run endometrium discovery/transfer analysis and render Figure 3 | `results/fig3/`, `figures/Figure3.*` |
| 6 | `python scripts/supplementary_feature_weights.py` | Collate native feature weights and render the liver CT10 heatmap | `results/supp_fig1/`, `results/supplementary_data/`, `figures/supp_fig1_liver_CT10_feature_weights.*` |
| 7 | `python scripts/supplementary_figure2_benchmarks.py --plot` | Benchmark circular ordering and harmonic distortion | `results/supp_fig2/`, `figures/supp_fig2_ordering_benchmarks.*` |
| 8 | `python scripts/supplementary_chiral_comparison.py` | Benchmark CHIRAL on the exact current Figure 1 probe panels using shared R and error scoring | `results/chiral_comparison/` |
| 9 | `python scripts/build_supplement.py` | Assemble the supplementary-figure PDF | `results/supplementary_data/supplement.pdf` |

`make all` performs the complete sequence, including downloads. `make figures`
assumes the public inputs are already present.
`make supplements` runs steps 6, 7, and 9 from existing main-analysis results.
Step 6 also packages the saved comparison in the separate
`Supplementary_Data_4_CHIRAL_vs_clrcycle.xlsx`, with `Summary`,
`Tissue_comparison`, `Sample_phases`, and `Methods` tabs. Run step 8 first only
when the saved CHIRAL comparison needs refreshing. To rebuild just this workbook,
run `python scripts/build_chiral_workbook.py`.

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

The top-level scripts are executable entry points. Code shared across figures
lives in `scripts/lib/`:

- `lib/circadian.py` parses and analyzes the GSE54650 tissue atlas.
- `lib/fucci.py` performs the FUCCI clrcycle fit and embedding benchmark.
- `lib/endometrium.py` performs RNA-seq discovery and frozen cross-platform
  transfer.
- `lib/common.py` contains general parsing, projection, circular-statistics,
  styling, and provenance helpers.

The sample panels use the core package's raw-coordinate ellipse fit. Ellipse
curves and vertices undergo only the existing plotting rotation/reflection.
Each analysis writes `results/fig1/ellipse_geometry.csv`,
`results/fig2/ellipse_geometry.csv`, or `results/fig3/training_ellipse_geometry.csv`.
Circadian also saves separate liver and atlas tables in its existing subfolders.
The first three supplementary workbooks group feature results by analysis: circadian
(GSE54650), FUCCI (GSE121265), and endometrium (GSE234354/GSE111976).
Each has an `Ellipse_features` tab with vertex coordinates and metadata.
Circadian feature weights have one tab per tissue; FUCCI and endometrium each
have a single `Feature_weights` tab. Combined ellipse CSV files are also saved.
Endometrium ellipses are computed only for the training projection in panels A and C.
Held-out and transfer panels have no ellipse fits. Anisotropy p-values use
1,000 radius-preserving random-angle draws with seed 17, without multiplicity
adjustment; leave-one-out columns report robustness, not significance.

Each analysis follows the same broad sequence:

1. Validate and parse source data.
2. Construct the analysis matrix and metadata tables.
3. Fit clrcycle or a comparison method using fixed seeds.
4. Write numerical results before rendering figures.
5. Render PDF, SVG, and high-resolution PNG outputs.
6. Record software versions alongside the numerical results.

## Inputs and generated files

Downloaded source files belong under `data/raw/<accession>/`. Reusable prepared
tables belong under `data/processed/`. Numerical analysis artifacts are written
under `results/<figure>/`, while final figures are written under `figures/`.

The scripts resolve these paths relative to the repository root, so they do not
depend on the shell's current working directory. They never modify files under
`data/raw/` after download.

## Fast reruns

- `python scripts/fig1_circadian.py --plot-only` redraws Figure 1 from existing
  result tables without refitting.
- `python scripts/fig2_fucci.py` runs the numerical analysis and prints the
  benchmark table; add `--plot` to render Figure 2.
- `python scripts/supplementary_figure2_benchmarks.py` writes numerical
  benchmark outputs; add `--plot` to render the supplementary figure.
- `python scripts/supplementary_feature_weights.py` collates the native
  per-sample feature weights and renders Supplementary Figure 1.

These modes are useful when adjusting presentation while preserving validated
numerical artifacts.

## Reproducibility notes

- Stochastic routines use fixed seeds defined in the relevant analysis module.
- Feature selection and any post-fit alignment are documented separately in
  the code because they answer different scientific questions.
- CSV and JSON outputs are the audit trail for each figure; plotting code reads
  those artifacts whenever a plot-only mode is available.
- Some FUCCI and endometrium preprocessing steps call `Rscript`, so the Conda
  environment and its R dependencies must be active.
- `supplementary_chiral_comparison.py` also calls `Rscript`, and additionally
  requires the CHIRAL R package, which is not on CRAN/conda-forge:
  `Rscript -e 'devtools::install_github("naef-lab/CHIRAL/Pkg/CHIRAL")'`.
- The public-data downloader records file sizes and SHA-256 checksums in
  `data/processed/download_manifest.json`.

For dataset citations and the high-level reproduction commands, see the
repository-level `README.md`.
