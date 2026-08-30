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
| 6 | `python scripts/supplementary_figure1_benchmarks.py --plot` | Benchmark circular ordering and harmonic distortion | `results/supp_fig1/`, `figures/Supplementary_Figure1.*` |

`make all` performs the complete sequence, including downloads. `make figures`
assumes the public inputs are already present.

## Script structure

The top-level scripts are executable entry points. Code shared across figures
lives in `scripts/lib/`:

- `lib/circadian.py` parses and analyzes the GSE54650 tissue atlas.
- `lib/fucci.py` performs the FUCCI clrcycle fit and embedding benchmark.
- `lib/endometrium.py` performs RNA-seq discovery and frozen cross-platform
  transfer.
- `lib/common.py` contains general parsing, projection, circular-statistics,
  styling, and provenance helpers.

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
- `python scripts/supplementary_figure1_benchmarks.py` writes numerical
  benchmark outputs; add `--plot` to render the supplementary figure.

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
- The public-data downloader records file sizes and SHA-256 checksums in
  `data/processed/download_manifest.json`.

For dataset citations and the high-level reproduction commands, see the
repository-level `README.md`.
