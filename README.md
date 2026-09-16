# clrcycle paper analyses

This repository contains the data-processing, analysis, benchmarking, and figure-generation code for the clrcycle paper.

> **Preprint:** link coming soon.

## Figures

| Output | Analysis | Script |
|---|---|---|
| Methods Figure | clrcycle workflow | `scripts/methods_figure.py` |
| Figure 1 | Circadian atlas analysis | `scripts/fig1_circadian.py` |
| Figure 2 | FUCCI cell-cycle analysis | `scripts/fig2_fucci.py` |
| Figure 3 | Endometrium menstrual-cycle analysis | `scripts/fig3_endometrium.py` |
| Supplementary Figure 1 | Liver CT10 native feature-weight decomposition | `scripts/supplementary_feature_weights.py` |
| Supplementary Figure 2 | Circular-ordering and harmonic-distortion benchmarks | `scripts/supplementary_figure2_benchmarks.py` |

Figures are written to `figures/` as PDF, SVG, and 300-dpi PNG files. Numerical outputs are written to the corresponding directory under `results/`.

`scripts/supplementary_chiral_comparison.py` is a standalone benchmark (not part of `make figures`/`make all`): it fits CHIRAL (Naef lab) on the exact 240-probe panels and samples from the current Figure 1 (`results/fig1/all_tissues/`), without gene-symbol deduplication. It verifies that the shared expression entries reconstruct the saved clrcycle coordinates, then scores both methods with Figure 1's phase-resultant R and mean absolute phase error in hours using the same rotation/reflection alignment. It requires the CHIRAL R package in addition to the Conda environment; see `scripts/README.md` for the install command. Output: `results/chiral_comparison/`.

The comparison replaces the earlier 96-gene benchmark, which used different R
definitions for the two methods. Raw and aligned phases and input fingerprints
are saved for audit; `--rescore-only` reuses saved CHIRAL phases with matching
inputs. Figure 1's probe selection uses collection times, so this comparison is
conditional on those selected panels. A single CHIRAL seed is evaluated; the
table does not establish robustness to initialization or a failure mechanism.

Four submission workbooks are written to `results/supplementary_data/`:

- `Supplementary_Data_1_GSE54650.xlsx`: circadian tissues.
- `Supplementary_Data_2_GSE121265.xlsx`: FUCCI cells.
- `Supplementary_Data_3_GSE234354_GSE111976.xlsx`: endometrium training biopsies.
- `Supplementary_Data_4_CHIRAL_vs_clrcycle.xlsx`: CHIRAL benchmark with `Overall_scores`, `Tissue_comparison`, `Sample_phases`, and `Methods` tabs.

All four use the same formatting and start with a brief, linked `Summary` contents tab.
`scripts/supplementary_workbooks.py` builds them from saved tables, separately
from analysis and figure generation. Formatting and Excel writers live in
`scripts/lib/workbooks.py`.

The first three workbooks have an `Ellipse_features` tab containing ellipse geometry,
diagnostics, and all four vertices with their metadata. Supplement 1 has a
separate feature-weight tab for each tissue; Supplements 2 and 3 each have one
`Feature_weights` tab. These contain complete sample-by-feature weights and annotations. Column-header
comments explain key quantities. Endometrium ellipse rows and feature weights describe the training biopsies only.
Supplement 4 compares both methods across the 12 circadian tissues, with
overall R, phase error, and tissue wins in `Overall_scores`.
It reads saved tables in `results/chiral_comparison/`; run `make chiral`
first if those tables are missing or need refreshing.
Major-vertex timing is reported only when the anisotropy test has $p<0.05$;
raw fitted coordinates remain available for diagnostics.
The corresponding compressed feature-weight CSVs and combined
`ellipse_geometry.csv` and `ellipse_vertices.csv` remain available.
Per-analysis ellipse tables are saved as `results/fig1/ellipse_geometry.csv`,
`results/fig2/ellipse_geometry.csv`, and `results/fig3/training_ellipse_geometry.csv`.

See [`scripts/README.md`](scripts/README.md) for the execution order, inputs,
outputs, fast rerun modes, shared-library map, and reproducibility notes for
each analysis.

## Installation

```bash
git clone <repository-url>
cd SEP_2026
conda env create -f environment.yml
conda activate sep_2026
```

## Download and prepare data

From the repository root, run:

```bash
python scripts/download_data.py
```

Public inputs are stored under `data/raw/`; prepared analysis inputs are stored under `data/processed/`.

## Reproduce all analyses and figures

```bash
make figures
```

To download the data first and then run the figure workflow:

```bash
make all
```

## Reproduce individual figures

```bash
# Methods Figure: clrcycle workflow
python scripts/methods_figure.py

# Figure 1: circadian atlas
python scripts/fig1_circadian.py

# Figure 2: FUCCI cell cycle
python scripts/fig2_fucci.py --plot

# Figure 3: endometrium menstrual cycle
python scripts/fig3_endometrium.py

# Supplementary Figure 1 and collated feature-weight data
python scripts/supplementary_feature_weights.py

# Supplementary Figure 2: ordering and distortion benchmarks
python scripts/supplementary_figure2_benchmarks.py --plot

# Optional: benchmark clrcycle against CHIRAL on the Figure 1 gene panels
python scripts/supplementary_chiral_comparison.py

# Assemble the supplement document
python scripts/build_supplement.py
```

`make supplements` rebuilds the supplementary figures, collated CSVs,
and supplement PDF from existing main-analysis results. It does not build workbooks.

## Generate supplementary workbooks

After `make supplements` has collated the feature-weight and ellipse tables:

```bash
make workbooks WORKBOOKS="1 2 3"  # feature-weight workbooks only
make chiral                     # generate/refresh the separate CHIRAL results
make workbooks                  # all four workbooks from saved tables
```

Equivalent direct command: `python scripts/supplementary_workbooks.py`.
Use `--workbooks 4` to rebuild only workbook 4 or `--out-dir PATH` to export elsewhere.
Missing inputs are reported before any workbook is written. `make figures` and
`make all` generate analyses and figures; workbook export is a separate step.

Figure 1 can also be redrawn from existing numerical results without rerunning its analysis:

```bash
python scripts/fig1_circadian.py --plot-only
```

## Repository structure

```text
SEP_2026/
├── data/
│   ├── raw/                 # downloaded source datasets
│   └── processed/           # prepared analysis inputs
├── figures/                 # publication-ready figure files
├── results/
│   ├── methods_figure/
│   ├── fig1/                # circadian outputs
│   ├── fig2/                # FUCCI outputs
│   ├── fig3/                # endometrium outputs
│   ├── supplementary_data/  # collated tables and Excel workbooks
│   ├── supp_fig1/           # CT10 heatmap data
│   ├── supp_fig2/           # benchmark outputs
│   └── chiral_comparison/   # CHIRAL vs clrcycle benchmark
├── scripts/
│   ├── lib/                 # shared analysis functions
│   ├── methods_figure.py
│   ├── fig1_circadian.py
│   ├── fig2_fucci.py
│   ├── fig3_endometrium.py
│   ├── supplementary_feature_weights.py
│   ├── supplementary_figure2_benchmarks.py
│   ├── supplementary_chiral_comparison.py
│   └── supplementary_workbooks.py
├── environment.yml
└── Makefile
```

## Reproducibility

Run all commands from the repository root with the `sep_2026` environment active. Stochastic analyses use fixed random seeds, and scripts write software-version records alongside their numerical outputs.

## Citation

Citation information and a preprint link will be added when available.
