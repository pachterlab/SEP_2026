# clrcycle paper analyses

This repository contains the data-processing, analysis, benchmarking, and figure-generation code for the clrcycle paper.

> **Preprint:** link coming soon.

## Figures

| Output | Analysis | Script |
|---|---|---|
| Figure 1 | clrcycle architecture | `scripts/fig1_architecture.py` |
| Figure 2 | Circadian atlas analysis | `scripts/fig2_circadian.py` |
| Figure 3 | FUCCI cell-cycle analysis | `scripts/fig3_fucci.py` |
| Figure 4 | Endometrium menstrual-cycle analysis | `scripts/fig4_endometrium.py` |
| Supplementary Figure 1 | Liver CT10 native feature-weight decomposition | `scripts/supplementary_feature_weights.py` |
| Supplementary Figure 2 | Circular-ordering and harmonic-distortion benchmarks | `scripts/supplementary_figure2_benchmarks.py` |

Main figure numbering is: 1 architecture, 2 circadian, 3 FUCCI, and 4 endometrium.
Script names match the current figure numbers. The existing numerical directories
`results/fig1/`–`results/fig3/` retain their historical names for circadian, FUCCI,
and endometrium results, respectively, to preserve saved analysis inputs.
Supplementary figures retain their existing numbering.

Figures are written to `figures/` as PDF, SVG, and 300-dpi PNG files. Numerical outputs are written to the corresponding directory under `results/`.

`scripts/supplementary_chiral_comparison.py` is a standalone benchmark (not part of `make figures`/`make all`): it fits CHIRAL (Naef lab) on the exact 240-probe panels and samples from the current Figure 2 (`results/fig1/all_tissues/`), without gene-symbol deduplication. It verifies that the shared expression entries reconstruct the saved clrcycle coordinates, then scores both methods with Figure 2's phase-resultant R and mean absolute phase error in hours using the same rotation/reflection alignment. The project environment setup installs its pinned CHIRAL R package. Output: `results/chiral_comparison/`.

The comparison replaces the earlier 96-gene benchmark, which used different R
definitions for the two methods. Raw and aligned phases and input fingerprints
are saved for audit; `--rescore-only` reuses saved CHIRAL phases with matching
inputs. Figure 2's probe selection uses collection times, so this comparison is
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
make system-deps  # Debian/Ubuntu: install the LaTeX tools and fonts
make environment
conda activate sep_2026
make check-env
```

The setup installs Python dependencies, R, Matrix, Biobase, and a pinned CHIRAL
revision. clrcycle is installed from a pinned Git revision, so no local `vendor/`
folder is required. R 3.6.3 and reference BLAS/LAPACK preserve the original
CHIRAL calculation; Python 3.11.8 permits these libraries in one environment.
`make environment` also updates an existing environment.
Use `make environment ENV_NAME=another_name` for a separate installation.
The methods figure uses actual LaTeX for its mathematical labels. System
requirements are recorded in `system-requirements.txt`; `make system-deps`
installs them on Debian/Ubuntu (administrator access is required). Other
systems need a TeX Live/MacTeX installation providing `latex`, `kpsewhich`,
`type1cm.sty`, and `type1ec.sty` (CM-Super fonts), plus `dvipng` and Ghostscript
on `PATH`. Non-math figure text requires genuine Helvetica. Install it on your
system or place your Helvetica TTF/OTF files in `fonts/Helvetica/`; see
[`fonts/README.md`](fonts/README.md). Font files are supplied separately from
Conda. The scripts report an error if Helvetica is missing; no Nimbus Sans,
Arial, or DejaVu substitution is permitted. These are external tools, so installing the Conda YAML alone is
not sufficient. `make check-env` verifies the tools and fonts and exports the
actual methods projection panel to PDF, SVG, and PNG using LaTeX.
The optional combined supplement-document builder remains local and ignored.

All current presentation changes are implemented in the tracked figure scripts:
consistent first Fourier coordinate labels with `(u)` and `(v)`, gene guide
lines in Figures 1 and 2, separated projection formulas, and the
“clrcycle architecture” title. The standard commands below regenerate them;
no temporary rendering scripts or manual image edits are required.

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
# Figure 1: clrcycle architecture
python scripts/fig1_architecture.py

# Figure 2: circadian atlas
python scripts/fig2_circadian.py

# Figure 3: FUCCI cell cycle
python scripts/fig3_fucci.py --plot

# Figure 4: endometrium menstrual cycle
python scripts/fig4_endometrium.py

# Supplementary Figure 1 and collated feature-weight data
python scripts/supplementary_feature_weights.py

# Supplementary Figure 2: ordering and distortion benchmarks
python scripts/supplementary_figure2_benchmarks.py --plot

# Optional: benchmark clrcycle against CHIRAL on the Figure 2 gene panels
python scripts/supplementary_chiral_comparison.py
```

`make supplements` rebuilds the supplementary figures and collated CSVs
from existing main-analysis results. Workbooks are generated separately.

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

Figure 2 can also be redrawn from existing numerical results without rerunning its analysis:

```bash
python scripts/fig2_circadian.py --plot-only
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
│   ├── fig1_architecture.py
│   ├── fig2_circadian.py
│   ├── fig3_fucci.py
│   ├── fig4_endometrium.py
│   ├── supplementary_feature_weights.py
│   ├── supplementary_figure2_benchmarks.py
│   ├── supplementary_chiral_comparison.py
│   └── supplementary_workbooks.py
├── environment.yml
└── Makefile
```

## Reproducibility

Run all commands from the repository root with the `sep_2026` environment active. Stochastic analyses use fixed random seeds, and scripts write software-version records alongside their numerical outputs.

The complete workflow was tested in a fresh Linux x86_64 environment: 144
numerical tables and all four workbooks matched the saved values exactly.
All nine PNGs were byte-identical; PDF/SVG content matched after excluding
timestamps and generated SVG IDs. Earlier 96-probe and other legacy artifacts
remain historical files; the current Make targets do not regenerate them.

## Citation

Citation information and a preprint link will be added when available.
