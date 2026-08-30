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
| Supplementary Figure 1 | Circular-ordering and harmonic-distortion benchmarks | `scripts/supplementary_figure1_benchmarks.py` |

Figures are written to `figures/` as PDF, SVG, and 300-dpi PNG files. Numerical outputs are written to the corresponding directory under `results/`.

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

The paper version of `clrcycle` is included under `vendor/clrcycle`.

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

To download the data first and then run every analysis:

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

# Supplementary Figure 1: ordering and distortion benchmarks
python scripts/supplementary_figure1_benchmarks.py --plot
```

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
│   └── supp_fig1/           # benchmark outputs
├── scripts/
│   ├── lib/                 # shared analysis functions
│   ├── methods_figure.py
│   ├── fig1_circadian.py
│   ├── fig2_fucci.py
│   ├── fig3_endometrium.py
│   └── supplementary_figure1_benchmarks.py
├── vendor/clrcycle/         # paper version of clrcycle
├── environment.yml
└── Makefile
```

## Reproducibility

Run all commands from the repository root with the `sep_2026` environment active. Stochastic analyses use fixed random seeds, and scripts write software-version records alongside their numerical outputs.

## Citation

Citation information and a preprint link will be added when available.
