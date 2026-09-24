.PHONY: data figures supplements workbooks chiral environment system-deps check-env all

PYTHON := python
WORKBOOKS := 1 2 3 4
CONDA := conda
ENV_NAME := sep_2026

data:
	$(PYTHON) scripts/download_data.py

figures:
	$(PYTHON) scripts/fig1_architecture.py
	$(PYTHON) scripts/fig2_circadian.py
	$(PYTHON) scripts/fig3_fucci.py --plot
	$(PYTHON) scripts/fig4_endometrium.py
	$(MAKE) supplements PYTHON="$(PYTHON)"

supplements:
	$(PYTHON) scripts/supplementary_feature_weights.py
	$(PYTHON) scripts/supplementary_figure2_benchmarks.py --plot

workbooks:
	$(PYTHON) scripts/supplementary_workbooks.py --workbooks $(WORKBOOKS)

chiral:
	$(PYTHON) scripts/supplementary_chiral_comparison.py

# Debian/Ubuntu only; run once before creating/updating the Conda environment.
system-deps:
	sudo apt-get update
	sudo apt-get install -y $(shell cat system-requirements.txt)

environment:
	$(CONDA) env update --name $(ENV_NAME) --file environment.yml
	$(CONDA) run --no-capture-output --name $(ENV_NAME) Rscript scripts/install_chiral.R
	$(CONDA) run --no-capture-output --name $(ENV_NAME) python scripts/check_environment.py

check-env:
	$(PYTHON) scripts/check_environment.py

all: data
	$(MAKE) figures PYTHON="$(PYTHON)"
