.PHONY: data figures all

PYTHON := python

data:
	$(PYTHON) scripts/download_data.py

figures:
	$(PYTHON) scripts/methods_figure.py
	$(PYTHON) scripts/fig1_circadian.py
	$(PYTHON) scripts/fig2_fucci.py --plot
	$(PYTHON) scripts/fig3_endometrium.py
	$(PYTHON) scripts/supplementary_feature_weights.py
	$(PYTHON) scripts/supplementary_figure2_benchmarks.py --plot

all: data figures
