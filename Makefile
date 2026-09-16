.PHONY: data figures supplements workbooks chiral all

PYTHON := python
WORKBOOKS := 1 2 3 4

data:
	$(PYTHON) scripts/download_data.py

figures:
	$(PYTHON) scripts/methods_figure.py
	$(PYTHON) scripts/fig1_circadian.py
	$(PYTHON) scripts/fig2_fucci.py --plot
	$(PYTHON) scripts/fig3_endometrium.py
	$(MAKE) supplements PYTHON="$(PYTHON)"

supplements:
	$(PYTHON) scripts/supplementary_feature_weights.py
	$(PYTHON) scripts/supplementary_figure2_benchmarks.py --plot
	$(PYTHON) scripts/build_supplement.py

workbooks:
	$(PYTHON) scripts/supplementary_workbooks.py --workbooks $(WORKBOOKS)

chiral:
	$(PYTHON) scripts/supplementary_chiral_comparison.py

all: data
	$(MAKE) figures PYTHON="$(PYTHON)"
