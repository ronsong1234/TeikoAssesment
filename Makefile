PYTHON ?= python3
NPM ?= npm

.PHONY: setup pipeline dashboard

setup:
	$(PYTHON) -m pip install -r requirements.txt
	$(NPM) --prefix dashboard ci --ignore-scripts --no-audit --no-fund

pipeline:
	$(PYTHON) load_data.py
	$(PYTHON) analysis.py
	$(PYTHON) statistical_analysis.py
	$(PYTHON) subset_analysis.py
	$(PYTHON) generate_dashboard_data.py

dashboard:
	$(NPM) --prefix dashboard run dev -- --host 0.0.0.0
