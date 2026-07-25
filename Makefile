PYTHON ?= python3

.PHONY: setup pipeline dashboard

setup:
	@echo "No external dependencies required."

pipeline:
	$(PYTHON) load_data.py
	$(PYTHON) analysis.py
	$(PYTHON) statistical_analysis.py
	$(PYTHON) subset_analysis.py
	$(PYTHON) generate_dashboard_data.py

dashboard:
	$(PYTHON) -m http.server 8000 --directory dashboard
