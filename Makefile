PYTHON ?= python

.PHONY: test install dev clean

install:
	$(PYTHON) -m pip install -e .

dev: install
	$(PYTHON) -m pip install pytest

test:
	$(PYTHON) -m pytest -q

clean:
	rm -rf .pytest_cache .mypy_cache out out_rnaseq __pycache__ */__pycache__
