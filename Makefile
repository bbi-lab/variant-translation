PYTHON ?= $(if $(wildcard .venv/bin/python),.venv/bin/python,python)

.PHONY: install-publish-deps test clean build-dist check-dist release-prep release release-testpypi

install-publish-deps:
	$(PYTHON) -m pip install .[publish]

test:
	$(PYTHON) -m pytest tests/test_reverse_translate_variants.py tests/test_compare_reverse_translated_variants.py -q

clean:
	rm -rf dist build *.egg-info

build-dist:
	$(PYTHON) -m build

check-dist:
	$(PYTHON) -m twine check dist/*

release-prep: test clean build-dist check-dist

release: release-prep
	$(PYTHON) -m twine upload dist/*

release-testpypi: release-prep
	$(PYTHON) -m twine upload --repository testpypi dist/*
