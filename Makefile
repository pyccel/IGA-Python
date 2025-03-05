BUILD_DIR = ./_build
SRC_DIR = .

.PHONY: docs docs-strict clean clean-notebooks

docs:
	@echo "Building IGA-Python docs..."
	@jupyter-book build .
	@echo "Done."

docs-strict:
	@echo "Building IGA-Python docs with strict rules on..."
	@jupyter-book build --warningiserror --nitpick --keep-going .
	@echo "Done."

clean:
	@echo "Removing previous IGA-Python build artifacts..."
	@jupyter-book clean .
	@echo "Done."

clean-notebooks:
	@echo "Running 'nb-clean --remove-empty-cells --remove-all-notebook-metadata' on all '*.ipynb' files..."
	@find . -type f -iname '*.ipynb' -exec nb-clean clean --remove-empty-cells --remove-all-notebook-metadata {} \+
	@echo "Done."

