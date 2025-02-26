BUILD_DIR = ./_build
SRC_DIR = .

docs:
	jupyter-book build .

docs-strict:
	jupyter-book build --warningiserror --nitpick --keep-going .

.PHONY: clean

clean:
	jupyter-book clean .

