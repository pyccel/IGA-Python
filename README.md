# Welcome to IGA-Python

This project provides a tutorial for isogeometric analysis (IGA) using Python and the Psydac library (pyccel/psydac). The numerical examples can be consulted online at pyccel.github.io/IGA-Python, or run with JupyterNotebook on a personal computer.

## Editing and building IGA-Python locally

1. Clone this repository and then install the required dependencies.

```bash
git clone https://github.com/pyccel/IGA-Python.git
cd IGA-Python
IGA_PYTHON_DIR=$(pwd)

# Install dependencies on a virtual environment
python3 -m venv iga-python-env
source iga-python-env/bin/activate
pip3 install -r requirements_ntbk.txt
```

2. Install Psydac. Skip this step if Psydac is already installed.

```bash
# Modify these variables if you're using your own psydac fork/branch
PSYDAC_REMOTE="https://github.com/pyccel/psydac.git"
BRANCH="devel"

# Install psydac
pip install git+${PSYDAC_REMOTE}@${BRANCH}
```

3. Run or modify the desired Python notebooks (`*.ipynb`) and Markdown files (`*.md`). Check the [MyST syntax cheat sheet](https://jupyterbook.org/en/stable/reference/cheatsheet.html) for reference.

> [!NOTE]  
> Before committing changes to git, the notebooks have to be cleaned first. See ["Committing your changes to git"](#committing-your-changes-to-git) for more information.

4. Build the docs by running `make` under the IGA-Python folder. This involves running all `*.ipynb` files in the background to make sure they are functional.

```bash
cd ${IGA_PYTHON_DIR}
make
```

5. View your changes on the browser.

```bash
open ${IGA_PYTHON_DIR}/_build/html/index.html
```

## Committing your changes to git

Running the Python notebooks embeds extra information like the name of virtual environment, Python version used, cell outputs, etc. The notebooks should be free of system- and runtime-specific information before committing them to source control. We suggest the following commit workflow:

1. Run [`nb-clean`](https://github.com/srstevenson/nb-clean) on the modified notebook/s. There are different ways to run this command:

```bash
# Clean a single notebook
nb-clean clean --remove-empty-cells --remove-all-notebook-metadata chapter1/poisson.ipynb

# Clean all '.ipynb' files under IGA-Python
make clean-notebooks

# Automatically run nb-clean on `git add`-ed *.ipynb files
nb-clean add-filter --remove-empty-cells --remove-all-notebook-metadata

# Undo previous command
nb-clean remove-filter
```

2. Create a branch for your local changes, e.g. `git checkout -b my-local-fixes`.
3. Stage the modified files with `git add`. Then run `git diff` to check if the changed `*.ipynb` notebooks doesn't include unnecessary diffs (e.g. notebook metadata).
4. `git commit` your changes.
5. *OPTIONAL*. Share your changes to this repo via a [pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request).
