# Welcome to IGA-Python

Contains the source code of the [IGA-Python](https://pyccel.github.io/IGA-Python/) site.

## Editing and building IGA-Python locally

1. Clone this repository and then install the required dependencies.

```bash
git clone https://github.com/pyccel/IGA-Python.git
cd IGA-Python

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

3. Edit the desired Python notebooks (`*.ipynb`) and/or Markdown files (`*.md`). Check the [MyST syntax cheat sheet](https://jupyterbook.org/en/stable/reference/cheatsheet.html) for reference.
4. Build the docs by running `make`.
5. View your changes on the browser: `open _build/html/index.html`
