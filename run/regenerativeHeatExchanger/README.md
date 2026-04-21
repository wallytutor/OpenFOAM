# Regenerative heat exchanger

Conjugate heat transfer models for regenerative heat exchangers.

---

## Preparing the environment

- Models developed and tested under [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) running Ubuntu 24.04; if you are working under Windows, you can install the referred distribution with `wsl --install Ubuntu-24.04`.

- Install [OpenFOAM (v13)](https://openfoam.org/download/13-ubuntu/) for your distribution (click link for instructions related to the referenced version).

- Source the bootstrap script `source bootstrap.sh` to set up the environment; alternatively (if you have OpenFOAM and uv installed), you can manually start a virtual environment with `uv`:

```bash
# Create a virtual environment
uv venv --python 3.12 .venv

# Activate the virtual environment
source .venv/bin/activate

# Install the required packages
uv pip install -r requirements.txt
```

## Generating the report

```bash
export QUARTO_PYTHON="${PWD}/.venv/bin/python"

quarto install tinytex
quarto render report/
```

---

## Running the models

Each model has its own *README.md* file which can be read by [Jupytext](https://jupytext.readthedocs.io/en/latest/) within a Jupyter environment to manage the models. Please refer to the respective files for instructions on how to run each model.
