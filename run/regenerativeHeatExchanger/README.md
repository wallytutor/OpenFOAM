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

---

## Generating the report

```bash
export QUARTO_PYTHON="${PWD}/.venv/bin/python"

quarto install tinytex
quarto render tutorial.qmd
```

---

## Running the models

Each model has its own pair of `Allrun` and `Allclean` scripts, which can be used to run the case and clean the generated files, respectively. The `Allrun` script is organized in a way that allows you to easily switch between different modes of operation (e.g., starting from scratch vs. restarting from a previous case, including buoyancy effects or not). You can edit the script to set the desired mode before running it.
