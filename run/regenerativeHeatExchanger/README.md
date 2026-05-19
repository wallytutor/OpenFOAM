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

---

## Running with Podman

As an alternative to installing OpenFOAM 13 locally, you can build and run the models using [Podman](https://podman.io/) (or Docker) with the provided `Containerfile`.

### 1. Build the Container Image

Run the following command in this directory to build the minimal Ubuntu 24.04 image with OpenFOAM 13:

```bash
podman build -t openfoam13 -f Containerfile .
```

### 2. Run the Cases from the Current Directory

Mount the root directory of the project to `/workspace` inside the container, *i.e.* run the following command from `OpenFOAM` folder in your host machine, as the project depends on scripts located in this directory:

```bash
podman run --rm -it -v ".:/workspace:Z" openfoam13 bash
```

Using `cd` navigate to this directory and run as in a local environment. Since the container has an entrypoint that automatically sources the OpenFOAM 13 environment, you can run any OpenFOAM or local script directly.
