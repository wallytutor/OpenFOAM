# Cantera container and Apptainer environment

This directory contains the [`Containerfile`](Containerfile) and [`Makefile`](Makefile) for compiling **Cantera 3.2.0** with multithreading, OpenMPI, custom SUNDIALS & Eigen, and a dedicated Python/PyO3/Rust toolchain.

---

## Building the Apptainer image

From this directory (`cantera/`), build the container and convert it into an Apptainer `.sif` image. This is automated with the following, which builds the container via `podman`, exports the archive, and builds the `.sif` image:

```bash
make image
```

This generates `cantera-$(whoami).sif`.

---

## Running with Apptainer instance

You can start the image as a background instance and connect to it. To start the background instance with writable *tmpfs* and current directory bound:

```bash
# Assuming you are one level above the repository directory:
apptainer instance start     \
    -B "$(realpath $PWD)" --writable-tmpfs \
    "OpenFOAM/cantera/cantera-$(whoami).sif"  \
    "cantera"
```

> If working in a symbolic link directory, it is important to use `realpath` instead of `$PWD` directly, as Apptainer needs the actual resolved path to work.

---

### Entering the running instance

Enter the active container instance in interactive shell mode:

```bash
apptainer shell instance://cantera
```

Inside the shell, verify that the environment and Cantera are functional:

```bash
python -c "import cantera as ct; print(f'Cantera version: {ct.__version__}')"
```

To be able to connect a notebook (or VS Code Jupyter extension) to the running instance, you need to start the Jupyter server **on the node** where the instance is running:

```bash
# Start Jupyter server to connect from case notebook:
uv run jupyter-notebook --no-browser --ip=0.0.0.0 \
    --ServerApp.token='' --ServerApp.password=''
```

> It is practical to disable the token if you know what you are doing or working internally (e.g. not exposed to internet).

---

### Stopping the instance

Once finished, exit the shell and stop the running instance:

```bash
apptainer instance stop cantera
```
