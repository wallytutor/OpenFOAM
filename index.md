# OpenFOAM Extensions

Curated OpenFOAM 13 directory with additional workflow tools.

## Prerequisites

This respository aims at being homogeneous in the sense it tries to keep a consistent toolset across all cases. Shell automation is done solely with bash and indirectly through Python in some cases. Unless otherwise specified, the following tools are used:

- OpenFOAM v13
- Python managed by uv
- Quarto with Typst support
- Apptainer with Docker/Podman

## Before starting

Many features available in this repository rely on the environment configuration provided by `etc/bashrc`. In general the scripts will source that file for you, but environment configuration will not remain in CLI after the script finishes. To experience a smooth experience, it is recommended to source that file in your current shell from this directory:

```bash
source $PWD/etc/bashrc
```

If you want this to be done automatically for every new bash session, you can append that to your `~/.bashrc` by calling `foam_source_append`. This function only appends the configuration if not already present and tests for the existence of the configuration file so that your shell will not raise error messages if you delete this repository later.

## Compilation methods

It is possible to build the whole extensions library or isolated features, depending on your needs. The following logic is applied throughout this sources directory:

- If a directory provides an `Allwmake` file, it allows to compile the whole tree below it; this is the case for the root directory, for instance.

- If a directory of `src/` provides a subdirectory called `Make/`, then it supports the individual library build which is done by running `wmake libso` from that directory.

All builds are written to `$FOAM_USER_APPBIN` and `$FOAM_USER_LIBBIN`.

The simplest way it to build everything from the repository root:

```bash
./Allwmake
```

> If you encounter failures when loading a library, please consider running `./Allwclean` first to clear the build cache. Sometimes previous build artifacts may cause unexpected behavior.

## Running tests

Some libraries may provide test programs under a `test/` sub-directory, and all of them are structured in the same way. The following instructions indicate how to build, run, and clean those directories.

To build and run the test suite:

```bash
cd test
./Allrun
```

To clean test build artifacts:

```bash
cd test
./Allclean
```

## Running cases

Whenever possible, cases are meant to be run from within an Apptainer instance. The SIF file can be generated with script `Containerfile.sh` located at the root of this repository. It can be later transferred to any HPC environment where Apptainer is supported, ensuring a consisten environment. Assuming you are connected to the target compute node, the following elements indicate how to instantiate an environment, run the simulation, and inspect partial results on demand. Notice that we will make use of the `screen` utility to keep services alive even if one disconnects from the node.

- Connecting to a session:

```bash
# Declare variables to set the case to be run:
APP_NAME="methaneAirSensitivity"
APP_PATH="tutorials/multicomponentFluid/$APP_NAME"

# A SIF generated Containerfile.sh should be named like this:
SIF_FILE="${name}-$(whoami).sif"

# Start a named screen session:
screen -S $APP_NAME

# Start apptainer instance on the background:
apptainer instance start -B $PWD --writable-tmpfs $SIF_FILE $APP_NAME

# Enter running instance in shell mode:
apptainer shell instance://$APP_NAME
```

- From within a session, setup the base environment and run the case(s):

```bash
# Instantiate the environment from instance:
(cd $APP_PATH && uv sync)

# Move into the simulation directory and run:
(cd $APP_PATH && ./Allrun &)

# Ctrl + A, D to detach screen
```

- Inspecting a running simulation (from the same node as above):

```bash
# Attach to the running screen session
screen -r $APP_NAME

# Enter running instance in shell mode:
apptainer shell instance://$APP_NAME

# Move into the simulation directory and check status:
cd $APP_PATH

# Start Jupyter server to connect from case notebook:
uv run jupyter-notebook --no-browser --ip=0.0.0.0 \
    --ServerApp.token='' --ServerApp.password=''
```

- Stop instance after working:

```bash
# Stop instance after working:
apptainer instance stop $APP_NAME
```
