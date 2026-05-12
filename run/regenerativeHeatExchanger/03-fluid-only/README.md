# Case-specific instructions

**Goal:** include turbulence in initialized fluid-only case.

## Case preparation

```bash
export NUM_PROCS=4

# Copy mesh from previous case:
cp -avr ../01-fluid-only/constant/polyMesh constant/polyMesh

# Create initial conditions by one of the following options:
mkdir -p 0.00000000e+00

# ... if copying, remove the time directory under uniform:
cp -avr ../02-fluid-only/2.00000000e+00 0.00000000e+00
rm -rf 0.00000000e+00/uniform/

# If a directory was created, expand original conditions to actual mesh:
foamDictionary 0.orig/p      -expand > 0.00000000e+00/p
foamDictionary 0.orig/T      -expand > 0.00000000e+00/T
foamDictionary 0.orig/U      -expand > 0.00000000e+00/U
foamDictionary 0.orig/alphat -expand > 0.00000000e+00/alphat
foamDictionary 0.orig/nut    -expand > 0.00000000e+00/nut
foamDictionary 0.orig/k      -expand > 0.00000000e+00/k
foamDictionary 0.orig/omega  -expand > 0.00000000e+00/omega

# Otherwise only the buoyancy related fields:
foamDictionary 0.orig/p_rgh  -expand > 0.00000000e+00/p_rgh

# XXX just for test, then stop and run in parallel
foamRun
```

## Run simulation

```bash
decomposePar

mpiexec -n $NUM_PROCS foamRun -parallel > log.foamMultiRun &

reconstructPar -latestTime
```