# Case-specific instructions

**Goal:** establish a fluid-only case for initialization of flow pattern.

## Preliminary mesh generation

```bash
python3 model.py

export NUM_PROCS=4

gmshToFoam mesh.msh

createPatch

renumberMesh

splitMeshRegions -cellZones

rm -rf constant/polyMesh
```

## Transform into fluid-only case

```bash
mv constant/fluid/polyMesh constant/polyMesh

rm -rf constant/fluid constant/solid constant/cellToRegion
rm -rf constant/polyMesh/*RegionAddressing

fmesh='constant/polyMesh/boundary'
where='entry0/fluid_to_solid'

foamDictionary $fmesh -entry $where/type -set wall
foamDictionary $fmesh -entry $where/inGroups -remove
foamDictionary $fmesh -entry $where/neighbourRegion -remove
foamDictionary $fmesh -entry $where/neighbourPatch -remove

checkMesh
```

## Prepare conditions and run

Clean automatically generated initial conditions and replace by actual ones:

```bash
rm -rf 0.00000000e+00
cp -avr 0.orig 0.00000000e+00

# XXX just for test, then stop and run in parallel
foamRun
```

Optionally, initialize velocity field by potentialFoam:

```bash
potentialFoam -initialiseUBCs
```

Create handle for ParaView and run in parallel:

```bash
touch case.foam

decomposePar
rm -rf 0.00000000e+00

mpiexec -n $NUM_PROCS foamRun -parallel > log.foamMultiRun &

reconstructPar -latestTime

# XXX only if VTK rendering is required/bugs with PyVista
foamToVTK -latestTime
```