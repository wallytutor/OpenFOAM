# Case-specific instructions

**Goal:** establish a fluid-only case for initialization of flow pattern.

## Mesh generation

- Create the mesh for the hole system:

```bash
export NUM_PROCS=4

python3 model.py

gmshToFoam mesh.msh

createPatch

renumberMesh

splitMeshRegions -cellZones

rm -rf constant/polyMesh
```

- Transform into fluid-only case by extracing the fluid region:

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
