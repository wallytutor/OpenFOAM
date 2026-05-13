# Case-specific instructions

**Goal:** full-coupling of fluid-solid conjugate heat transfer.

## Mesh generation

```bash
python3 model.py

export NUM_PROCS=4

gmshToFoam mesh.msh

createPatch

renumberMesh

splitMeshRegions -cellZones

rm -rf constant/polyMesh

checkMesh -region fluid
checkMesh -region solid
```