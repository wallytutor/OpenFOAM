# Sandia Propane Jet

Reference case for propane jet dilution provided by [Sandia](https://tnfworkshop.org/data-archives/simplejet/propanejet/).

## 🔨 Requirements

- Python managed by uv
- Quarto with Typst support
- OpenFOAM v13

## 🤷‍♂️ Usage

Before running anything for the first time consider syncing the environment:

```bash
uv sync
```

For generating the STL files required for OpenFOAM meshing run:

```bash
uv run ipython model/constant/geometry/domain.py
```

## 📃 Generating the report

```bash
uv run majordome-build-qmd --file docs/tutorial.qmd
```

## 🛰️ Download reference

```powershell
$url = "https://tnfworkshop.org/wp-content/uploads/2019/02"
$file = "SandiaPropaneJet.zip"

Start-BitsTransfer -Source "$url/$file" -Destination "$file"
Expand-Archive -Path "$file" -DestinationPath "reference/"
```

To understand the above, the following are also required:

- [Dibble (1984)](https://doi.org/10.1016/0010-2180(84)90170-6)
- [Dibble (1987)](https://doi.org/10.1007/BF00776180)

## 🛖 Concept and files

### Meshing

An automation script **[Allmesh](model/Allmesh)** is provided for meshing generation. If can be provided a number of cores so that meshing can be run in parallel:

```bash
./Allmesh 4
```

This script orchestrates the full sequence of meshing:

- Cleans the case.

- Runs `surfaceFeatures` and `blockMesh`.

- Decomposes the domain if parallel cores (4 or 20) are specified.

- Executes `snappyHexMesh` (in parallel or sequentially).

- Reconstructs the mesh using `reconstructPar -constant`.

- Runs `topoSet` and `createPatch -overwrite` to split the inlet boundaries.

- Runs `renumberMesh` to reduce bandwidth.

- Runs `checkMesh` to verify the mesh.

---

The following files are used for preliminary meshing:

- **[blockMeshDict](model/system/blockMeshDict)**: defines the background grid (required as starting point for `snappyHexMesh`) and pre-defines patches.

- **[surfaceFeaturesDict](model/system/surfaceFeaturesDict)**: extracts feature edges from all STL files for explicit snapping.

- **[meshQualityDict](/model/system/meshQualityDict)**: sources the default values for mesh quality controls required by `snappyHexMesh`.

- **[snappyHexMeshDict](/model/system/snappyHexMeshDict)**: the main meshing script, it sources all geometry files and defines the required refinement/boundary layer controls.

---

Because `inletFuel` and `inletCoflow` lie on the outer boundary where `snappyHexMesh` does not perform cuts, they were initially removed as zero-sized. A standard post-meshing splitting phase was introduced:

- **[topoSetDict](model/system/topoSetDict)**: selects `faceSets` on the `inlet` patch using `cylinderToFace` geometry sources to isolate `inletFuelFaces` and `inletCoflowFaces`.

- **[createPatchDict](model/system/createPatchDict)**: reconstructs the new physical patches `inletFuel` and `inletCoflow` from the `faceSets`.

---

After running, one must consider checking the following logs and at least check for the illustrated example outputs:

- **log.snappyHexMesh**: snapping and layer addition phases finished successfully.

  ```plaintext
  Finished meshing without any errors
  Finished meshing in = 191.008 s.
  ```

- **log.checkMesh**: check mesh validity and boundary patches.

  ```plaintext
  Mesh stats
      points:           2867496
      faces:            8268812
      internal faces:   8159386
      cells:            2701398
      faces per cell:   6.08137
      boundary patches: 6

  Checking patch topology for multiply connected surfaces...
                     Patch    Faces   Points                  Surface topology
                    outlet     4080     4201  ok (non-closed singly connected)
                 ductWalls    98150    98640  ok (non-closed singly connected)
               pipeEndWall      112      156  ok (non-closed singly connected)
             pipeInnerWall     1600     1632  ok (non-closed singly connected)
             pipeOuterWall     3200     3576  ok (non-closed singly connected)
                 inletFuel      152      157  ok (non-closed singly connected)
               inletCoflow     2132     2376  ok (non-closed singly connected)

  Checking geometry...
      Overall domain bounding box (-0.15 -0.15 -0.05) (0.15 0.15 2)
      Max cell openness = 3.72574e-16 OK.
      Max aspect ratio = 16.2085 OK.
      Mesh non-orthogonality Max: 54.9018 average: 4.92425
      Non-orthogonality check OK.
      Face pyramids OK.
      Max skewness = 2.69949 OK.
      Mesh OK.
  ```
