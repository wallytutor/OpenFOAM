# Simple Geometries Database

## 🔨 Requirements

- Python managed by uv
- OpenFOAM v13

## 🤷‍♂️ Usage

- Run the desired geometry script:

```bash
uv run python wedge-straight-pipe.py

# ... or interactively if modifying ...
uv run ipython -i wedge-straight-pipe.py
```

- Call `gmshToFoam` for the generated mesh:

```bash
gmshToFoam wedge-straight-pipe.msh
```

- Verify mesh generation with `checkMesh`:

```bash
checkMesh
```

> Some specific geometries have an associated shell script with additional details regarding proper conversion into an OpenFOAM polyMesh. Please refer to the respective script for more information.
