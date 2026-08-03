# Simple Geometries Database

## 🔨 Requirements

- Python managed by uv
- OpenFOAM v13

## 🤷‍♂️ Usage

- Run the desired geometry script:

```bash
uv run wedge-straight-pipe.py
```

- Call `gmshToFoam` for the generated mesh:

```bash
gmshToFoam wedge-straight-pipe.msh
```

- Verify mesh generation with `checkMesh`:

```bash
checkMesh
```
