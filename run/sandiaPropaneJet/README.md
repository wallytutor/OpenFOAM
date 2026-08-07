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

> If updating the mesh is required, please edit `domain.py` and then run `uv run python domain.py` before running `Allrun`. The validated grid is commited within the project for reproducibility.

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
