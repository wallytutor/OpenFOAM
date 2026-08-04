# Propane Oxicombustion

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