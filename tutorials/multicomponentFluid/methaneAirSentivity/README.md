# Methane-Air Sensitivity Analysis

In this study we compare the effect of air temperature over the flame structure for stoichiometric methane-air mixture.

## 🔨 Requirements

- Python managed by uv
- Quarto with Typst support
- OpenFOAM v13

## 🤷‍♂️ Usage

Before running anything for the first time consider syncing the environment:

```bash
uv sync
```

Generate the base mesh (in interactive mode, if required) by running:

```bash
uv run ipython -i domain.py
```

## 📃 Generating the report

After running the simulations, generate the report with the following:

```bash
uv run majordome-build-qmd --file report.ipynb
```