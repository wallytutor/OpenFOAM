# Methane-Air Sensitivity Analysis

In this study we compare the effect of air temperature over the flame structure for stoichiometric methane-air mixture.

Before running anything for the first time consider syncing the environment:

```bash
uv sync
```

The full simulation workflow is orchestrated by the `Allrun` script. It covers meshing (if needed), setting up the simulation, and running it.

> Generate the base mesh (in interactive mode, if modifying the case) by running `uv run ipython -i domain.py`.

After running the simulations, and only after that, generate the report with the following:

```bash
uv run majordome-build-qmd --file report.ipynb
```