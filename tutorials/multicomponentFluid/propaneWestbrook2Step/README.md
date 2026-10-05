# Propane-Air Westbrook 2-Step Flame Simulation

In this study we setup a non-premixed stoichiometric propane-air flame using the Westbrook 2-step mechanism in a coaxial jet geometry.

Before running anything for the first time consider syncing the environment with `uv sync`. The full simulation workflow is orchestrated by the `Allrun` script. It covers meshing (if needed), setting up the simulation, and running it.

After running the simulations, and only after that, generate the report with:

```bash
uv run majordome-build-qmd --file report.ipynb
```