# Acetylene Pyrolysis

In this study we analyze the effect of reactor temperature over the decomposition of acetylene under low pressure conditions.

Before running anything for the first time consider syncing the environment:

```bash
uv sync
```

The full simulation workflow is orchestrated by the `Allrun` script.

After running the simulations, and only after that, generate the report with the following:

```bash
uv run majordome-build-qmd --file report.ipynb
```
