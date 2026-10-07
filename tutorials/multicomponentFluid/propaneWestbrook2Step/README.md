# Propane-Air Westbrook 2-Step Flame Simulation

In this study we setup a non-premixed stoichiometric propane-air flame using the Westbrook 2-step mechanism in a coaxial jet geometry. The project structure is documented below:

- [Allrun](./Allrun): main workflow script, orchestrates running the pre-set simulation cases.

- [Allclean](./Allclean): main clean script, removes all generated files.

- [notebook.sh](./notebook.sh): helper script to run jupyter notebooks and generate reports; run `./notebook.sh -h` for details.

- [project.py](./project.py): implements the pre-/post-processing of the model.

- [domain.py](./domain.py): provides a (partially) parametric domain relying on constants provided in `project.py`.

- [report.ipynb](./report.ipynb): notebook for the generation of the report of the default simulations.

- [watcher.ipynb](./watcher.ipynb): notebook for watching the results of the simulation.

- [chemkin/](./chemkin/): contains the files and scripts used to generated the kinetics model used in `model/`.

- [model/](./model/): template OpenFOAM case; not to be run directly but as a basis for new simulations.

> Before running anything for the first time consider syncing the environment with `uv sync`. After running the simulations, and only after that, generate the report with `./notebook.sh -report`
