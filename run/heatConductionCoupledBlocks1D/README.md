# Heat conduction between two blocks in 1D

- Start the project by sourcing `source bootstrap.sh` in the terminal.

- Run the case with `./Allrun`, this will handle the whole workflow.

- Once finished, run `python3 model.py` for post-processing of the results.

---

## Sample results

![Temperature evolution](evolution.png)

---

## Manual post-processing

Below we provide a snippet for manual post-processing of the results, which can be used as a template for more complex cases. This should be run in the same environment as `model.py` (i.e. with the same dependencies installed).

```python
import pandas as pd
from pathlib import Path
from majordome_utilities.plotting import plot2d

# Select a file from the `postProcessing` directory:
start = "0.00000000e+00"
block = Path("postProcessing/block/")
fname = block / "blockToMediaHeatFlux" / start / "wallHeatFlux.dat"

# This should be able to read OpenFOAM generated tables:
df = pd.read_csv(fname, sep=r"\s+", header=None, comment="#")

# Select the right columns/unit conversions:
t = df[0].to_numpy() / 60
y = df[5].to_numpy() * -1

# Plot the results:
p = plot2d(t, y)
p.axes[0].set_xlabel("time [min]")
p.axes[0].set_ylabel("heat flux [W/m²]")
p.axes[0].set_title("Heat flux at the block-media interface")
p.figure.savefig("heat_flux.png", dpi=300)
```
