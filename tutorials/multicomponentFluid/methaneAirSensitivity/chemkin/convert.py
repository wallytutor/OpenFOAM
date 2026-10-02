# -*- coding: utf-8 -*-

import majordome as mj
import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent


def fit_transport(
        mech: str,
        phase: str,
        T_space: tuple[float, float, int] = (400, 2300, 100)
    ) -> mj.SutherlandFitting:
    """ Fit Sutherland coefficients for a mechanism and phase. """
    sutherland = mj.SutherlandFitting(mech, name=phase)
    sutherland.fit(np.linspace(*T_space))

    transport = HERE / "transport"
    transport.unlink(missing_ok=True)
    sutherland.as_openfoam_dict(transport)

    # Use the following for inspection:
    figs = HERE / "media"
    figs.mkdir(exist_ok=True)

    for s in sutherland.coefs_table["species"]:
        p = sutherland.plot_species(s)
        p.savefig(figs / f"{s}.png")


if __name__ == "__main__":
    fit_transport("BFER_methane.yaml", "CH4_BFER_multi")
