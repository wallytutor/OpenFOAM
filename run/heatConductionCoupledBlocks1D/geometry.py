# -*- coding: utf-8 -*-
from pathlib import Path
from majordome_simulation.meshing import GmshOCCModel
from majordome_simulation.meshing import GeometricProgression
from screeninfo import get_monitors
import numpy as np

RENDER = False

MONITOR = get_monitors()[0]

DEFAULT_OPTIONS = {
    "General.GraphicsWidth": int(MONITOR.width*0.8),
    "General.GraphicsHeight": int(MONITOR.height*0.8),
    "General.Verbosity": 1,
    "Geometry.Points": True,
    "Geometry.Lines": True,
    "Geometry.Surfaces": True,
    "Mesh.ColorCarousel": 2,
    "Mesh.SurfaceFaces": True,
    "Mesh.SaveAll": False,
    "Mesh.SaveGroupsOfNodes": True,
    "Mesh.Algorithm": 6,
    "Mesh.ElementOrder": 1,
    "Mesh.MshFileVersion": 2.2
}


def generate_domain(
        layer: float,
        saveas: str | None = None,
        overwrite: bool = False,
        options: dict = DEFAULT_OPTIONS
    ) -> None:
    """ Generates the mesh for the heat conduction 1D problem. """
    grid_fine: float = 0.001
    grid_coarse: float = 0.05

    side  = 0.1
    depth = 2.0

    thickness1 = 0.0254 * layer
    thickness2 = depth - thickness1

    n, q = GeometricProgression.fit(thickness2, grid_fine, grid_coarse)
    gp = GeometricProgression(n, grid_fine, q=q)

    layers1 = int(thickness1 / grid_fine)

    origin = (0.0, 0.0, 0.0)

    with GmshOCCModel(render=RENDER, **options) as model:
        #region: base
        tag_base = model.add_rectangle(*origin, side, side)
        model.synchronize()

        bounds = model.get_boundary([(2, tag_base)])
        bounds = [tag for dim, tag in bounds if dim == 1]

        for tag in bounds:
            model.set_transfinite_curve(tag, 2)

        model.set_recombine(2, tag_base)
        model.synchronize()
        #endregion: base

        #region: extrude1
        [tag_inter, vol1, *sides1] = model.extrude(
            dimTags     = [(2, tag_base)],
            dx          = 0.0,
            dy          = 0.0,
            dz          = thickness1,
            numElements = [layers1],
            recombine   = True
        )
        model.synchronize()

        tag_inter = tag_inter[1]
        vol1      = vol1[1]
        sides1    = [side[1] for side in sides1]
        #endregion: extrude1

        #region: extrude2
        # XXX this freezes, so I am brute-forcing the extrusion:
        # [tag_infty, vol2, *sides2] = model.extrude(
        #     dimTags     = [(2, tag_inter)],
        #     dx          = 0.0,
        #     dy          = 0.0,
        #     dz          = thickness2,
        #     numElements = [layers2],
        #     recombine   = True
        # )
        # model.synchronize()

        previous   = vol1
        to_extrude = (2, tag_inter)
        vol2   = []
        sides2 = []

        for dh in gp.heights:
            [tag_inter, vol, *sides] = model.extrude(
                dimTags     = [to_extrude],
                dx          = 0.0,
                dy          = 0.0,
                dz          = dh,
                numElements = [1],
                recombine   = True
            )
            model.synchronize()

            model.fragment([(3, previous)], [vol])
            model.synchronize()
            previous = vol[1]

            to_extrude = tag_inter
            vol2.append(vol[1])
            sides2.extend([side[1] for side in sides])
        #endregion: extrude2

        #region: physical groups

        surface     = [tag_base]
        sides_block = sides1
        sides_media = sides2
        infinite    = [to_extrude[1]]

        block = [vol1]
        media = vol2
        model.add_physical_groups(
            surfaces = [
                {"tags": surface,     "name": "surface",     "tag_id": 1},
                {"tags": sides_block, "name": "sides_block", "tag_id": 2},
                {"tags": sides_media, "name": "sides_media", "tag_id": 3},
                {"tags": infinite,    "name": "infinite",    "tag_id": 4},
            ],
            volumes = [
                {"tags": block, "name": "block", "tag_id": 10},
                {"tags": media, "name": "media", "tag_id": 20},
            ],
        )
        model.synchronize()
        #endregion: physical groups

        model.generate_mesh(dim=2)
        model.generate_mesh(dim=3)
        model.synchronize()

        if saveas and (overwrite or not Path(saveas).exists()):
            model.dump(saveas)


if __name__ == "__main__":
    here = Path(__file__).parent

    for layer in [4, 5, 6]:
        micron = f"{layer * 25.4:05.3f}".replace(".", "")
        saveas = str(here / f"mesh_layers_{micron}um.msh")
        generate_domain(layer, saveas=saveas, overwrite=False)
