# -*- coding: utf-8 -*-
from pathlib import Path
from heat_recovery.cowper_like  import get_model_dump_mesh

RENDER: bool = False
""" Whether to render the Gmsh GUI during mesh generation. """


def generate_domain():
    get_model_dump_mesh(
        config = "dimensioning.yaml",
        mesh   = "mesh.msh",
        render = RENDER
    )


def post_process():
    pass


if __name__ == "__main__":
    if not Path("mesh.msh").exists():
        generate_domain()
    else:
        post_process()
