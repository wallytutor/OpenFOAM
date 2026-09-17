# -*- encoding: utf-8 -*-
""" Provides geometry and meshing of wedge nozzle.

This is a 2D simplification (round instead of square) domain for fast
study of the Sandia Propane Jet experiment. It is provided to be used
during model setup and variant cases before moving into full 3D with
the actual system geometry (square duct).

Important: this is conceived with Y as the vertical axis; the full 3D
case is conceived with Z vertically, so consider changing `g` vector.
"""

import math
import gmsh

from pathlib import Path

#region: 0. Configuration
# Small angle for wedge (<5°):
wedge_angle: float = math.radians(2.0)

# Inner diameter of pipe:
d_in: float = 0.00526

# Outer diameter of pipe:
d_out: float = 0.00900

# Arbitrary pipe length:
l_pipe: float = 0.05

# Width of duct (radius here):
s_duct: float = 0.30

# Lenght of duct:
l_duct: float = 2.0

# Refinement length at outlet:
l_exp: float = 3 * d_in
#endregion

#region: 1. Initialize the model
# Path of the current file:
path = Path(__file__)

# Actual gmsh initialization
gmsh.initialize()
gmsh.model.add(path.stem)

# Get practical aliases
opt = gmsh.option
mod = gmsh.model
occ = gmsh.model.occ
msh = gmsh.model.mesh

# Note: using 3 during concept design is helpful!
opt.set_number("General.Axes", 2)
opt.set_number("Geometry.Points", 0)
opt.set_number("Geometry.Lines", 1)
opt.set_number("Geometry.Surfaces", 1)
opt.set_number("Mesh.Algorithm", 6)
opt.set_number("Mesh.SurfaceFaces", 1)
opt.set_number("Mesh.ColorCarousel", 2)
opt.set_number("Mesh.ElementOrder", 1)      # XXX: important!
opt.set_number("Mesh.MshFileVersion", 2.2)  # XXX: important!
#endregion

#region: 2. Create base geometry
dx1 = d_in / 2.0
dx2 = (d_out - d_in) / 2.0
dx3 = (s_duct - d_out) / 2

dy1 = l_pipe
dy2 = l_exp
dy3 = l_duct - l_exp

# (1) Inlet pipe
occ.add_rectangle(
    x   = 0.0,
    y   = -dy1,
    z   = 0.0,
    dx  = dx1,
    dy  = dy1,
    tag = 1,
)

# (2) Inlet co-flow
occ.add_rectangle(
    x   = dx1 + dx2,
    y   = -dy1,
    z   = 0.0,
    dx  = dx3,
    dy  = dy1,
    tag = 2,
)

# (3) Above inlet pipe
occ.add_rectangle(
    x   = 0.0,
    y   = 0.0,
    z   = 0.0,
    dx  = dx1,
    dy  = dy2,
    tag = 3,
)

# (4) Above pipe wall
occ.add_rectangle(
    x   = dx1,
    y   = 0.0,
    z   = 0.0,
    dx  = dx2,
    dy  = dy2,
    tag = 4,
)

# (5) Above inlet co-flow
occ.add_rectangle(
    x   = dx1 + dx2,
    y   = 0.0,
    z   = 0.0,
    dx  = dx3,
    dy  = dy2,
    tag = 5,
)

# (6) Above inlet pipe
occ.add_rectangle(
    x   = 0.0,
    y   = dy2,
    z   = 0.0,
    dx  = dx1,
    dy  = dy3,
    tag = 6,
)

# (7) Above pipe wall
occ.add_rectangle(
    x   = dx1,
    y   = dy2,
    z   = 0.0,
    dx  = dx2,
    dy  = dy3,
    tag = 7,
)

# (8) Above inlet co-flow
occ.add_rectangle(
    x   = dx1 + dx2,
    y   = dy2,
    z   = 0.0,
    dx  = dx3,
    dy  = dy3,
    tag = 8,
)
occ.synchronize()

# Get list of all entities:
all_surf = [(2, i) for i in range(1, 8+1)]

# XXX if rotation is done before fragmentation we end up with the same
# line indices of the original geometry, what is more logical for what
# comes next (easier to remember numbers).

# Rotate all planes for wedge concept:
occ.rotate(
    dimTags = all_surf,
    x       = 0,
    y       = 0,
    z       = 0,
    ax      = 0,
    ay      = 1,
    az      = 0,
    angle   = -wedge_angle / 2.0
)
occ.synchronize()

# Ensure shared boundaries are merged:
occ.fragment(all_surf, [])
occ.synchronize()
#endregion

#region: 3. Discretization before 3-D
nx1 = 6
rx1 = 0.9

nx2 = 4
rx2 = 0.25

nx3 = 40
rx3 = 0.1

ny1 = 20
ry1 = 0.9

ny2 = int(1000 * dy2)
ry2 = 1.0

ny3 = int(500 * dy3)
ry3 = 1.0

msh.set_transfinite_curve(1,  nx1, "Progression", rx1)
msh.set_transfinite_curve(3,  nx1, "Progression", 1/rx1)
msh.set_transfinite_curve(10, nx1, "Progression", 1.0)
msh.set_transfinite_curve(18, nx1, "Progression", 1.0)

msh.set_transfinite_curve(12, nx2, "Bump", rx2)
msh.set_transfinite_curve(14, nx2, "Progression", 1.0)
msh.set_transfinite_curve(21, nx2, "Progression", 1.0)

msh.set_transfinite_curve(5,  nx3, "Bump", rx3)
msh.set_transfinite_curve(7,  nx3, "Bump", rx3)
msh.set_transfinite_curve(16, nx3, "Bump", rx3)
msh.set_transfinite_curve(23, nx3, "Progression", 1.0)

msh.set_transfinite_curve(2,  ny1, "Progression", ry1)
msh.set_transfinite_curve(4,  ny1, "Progression", 1/ry1)
msh.set_transfinite_curve(6,  ny1, "Progression", ry1)
msh.set_transfinite_curve(8,  ny1, "Progression", 1/ry1)

msh.set_transfinite_curve(9,  ny2, "Progression", ry2)
msh.set_transfinite_curve(11, ny2, "Progression", ry2)
msh.set_transfinite_curve(13, ny2, "Progression", ry2)
msh.set_transfinite_curve(15, ny2, "Progression", ry2)

msh.set_transfinite_curve(17, ny3, "Progression", ry3)
msh.set_transfinite_curve(19, ny3, "Progression", ry3)
msh.set_transfinite_curve(20, ny3, "Progression", ry3)
msh.set_transfinite_curve(22, ny3, "Progression", ry3)

for dim, tag in all_surf:
    msh.set_transfinite_surface(tag=tag)
    msh.set_recombine(dim=dim, tag=tag)
#endregion

#region: 4. Revolution to get 3-D
occ.revolve(
    dimTags     = all_surf,
    x           = 0,
    y           = 0,
    z           = 0,
    ax          = 0,
    ay          = 1,
    az          = 0,
    angle       = wedge_angle,
    numElements = [1],
    recombine   = True
)
occ.synchronize()
#endregion

#region: 5. Add physical groups
config_groups = [
    {
        "name": "inletFuel",
        "tags": [9],
        "id": 10
    },
    {
        "name": "inletCoflow",
        "tags": [13],
        "id": 11
    },
    {
        "name": "outlet",
        "tags": [29, 32, 35],
        "id": 12
    },
    {
        "name": "pipeOuterWall",
        "tags": [16],
        "id": 20
    },
    {
        "name": "pipeInnerWall",
        "tags": [10],
        "id": 21
    },
    {
        "name": "pipeEndWall",
        "tags": [21],
        "id": 22
    },
    {
        "name": "ductWalls",
        "tags": [14, 25, 34],
        "id": 23
    },
    {
        "name": "front",
        "tags": [1, 2, 3, 4, 5, 6, 7, 8],
        "id": 30,
    },
    {
        "name": "back",
        "tags": [12, 17, 20, 24, 27, 30, 33, 36],
        "id": 31,
    },
]

for group in config_groups:
    mod.add_physical_group(
        dim  = 2,
        tags = group["tags"],
        tag  = group["id"],
        name = group["name"]
    )

mod.add_physical_group(
    dim  = 3,
    tags = list(range(1, 8+1)),
    tag  = 100,
    name = "fluid"
)

msh.generate(dim=3)
gmsh.write(str(path.with_suffix(".msh")))
#endregion

#region: 6. Display and finalize
gmsh.fltk.run()
gmsh.finalize()
#endregion