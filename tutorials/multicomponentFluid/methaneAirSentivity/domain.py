# -*- encoding: utf-8 -*-

import math
import gmsh

from pathlib import Path

#region: 0. Configuration
path = Path(__file__)
wedge_angle = math.radians(2.0)

len_inlet_ax = 0.05
len_inlet_an = 0.03
len_disperse = 1.00
len_flue_out = 0.50

dia_inlet_ax = 0.025
thk_inlet_wl = 0.003
thk_inlet_an = 0.003
thk_side_box = 0.100

extended = False
#endregion

#region: 1. Initialize the model
gmsh.initialize()
gmsh.model.add(path.stem)

# Get practical aliases
opt = gmsh.option
mod = gmsh.model
occ = gmsh.model.occ
msh = gmsh.model.mesh
#endregion

#region: 2. Configure gmsh options
# Note: using 3 during concept design is helpful!
opt.set_number("General.Axes", 0)

opt.set_number("Geometry.Points", 0)
opt.set_number("Geometry.Lines", 1)
opt.set_number("Geometry.Surfaces", 1)

opt.set_number("Mesh.Algorithm", 6)
opt.set_number("Mesh.SurfaceFaces", 1)

# Show each named surface with a different color
opt.set_number("Mesh.ColorCarousel", 2)

# XXX: mandatory for OpenFOAM meshing!
opt.set_number("Mesh.ElementOrder", 1)
opt.set_number("Mesh.MshFileVersion", 2.2)

# XXX: mesh-specific settings
opt.set_number("Mesh.MeshSizeMax", 0.01)
#endregion

#region: 3. Create base geometry
# - Axial inlet (negative x-direction)
occ.add_rectangle(
    x   = 0.0,
    y   = 0.0,
    z   = 0.0,
    dx  = -len_inlet_ax,
    dy  = dia_inlet_ax / 2.0,
    tag = 1
)

# - Dispersion length for axial inlet
occ.add_rectangle(
    x   = 0.0,
    y   = 0.0,
    z   = 0.0,
    dx  = len_disperse,
    dy  = dia_inlet_ax / 2.0,
    tag = 2
)

# - Dispersion length over wall
occ.add_rectangle(
    x   = 0.0,
    y   = dia_inlet_ax / 2.0,
    z   = 0.0,
    dx  = len_disperse,
    dy  = thk_inlet_wl,
    tag = 3
)

# - Annular inlet (negative x-direction)
occ.add_rectangle(
    x   = 0.0,
    y   = dia_inlet_ax / 2.0 + thk_inlet_wl,
    z   = 0.0,
    dx  = -len_inlet_an,
    dy  = thk_inlet_an,
    tag = 4
)

# - Dispersion length for annular inlet
occ.add_rectangle(
    x   = 0.0,
    y   = dia_inlet_ax / 2.0 + thk_inlet_wl,
    z   = 0.0,
    dx  = len_disperse,
    dy  = thk_inlet_an,
    tag = 5
)

# - Side expansion width
occ.add_rectangle(
    x   = 0.0,
    y   = dia_inlet_ax / 2.0 + thk_inlet_wl + thk_inlet_an,
    z   = 0.0,
    dx  = len_disperse,
    dy  = thk_side_box,
    tag = 6
)

if not extended:
    all_surf = [(2, k) for k in range(1, 6+1)]
else:
    # - Flue outlet (4 regions)
    occ.add_rectangle(
        x   = len_disperse,
        y   = 0.0,
        z   = 0.0,
        dx  = len_flue_out,
        dy  = dia_inlet_ax / 2.0,
        tag = 7
    )
    occ.add_rectangle(
        x   = len_disperse,
        y   = dia_inlet_ax / 2.0,
        z   = 0.0,
        dx  = len_flue_out,
        dy  = thk_inlet_wl,
        tag = 8
    )
    occ.add_rectangle(
        x   = len_disperse,
        y   = dia_inlet_ax / 2.0 + thk_inlet_wl,
        z   = 0.0,
        dx  = len_flue_out,
        dy  = thk_inlet_an,
        tag = 9
    )
    occ.add_rectangle(
        x   = len_disperse,
        y   = dia_inlet_ax / 2.0 + thk_inlet_wl + thk_inlet_an,
        z   = 0.0,
        dx  = len_flue_out,
        dy  = thk_side_box,
        tag = 10
    )

    all_surf = [(2, k) for k in range(1, 10+1)]

occ.rotate(
    dimTags = all_surf,
    x       = 0,
    y       = 0,
    z       = 0,
    ax      = 1,
    ay      = 0,
    az      = 0,
    angle   = -wedge_angle / 2.0
)
occ.synchronize()

# Ensure shared boundaries are merged
occ.fragment(all_surf, [])
occ.synchronize()
#endregion

#region: 4. Discretize before going 3-D
# ---------------------------------------------------------------------
# AXIAL INLET ZONE
# ---------------------------------------------------------------------

# - Axial inlet zone over x-axis
nx1 = 20
qx1 = 1.0
msh.set_transfinite_curve(1, nx1, "Progression", qx1)
msh.set_transfinite_curve(3, nx1, "Progression", qx1)

# - Axial inlet zone over y-axis
ny1 = 10
qy1 = 0.9
msh.set_transfinite_curve(2,  ny1, "Progression", qy1)
msh.set_transfinite_curve(4,  ny1, "Progression", 1/qy1)
msh.set_transfinite_curve(6,  ny1, "Progression", qy1)

# ---------------------------------------------------------------------
# ANNULAR INLET ZONE
# ---------------------------------------------------------------------

# - Annular inlet zone over x-axis
nx2 = 10
qx2 = 1.0
msh.set_transfinite_curve(11, nx2, "Progression", qx2)
msh.set_transfinite_curve(13, nx2, "Progression", qx2)

# - Annular inlet zone over y-axis
ny2 = 7
qy2 = 0.4
msh.set_transfinite_curve(12, ny2, "Bump", qy2)
msh.set_transfinite_curve(14, ny2, "Bump", qy2)
msh.set_transfinite_curve(15, ny2, "Bump", qy2)

# ---------------------------------------------------------------------
# DISPERSION ZONE
# ---------------------------------------------------------------------

# - Dispersion zone over x-axis
nx3 = 200
qx3 = 1.0
msh.set_transfinite_curve(5,  nx3, "Progression", qx3)
msh.set_transfinite_curve(7,  nx3, "Progression", qx3)
msh.set_transfinite_curve(9,  nx3, "Progression", qx3)
msh.set_transfinite_curve(16, nx3, "Progression", qx3)
msh.set_transfinite_curve(18, nx3, "Progression", qx3)

# - Inlet wall thickness zone over y-axis
ny3 = 5
qy3 = 0.9
msh.set_transfinite_curve(8,  ny3, "Progression", qy3)
msh.set_transfinite_curve(10, ny3, "Progression", 1/qy3)

# ---------------------------------------------------------------------
# SIDE EXPANSION ZONE
# ---------------------------------------------------------------------

ny4 = 25
qy4 = 0.85
msh.set_transfinite_curve(17, ny4, "Progression", 1/qy4)
msh.set_transfinite_curve(19, ny4, "Progression", qy4)

# ---------------------------------------------------------------------
# OTHER
# ---------------------------------------------------------------------

if extended:
    # - Flue outlet over x-axis
    nx = 100
    qx = 1.0
    msh.set_transfinite_curve(20, nx, "Progression", qx)
    msh.set_transfinite_curve(22, nx, "Progression", qx)
    msh.set_transfinite_curve(24, nx, "Progression", qx)
    msh.set_transfinite_curve(26, nx, "Progression", qx)
    msh.set_transfinite_curve(28, nx, "Progression", qx)

    # Axial inlet zone
    msh.set_transfinite_curve(21, ny1, "Progression", qy1)

    # Annular inlet zone
    msh.set_transfinite_curve(25, ny2, "Bump", qy2)

    # Inlet wall thickness
    msh.set_transfinite_curve(23, ny3, "Progression", 1/qy3)

    # Side expansion
    msh.set_transfinite_curve(27, ny4, "Progression", 1/qy4)

for _, tag in all_surf:
    msh.set_transfinite_surface(tag=tag)
    msh.set_recombine(dim=2, tag=tag)
#endregion

#region: 5. Key step of revolution to get 3-D
occ.revolve(
    dimTags     = all_surf,
    x           = 0,
    y           = 0,
    z           = 0,
    ax          = 1,
    ay          = 0,
    az          = 0,
    angle       = wedge_angle,
    numElements = [1],
    recombine   = True
)
occ.synchronize()
#endregion

#region: 6. Provide labels to the domain:
s_back = [k for _, k in all_surf]
v_tags = s_back.copy()

if extended:
    s_front    = [14, 17, 21, 26, 29, 33, 36, 39, 42, 45]
    s_inlet_ax = [11]
    s_inlet_an = [23]
    s_outlet   = [34, 37, 40, 43]
    s_walls    = [12, 20, 22, 24, 31, 32, 44]
else:
    s_front    = [10, 13, 17, 22, 25, 29]
    s_inlet_ax = [7]
    s_inlet_an = [19]
    s_outlet   = [11, 14, 23, 26]
    s_walls    = [8, 16, 18, 20, 27, 28]

mod.add_physical_group(
    dim  = 2,
    tags = s_back,
    tag  = 10,
    name = "back"
)
mod.add_physical_group(
    dim  = 2,
    tags = s_front,
    tag  = 11,
    name = "front"
)
mod.add_physical_group(
    dim  = 2,
    tags = s_inlet_ax,
    tag  = 20,
    name = "inletAxial"
)
mod.add_physical_group(
    dim  = 2,
    tags = s_inlet_an,
    tag  = 21,
    name = "inletAnnular"
)
mod.add_physical_group(
    dim  = 2,
    tags = s_outlet,
    tag  = 22,
    name = "outlet"
)
mod.add_physical_group(
    dim  = 2,
    tags = s_walls,
    tag  = 23,
    name = "walls"
)
mod.add_physical_group(
    dim  = 3,
    tags = v_tags,
    tag  = 100,
    name = "fluid"
)
#endregion

#region: 7. Generate, dump, display and finalize
msh.generate(dim=3)
gmsh.write(fileName=str(path.with_suffix(".msh")))

gmsh.fltk.run()
gmsh.finalize()
#endregion
