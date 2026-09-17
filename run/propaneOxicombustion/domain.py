# -*- encoding: utf-8 -*-

import math
import gmsh

from pathlib import Path

#region: 0. Configuration
path = Path(__file__)
wedge_angle = math.radians(2.0)

len_inlet_ax = 0.05
len_inlet_an = 0.03
len_disperse = 0.50
len_flue_out = 0.50

dia_inlet_ax = 0.025
thk_inlet_wl = 0.003
thk_inlet_an = 0.003
thk_side_box = 0.100

f_side_struc = 0.3
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

# 5 is better for mesh expansion
opt.set_number("Mesh.Algorithm", 5)
opt.set_number("Mesh.SurfaceFaces", 1)

# Show each named surface with a different color
opt.set_number("Mesh.ColorCarousel", 2)

# Turn off automatic interpolation/constraints from boundary entities
opt.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
opt.setNumber("Mesh.MeshSizeFromPoints", 0)
opt.setNumber("Mesh.MeshSizeFromCurvature", 0)

# XXX: mandatory for OpenFOAM meshing!
opt.set_number("Mesh.ElementOrder", 1)
opt.set_number("Mesh.MshFileVersion", 2.2)

# XXX: mesh-specific settings
# opt.set_number("Mesh.MeshSizeMax", 0.008)
#endregion

#region: 3. Create base geometry
# Radius of axial inlet:
r_ax = dia_inlet_ax / 2.0

# Inner radius of annular inlet:
r_an = dia_inlet_ax / 2.0 + thk_inlet_wl

# Position and step of structured side:
h_1  = dia_inlet_ax / 2.0 + thk_inlet_wl + thk_inlet_an
dh_1 = thk_side_box * f_side_struc

# Position and step of unstructured side:
h_2  = h_1 + dh_1
dh_2 = thk_side_box * (1 - f_side_struc)

# Total height of slice:
dh_t = h_2 + dh_2

# - Axial inlet (negative x-direction)
occ.add_rectangle(
    x   = 0.0,
    y   = 0.0,
    z   = 0.0,
    dx  = -len_inlet_ax,
    dy  = r_ax,
    tag = 1
)

# - Dispersion length for axial inlet
occ.add_rectangle(
    x   = 0.0,
    y   = 0.0,
    z   = 0.0,
    dx  = len_disperse,
    dy  = r_ax,
    tag = 2
)

# - Dispersion length over wall
occ.add_rectangle(
    x   = 0.0,
    y   = r_ax,
    z   = 0.0,
    dx  = len_disperse,
    dy  = thk_inlet_wl,
    tag = 3
)

# - Annular inlet (negative x-direction)
occ.add_rectangle(
    x   = 0.0,
    y   = r_an,
    z   = 0.0,
    dx  = -len_inlet_an,
    dy  = thk_inlet_an,
    tag = 4
)

# - Dispersion length for annular inlet
occ.add_rectangle(
    x   = 0.0,
    y   = r_an,
    z   = 0.0,
    dx  = len_disperse,
    dy  = thk_inlet_an,
    tag = 5
)

# - Side expansion width (structured)
occ.add_rectangle(
    x   = 0.0,
    y   = h_1,
    z   = 0.0,
    dx  = len_disperse,
    dy  = dh_1,
    tag = 6
)

# - Side expansion width (unstructured)
occ.add_rectangle(
    x   = 0.0,
    y   = h_2,
    z   = 0.0,
    dx  = len_disperse,
    dy  = dh_2,
    tag = 7
)

# - Flue outlet
occ.add_rectangle(
    x   = len_disperse,
    y   = 0.0,
    z   = 0.0,
    dx  = len_flue_out,
    dy  = dh_t,
    tag = 8
)

# Get list of all surfaces:
all_surf = [(2, k) for k in range(1, 8+1)]

# Rotate the whole system by half wedge angle:
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
qy2 = 0.25
msh.set_transfinite_curve(12, ny2, "Bump", qy2)
msh.set_transfinite_curve(14, ny2, "Bump", qy2)
msh.set_transfinite_curve(15, ny2) # XXX no bump

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
qy4 = 0.93
msh.set_transfinite_curve(17, ny4, "Progression", 1.0)
msh.set_transfinite_curve(19, ny4, "Progression", qy4)

ny = 10
qy = 0.9
msh.set_transfinite_curve(22, ny, "Progression", qy)

# ---------------------------------------------------------------------
# EXTENSION
# ---------------------------------------------------------------------

# - Flue outlet over x-axis
nx = 100
qx = 1.0
# msh.set_transfinite_curve(23, nx, "Progression", qx)

# ---------------------------------------------------------------------
# BACKGROUND FIELDS
# ---------------------------------------------------------------------

def threshold_restricted(
        curves_list,
        size_min = 0.003,
        size_max = 0.015,
        dist_min = 0.005,
        dist_max = 0.100,
        sampling = 100,
        sigmoid  = 0,
        restrict = None
    ):
    dist_field = msh.field.add("Distance")
    msh.field.setNumbers(dist_field, "CurvesList", curves_list)
    msh.field.setNumber(dist_field, "Sampling", sampling)

    thresh_field = msh.field.add("Threshold")
    msh.field.setNumber(thresh_field, "InField", dist_field)
    msh.field.setNumber(thresh_field, "SizeMin", size_min)
    msh.field.setNumber(thresh_field, "SizeMax", size_max)
    msh.field.setNumber(thresh_field, "DistMin", dist_min)
    msh.field.setNumber(thresh_field, "DistMax", dist_max)
    msh.field.setNumber(thresh_field, "Sigmoid", sigmoid)

    if restrict is None:
        return thresh_field

    field = msh.field.add("Restrict")
    msh.field.setNumber(field, "InField", thresh_field)
    msh.field.setNumbers(field, "SurfacesList", restrict)
    return field


field1 = threshold_restricted(
    curves_list = [6, 8, 15, 17, 23],
    restrict    = [8]
)
field2 = threshold_restricted(
    curves_list = [18],
    restrict    = [7],
)

min_field = msh.field.add("Min")
msh.field.setNumbers(min_field, "FieldsList", [field1, field2])
msh.field.setAsBackgroundMesh(min_field)

# ---------------------------------------------------------------------
# TRANSFINITE/RECOMBINE
# ---------------------------------------------------------------------

for _, tag in all_surf[:6]:
    msh.set_transfinite_surface(tag=tag)

for _, tag in all_surf[:6]:
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

s_front    = [12, 15, 19, 24, 27, 31, 35, 38]
s_inlet_ax = [9]
s_inlet_an = [21]
s_outlet   = [36]
s_walls    = [10, 18, 20, 22, 30, 33, 34, 37]

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
