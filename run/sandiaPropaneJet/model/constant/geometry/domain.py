# -*-  coding: utf -*-

import gmsh

from pathlib import Path

#region: dimensions
d_in   = 0.00526
d_out  = 0.00900
l_pipe = 0.05

s_duct = 0.30
l_duct = 2.0
#endregion: dimensions

#region: gmsh initialization
path = Path(__file__)

gmsh.initialize()
gmsh.model.add(path.stem)

# Get practical aliases
opt = gmsh.option
mod = gmsh.model
occ = gmsh.model.occ
msh = gmsh.model.mesh

# Note: using 3 during concept design is helpful!
opt.set_number("General.Axes", 3)
opt.set_number("Geometry.Points", 0)
opt.set_number("Geometry.Lines", 1)
opt.set_number("Geometry.Surfaces", 1)

opt.set_number("Mesh.StlOneSolidPerSurface", 2)
#endregion

#region: create pipe
x, y, z = (0.0, 0.0, -l_pipe)

# Create inner/outer diameter circles:
circ_inner = occ.add_circle(x, y, z, 0.5 * d_in)
circ_outer = occ.add_circle(x, y, z, 0.5 * d_out)

# Create curve loops from the circles:
loop_outer = occ.add_curve_loop([circ_outer])
loop_inner = occ.add_curve_loop([circ_inner])

# Create a plane surface with a hole (the profile):
pipe_profile = occ.add_plane_surface([loop_outer, loop_inner])
pipe_inlet   = occ.add_plane_surface([loop_inner])

# Create/cnvert the line into a Wire (required by add_pipe):
p_start   = occ.add_point(x, y, z)
p_end     = occ.add_point(x, y, z + l_pipe)
pipe_wire = occ.add_line(p_start, p_end)
pipe_wire = occ.add_wire([pipe_wire])

# Create the pipe volume using the profile and the wire:
pipe_volume = occ.add_pipe([(2, pipe_profile)], pipe_wire)
occ.synchronize()

# Retrive the tags of generated pipe walls:
pipe_walls = [t for d, t in mod.get_boundary(pipe_volume) if d == 2]
pipe_outer, pipe_inner, pipe_end = pipe_walls[1:]
#endregion

#region: create base geometry
base_duct = occ.add_rectangle(
    x  = -s_duct / 2,
    y  = -s_duct / 2,
    z  = -l_pipe,
    dx = s_duct,
    dy = s_duct,
)

outlet = occ.add_rectangle(
    x  = -s_duct / 2,
    y  = -s_duct / 2,
    z  = l_duct,
    dx = s_duct,
    dy = s_duct,
)

# Remove the pipe elements from base_duct:
tools = [(2, pipe_profile), (2, pipe_inlet)]
coflow_inlet = occ.fragment([(2, base_duct)], tools)[0][0][1]
occ.synchronize()

# Extrude duct walls around the system:
coflow_bounds = mod.get_boundary([(2, coflow_inlet)])[:4]
duct_walls = occ.extrude(coflow_bounds, 0.0, 0.0, l_pipe + l_duct)
duct_walls = [t for d, t in duct_walls if d == 2]
occ.synchronize()
#endregion

#region: add physical surfaces
mod.add_physical_group(
    dim  = 2,
    tags = [pipe_inlet],
    name = "inletFuel"
)
mod.add_physical_group(
    dim  = 2,
    tags = [coflow_inlet],
    name = "inletCoflow"
)

mod.add_physical_group(
    dim  = 2,
    tags = [pipe_outer],
    name = "pipeOuterWall"
)
mod.add_physical_group(
    dim  = 2,
    tags = [pipe_inner],
    name = "pipeInnerWall"
)
mod.add_physical_group(
    dim  = 2,
    tags = [pipe_end],
    name = "pipeEndWall"
)

mod.add_physical_group(
    dim  = 2,
    tags = duct_walls,
    name = "ductWalls"
)
mod.add_physical_group(
    dim  = 2,
    tags = [outlet],
    name = "outlet"
)
#endregion

#region: generate and dump
gmsh.write(fileName=str(path.with_suffix(".stl")))

gmsh.fltk.run()
gmsh.finalize()
#endregion