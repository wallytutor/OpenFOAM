# -*- encoding: utf-8 -*-

import math
import gmsh

from pathlib import Path

#region: 0. Configuration
name = Path(__file__).stem
wedge_angle = math.radians(5.0)

melt_radius = 0.020
melt_height = 0.040
base_height = 0.015
wall_radius = 0.027 - melt_radius
#endregion: 0. Configuration

#region: 1. Initialize the model:
gmsh.initialize()
gmsh.model.add(name)
#endregion: 1. Initialize the model:

#region: 2. Get practical aliases:
opt = gmsh.option
mod = gmsh.model
occ = gmsh.model.occ
msh = gmsh.model.mesh
#endregion: 2. Get practical aliases:

#region: 3. Configure gmsh options:
opt.set_number("Geometry.Points", 0)
opt.set_number("Geometry.Lines", 1)
opt.set_number("Geometry.Surfaces", 1)
opt.set_number("Mesh.Algorithm", 6)
opt.set_number("Mesh.SurfaceFaces", 1)
opt.set_number("Mesh.ColorCarousel", 2)
opt.set_number("Mesh.ElementOrder", 1)      # XXX: important!
opt.set_number("Mesh.MshFileVersion", 2.2)  # XXX: important!
#endregion: 3. Configure gmsh options:

#region: 4. Create base geometry:
occ.add_rectangle(0,                     0, 0, melt_height, melt_radius, tag=1)
occ.add_rectangle(0,           melt_radius, 0, melt_height, wall_radius, tag=2)
occ.add_rectangle(melt_height,           0, 0, base_height, melt_radius, tag=3)
occ.add_rectangle(melt_height, melt_radius, 0, base_height, wall_radius, tag=4)

all_surf = [(2, 1), (2, 2), (2, 3), (2, 4)]

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
#endregion: 4. Create base geometry:

#region: 5. Ensure shared boundaries are merged:
occ.fragment(all_surf, [])
occ.synchronize()
#endregion: 5. Ensure shared boundaries are merged:

#region: 6. Discretize before going 3-D:
msh.set_transfinite_curve(2, 20, "Progression", 0.9)
msh.set_transfinite_curve(4, 20, "Progression", 1/0.9)
msh.set_transfinite_curve(9, 20, "Progression", 0.9)

msh.set_transfinite_curve(1, 40, "Progression", 1.0)
msh.set_transfinite_curve(3, 40, "Progression", 1.0)
msh.set_transfinite_curve(6, 40, "Progression", 1.0)

msh.set_transfinite_curve(5,  14, "Bump", 0.25)
msh.set_transfinite_curve(7,  14, "Bump", 0.25)
msh.set_transfinite_curve(11, 14, "Bump", 0.25)

msh.set_transfinite_curve(8,  15, "Progression", 1/0.9)
msh.set_transfinite_curve(10, 15, "Progression", 0.9)
msh.set_transfinite_curve(12, 15, "Progression", 0.9)

for tag in range(1, 4+1):
    msh.set_transfinite_surface(tag=tag)
    msh.set_recombine(dim=2, tag=tag)
#endregion: 6. Discretize before going 3-D:

#region: 7. Key step of revolution to get 3-D:
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
#endregion: 7. Key step of revolution to get 3-D:

#region: 8. Provide labels to the domain:
wall_base = [10, 11, 13, 16, 17]

mod.add_physical_group(dim=2, tags=[1],          tag=10, name="fluidFront")
mod.add_physical_group(dim=2, tags=[8],          tag=11, name="fluidBack")
mod.add_physical_group(dim=2, tags=[7],          tag=12, name="surfFluid")

mod.add_physical_group(dim=2, tags=[2, 3, 4],    tag=20, name="solidFront")
mod.add_physical_group(dim=2, tags=[12, 15, 18], tag=21, name="solidBack")
mod.add_physical_group(dim=2, tags=wall_base,    tag=22, name="wallSolid")

mod.add_physical_group(dim=3, tags=[1],          tag=100, name="fluid")
mod.add_physical_group(dim=3, tags=[2, 3, 4],    tag=200, name="solid")
#endregion: 8. Provide labels to the domain:

#region: 9. Generate mesh and dump to file:
msh.generate(dim=3)
gmsh.write(f"{name}.msh")
#endregion: 9. Generate mesh and dump to file:

#region: 10. Display and finalize:
gmsh.fltk.run()
gmsh.finalize()
#endregion: 10. Display and finalize:
