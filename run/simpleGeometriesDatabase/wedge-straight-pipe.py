# -*- encoding: utf-8 -*-

import math
import gmsh

from pathlib import Path

name = Path(__file__).stem
wedge_angle = math.radians(5.0)

# 1. Initialize the model:
gmsh.initialize()
gmsh.model.add(name)

# 2. Get practical aliases:
opt = gmsh.option
mod = gmsh.model
occ = gmsh.model.occ
msh = gmsh.model.mesh

# 3. Configure gmsh options:
opt.set_number("Geometry.Points", 0)
opt.set_number("Geometry.Lines", 1)
opt.set_number("Geometry.Surfaces", 1)
opt.set_number("Mesh.Algorithm", 6)
opt.set_number("Mesh.SurfaceFaces", 1)
opt.set_number("Mesh.ColorCarousel", 2)
opt.set_number("Mesh.ElementOrder", 1)      # XXX: important!
opt.set_number("Mesh.MshFileVersion", 2.2)  # XXX: important!

# 4. Create base geometry:
occ.add_rectangle(0, 0, 0, 0.5, 0.025, tag=1)
occ.rotate([(2, 1)], 0, 0, 0, 1, 0, 0, -wedge_angle / 2.0)
occ.synchronize()

# 5. Discretize before going 3-D:
msh.set_transfinite_curve(1, 100, "Progression", 1.0)
msh.set_transfinite_curve(3, 100, "Progression", 1.0)
msh.set_transfinite_curve(2, 25, "Progression", 0.9)
msh.set_transfinite_curve(4, 25, "Progression", 1/0.9)
msh.set_transfinite_surface(tag=1)
msh.set_recombine(dim=2, tag=1)

# 6. Key step of revolution to get 3-D:
occ.revolve(
    dimTags     = [(2, 1)],
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

# 7. Provide labels to the domain:
mod.add_physical_group(dim=2, tags=[4], tag=10,  name="inlet")
mod.add_physical_group(dim=2, tags=[2], tag=20,  name="outlet")
mod.add_physical_group(dim=2, tags=[3], tag=30,  name="wall")
mod.add_physical_group(dim=2, tags=[1], tag=40,  name="front")
mod.add_physical_group(dim=2, tags=[5], tag=50,  name="back")
mod.add_physical_group(dim=3, tags=[1], tag=100, name="volume")

# 8. Generate mesh and dump to file:
msh.generate(dim=3)
gmsh.write(f"{name}.msh")

# 9. Display and finalize:
gmsh.fltk.run()
gmsh.finalize()
