= Generating STL files in Gmsh

== Options

Instead of post-processing a large STL file containing the full geometry, it is easier to write a loop for adding one surface at a time to the STL file. This ensures that each face has its own ID for snappyHexMesh to recognise as a patch. This is illustrated in @lst-gmsh-options-stl-single-patch.

#figure(
```python
gmsh.option.setNumber("Mesh.StlOneSolidPerSurface", 2)

# ...

gmsh.model.add_physical_group(
    dim=2,
    tags=tags,
    name=name
)
gmsh.write(f"{name}.stl")
gmsh.model.remove_physical_groups()
```,
caption: [Enforcing the export of one solid per surface.],
) <lst-gmsh-options-stl-single-patch>

One strategy to enforce finer STL resolution is to increase the accuracy of the CAD geometry triangulation process itself. One may do that by adjusting the following options:

#figure(
```python
gmsh.option.setNumber("Mesh.StlLinearDeflection", 0.0001)
gmsh.option.setNumber("Mesh.StlAngularDeflection", 0.05)
```,
caption: [Refining STL through geometry accuracy setup.],
) <lst-gmsh-options-stl-refinement>