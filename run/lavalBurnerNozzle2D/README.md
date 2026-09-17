# Burner nozzle 2D

The goal of this sample case is to simulate premixed combustion across a simple [de Laval](https://en.wikipedia.org/wiki/De_Laval_nozzle) nozzle using detailed kinetics. The model is conceived in axisymmetric 2D coordinates.

## Geometry

1. Use [laval.jl](laval.jl) to compute the required point set.

Alternative 1:
    2. Generate geometry and 2D mesh in GMSH, export as UNV.
    3. Import in Salome and rotate/extrude mesh for final formatting.

Alternative 2:
    2. Generate geometry in FreeCAD, export as STEP format.
    3. Import in Salome, rename surfaces, proceed with meshing.