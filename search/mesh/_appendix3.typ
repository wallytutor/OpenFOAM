= Troubleshooting <appendix-3>

== Dealing with `wrongFraces`

*Quick action:* To identify the problematic region, start by converting the face set to VTK with `foamToVTK -faceSet wrongFaces` so that it can be inspected in ParaView or PyVista.

A `faceSet` named `wrongFaces` generated after running `snappyHexMesh` indicates topological and orientation inconsistencies where face normal vectors disagree with neighbor-cell definitions or boundary connectivity. It may be caused by inconsistent of inverted STL normals (see @appendix-2-sec-paraview), self-intersecting or non-watertight geometry, over-aggressive snapping, improper baffle/face zone configuration, among others.

Generally, the easier to detect issue is the STL normals. If adjacent triangles on the input geometry have opposite normal directions or if internal feature edges create conflicting surface orientations, snappyHexMesh generates misoriented boundary faces. Other problems with geometry can be related to gaps, holes, or overlapping surfaces confuse the inside/outside determination and ray-tracing logic. Always thoroughly inspect the geometry before starting the meshing process.

If the problem is related to STL files, you can re-generate the geometry; for instance, if exporting from `gmsh`, inverting the direction vector used for extrusion of surfaces can in some cases fix the exported STL files. Otherwise, you can manipulate STL files with built-in tools as `surfaceOrient` as discussed in @appendix-2-sec-surfaceOrient.

The other group of root causes is related to the configuration. Extreme mesh deformation during the snap or layer addition stages can invert or collapse face pyramids, leading to invalid face-to-cell topology. Relaxing their inputs can be used to validate this hypothesis. To do so, start by increasing `nSolveIter` and `nRelaxIter`.

Finally, meshing interior sheets or non-manifold surfaces without correctly configuring `faceType` baffle; or internal face zones can be the origin; please refer to XXX for dealing with baffle setups
