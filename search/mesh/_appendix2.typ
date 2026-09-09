= Preprocessing geometry <appendix-2>

== Scaling the surfaces

When conceiving a geometry, one must always take care of the units used for dimensioning. Sometimes when working with engineering drawings it may be useful to keep the millimeters or inches used for design, for instance. After exporting the STL surfaces to be used with `snappyHexMesh`, it may be practical to scale dimensions into SI units (meters), as by default most of `OpenFOAM` environment and data is setup to work with them. Several transformations can be applied with `surfaceTransformPoints`, which transforms a surface geometry by translation, rotation and/or scaling. A sample of scaling from millimeters to meters is provided in @lst-surfaceTransformPoints-scaling. For more details regarding the tool and other types of transformations, consider running `surfaceTransformPoints -help` in your terminal.


#figure(
  ```bash
  surfaceTransformPoints           \
      -scale '(0.001 0.001 0.001)' \
      <input.stl>                  \
      <output.stl>
  ```,
  caption: [Use of `surfaceTransformPoints` for scaling an STL surface.],
) <lst-surfaceTransformPoints-scaling>

== Surface orientation <appendix-2-sec-surfaceOrient>

== Extracting surface features <appendix-2-sec-surfaceFeatures>

== ParaView <appendix-2-sec-paraview>

ParaView provides specialized filters and operations for geometric preparation, boundary manipulation, and mesh quality inspection. Its interface allow users to manually highlight individual elements, groups of faces, or polylines directly on the 3D model. It provides flexible visual control to identify localized geometric sections, such as a specific pipe ring, that require independent refinement levels.

When working interativelly in a geometry conception, it may be interesting to save the visualization state. By later using the load state functionality, it allows to restore a previously saved visualization environment, which in the context of geometry conception and meshing means applying automatically an workflow of filters (which can be cumbersome). It is also useful in post-processing, for recovering precise layouts and camera views, but that is beyond our scope here.

- *Normal Glyph:* This visualization filter generates vector arrows across a geometry to display the spatial orientation of its surface normals. It operates as a diagnostic check to verify that all normals consistently point outward, which is critical for accurate intersection calculations and robust mesh generation. This should always be your first step if using `snappyHexMesh`.

- *Crinkle Slice:* This filter cuts through a volumetric mesh while preserving the boundaries of the intersected cells, yielding a jagged rather than flat cross-section. It is primarily utilized to inspect internal mesh structures, identify abrupt transition rates between small and large cells, and visually evaluate regions prone to numerical diffusion. Because of its underlining octree structure, trying to interpret a slice of a mesh generated with `snappyHexMesh` is misleading and should not be done without recurring to this option.

// - *Feature Edges:* This filter identifies and extracts sharp geometric boundaries from a surface based on a user-defined feature angle. *Use Case:* It isolates specific geometric edges, which are subsequently saved as OBJ files and imported into meshing utilities to enforce targeted local mesh refinements.

// - *Extract Selection / Extract Surfaces:* These filters separate the discrete regions highlighted by the selection tools into independent datasets. *Use Case:* Once separated, these isolated subsets can be exported as individual files, enabling the assignment of unique meshing parameters to different portions of the computational domain.

// - *Generate Surface Normal:* This filter computes surface normals for a given geometry and includes a specific parameter to invert their orientation. *Use Case:* It is applied to repair geometric inconsistencies by flipping inward-pointing normals to point outward, ensuring a correctly oriented, watertight surface prior to meshing.

// - *Threshold:* This operation filters a dataset to isolate specific subsets based on discrete ID values. *Use Case:* It is used to separate distinct boundary patches—such as inlets, outlets, and specific pipe walls—from an assembled multi-surface geometry for individual isolation and manual renaming.

== PyVista
