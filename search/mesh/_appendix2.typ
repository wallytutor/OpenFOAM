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

== Extracting surface features

=== surfaceFeatures

=== ParaView and PyVista
