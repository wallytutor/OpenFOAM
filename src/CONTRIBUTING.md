# Contributing

## Directives

- Comments added by this repository to explain changes performed to base libraries are tagged with `// WWW:` so that they can be tracked in the future.

- Whenever possible/applicable, copy the original `files` list from the library you are modifying to use as a starting point. Remove or comment out files that are not implemented in this repository. Modify the base path of `LIB` to use `LIB = $(FOAM_USER_LIBBIN)/<your-library-name>` so that the path is writable.

- The same apply to `options`, which requires a few extra guidelines concerning its declared variables:

    - `EXE_INC`: keep the original list untouched, add the local include directories in the end. Sometimes the library will fail to compile because the original list was lacking something (in a library as OpenFOAM it probably worked because the missing files were sourced elsewhere during compilation and the authors never had to add then), append these in the end so that we can distinguish them from the ones already in the list.

    - `LIB_LIBS`: it the project compiles but on run-time there is a missing linking library, but proceed as for `EXE_INC` by listing the library in the end. It might be tricky to identify where the functionality comes from, maybe start by `ls $FOAM_LIBBIN` and look for possibly related files.

## Wishlist

- [ ] Unification/generalization of sub-classes of `solidProperties` as a new `extendedSolidProperties`; this can then be used as a replacement to `solidProperties` in the `extendedThermoParcel` framework (possible need to generalize the cloud definition).

## Building the docs

A `Makefile` is provided in the root directory: running `make` will generate a local documentation, `make publish` can be used by admins to publish to GitHub Pages, and `make clean` strips the generated outputs and intermediate files from the folder.
