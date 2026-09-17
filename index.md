# OpenFOAM Extensions

Curated OpenFOAM 13 directory with additional workflow tools.

## Compilation methods

It is possible to build the whole extensions library or isolated features, depending on your needs. The following logic is applied throughout this sources directory:

- If a directory provides an `Allwmake` file, it allows to compile the whole tree below it; this is the case for the root directory, for instance.

- If a directory of `src/` provides a subdirectory called `Make/`, then it supports the individual library build which is done by running `wmake libso` from that directory.

All builds are written to `$FOAM_USER_LIBBIN`.

## Running tests

Some libraries may provide test programs under a `test/` sub-directory, and all of them are structured in the same way. The following instructions indicate how to build, run, and clean those directories.

To build and run the test suite:

```bash
cd test
./Allrun
```

To clean test build artifacts:

```bash
cd test
./Allclean
```
