# OpenFOAM extensions

Run `Allwmake` from this directory to create the full library. For creating parts of it only, navigate to the required directory and run `Allmake` if present. Otherwise, if a `Make` directory is available, run `wmake libso`. These are the three possible library levels to be built.

## Lagrangian

### Extended Thermo Parcel

> *Provides a variant of `thermoCloud` with varying specific heat.*

Extending a parcel requires several templates to be adapted. The following building blocks are present in the implementation of `extendedThermoCloud`:

- Parcel template class in found in [ExtendedThermoParcel](lagrangian/parcel/parcels/Templates/ExtendedThermoParcel)

- The derived parcel instantiation [extendedThermoParcel](lagrangian/parcel/parcels/derived/extendedThermoParcel)

- The derived cloud instantiation [extendedThermoCloud](lagrangian/parcel/clouds/derived/extendedThermoCloud)

For using it, append the following to `controlDict` (it depends on `libextendedThermophysicalProperties` so that non-constant properties are supported):

```C
libs
(
    "libextendedThermophysicalProperties.so"
    "libextendedLagrangianParcel.so"
);
```

## Thermophysical Models

### Polynomial Solid Properties

### Tabulated Solid Properties

## Developping

- Whenever possible/applicable, copy the original `files` list from the library you are modifying to use as a starting point. Remove or comment out files that are not implemented in this repository. Modify the base path of `LIB` to use `LIB = $(FOAM_USER_LIBBIN)/<your-library-name>` so that the path is writable.

- The same apply to `options`, which requires a few extra guidelines concerning its declared variables:

    - `EXE_INC`: keep the original list untouched, add the local include directories in the end. Sometimes the library will fail to compile because the original list was lacking something (in a library as OpenFOAM it probably worked because the missing files were sourced elsewhere during compilation and the authors never had to add then), append these in the end so that we can distinguish them from the ones already in the list.

    - `LIB_LIBS`: it the project compiles but on run-time there is a missing linking library, but proceed as for `EXE_INC` by listing the library in the end. It might be tricky to identify where the functionality comes from, maybe start by `ls $FOAM_LIBBIN` and look for possibly related files.
