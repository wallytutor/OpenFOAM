# Lagrangian

## Parcels

### Extended Thermo Parcel

> *Provides a variant of `thermoCloud` with varying specific heat.*

**Important:** this extension currently support a single solid phase which needs to be a `tabulatedSolid` as provided in the thermophysical properties below. This will need to be generalized at some point.

Extending a parcel requires several templates to be adapted. The following building blocks are present in the implementation of `extendedThermoCloud`:

- Parcel template class in found in [ExtendedThermoParcel](lagrangian/parcel/parcels/Templates/ExtendedThermoParcel)

- The derived parcel instantiation [extendedThermoParcel](lagrangian/parcel/parcels/derived/extendedThermoParcel)

- The derived cloud instantiation [extendedThermoCloud](lagrangian/parcel/clouds/derived/extendedThermoCloud)

Or in a tree view:

```bash
parcel
├── clouds
│   └── derived
│       └── extendedThermoCloud
│           ├── extendedThermoCloud.C
│           └── extendedThermoCloud.H
└── parcels
    ├── Templates
    │   └── ExtendedThermoParcel
    │       ├── ExtendedThermoParcel.C
    │       ├── ExtendedThermoParcel.H
    │       ├── ExtendedThermoParcelI.H
    │       ├── ExtendedThermoParcelIO.C
    │       ├── ExtendedThermoParcelName.C
    │       └── ExtendedThermoParcelTrackingDataI.H
    └── derived
        └── extendedThermoParcel
            ├── extendedThermoParcel.H
            └── makeExtendedThermoParcelSubmodels.C
```

For using it, append the following to `controlDict` (it depends on extended thermophysical properties so that non-constant properties are supported):

```C
libs
(
    "libextendedThermophysicalProperties.so"
    "libextendedLagrangianParcel.so"
);
```