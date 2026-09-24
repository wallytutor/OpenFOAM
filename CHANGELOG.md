# Changelog

## 2026-09-24

- Enhanced `safe_run_ncores` in `etc/rc/coresrc` to accept an arbitrary sequence of fallback core counts. If the primary requested core count exceeds physical cores, fallback values are evaluated in order; if none are feasible or provided, execution falls back to half of the available physical cores. Also changed display of messages by sourcin `utilsrc` from `corerc` and improved re-sourcing time by testing for `WM_PROJECT` in `bashrc`.

## 2026-09-16

- Adapted `extendedThermoCloud` to upcast the solid object to a `tabulatedSolid` as a prototype for temperature dependent properties. Currently this supports a single species and need to be generalized in the future.

- Unified the build of all libraries under a single `Allwmake` at the sources root, following the full library standard. Added compilation instructions at multiple levels to `README.md`

- Added `extendedThermoCloud` support by extending `thermoCloud`; this extension if fact only needed working with `thermoParcel` as the physics was found to be at the parcel level. Added related library sources under `lagrangian/parcel/{cloud,parcels}`.