# Custom Temperature-Dependent Solid Properties Library for OpenFOAM v13

This library extends OpenFOAM's `solidProperties` base class using `Function1<scalar>` to support tabulated temperature-dependent density ($\rho(T)$), specific heat capacity ($C_p(T)$), and themal conductivity ($\kappa$).

---

## Directory Structure

```
tabulatedSolidProperties/
├── tabulatedSolid.H
├── tabulatedSolid.C
├── tabulatedSolidI.H
├── Make/
│   ├── files
│   └── options
└── test/
    ├── testTabulatedSolid.C
    ├── Allrun
    ├── Allclean
    └── Make/
        ├── files
        └── options
```

---

## Compilation Instructions

1. Navigate to the `tabulatedSolidProperties` directory:

   ```bash
   cd tabulatedSolidProperties
   ```

2. Compile the dynamic library using `wmake`:

   ```bash
   wmake libso
   ```

   This builds `libtabulatedSolidProperties.so` inside `$FOAM_USER_LIBBIN`.

---

## Usage in Case Files

### `system/controlDict`

Load the compiled library at startup:

```openfoam
libs
(
    "libtabulatedSolidProperties.so"
);
```

### `constant/physicalProperties`

Specify `tabulatedSolid` under `solids {}`:

```openfoam
solids
{
    tabulatedSolid
    {
        molWeight   60.084;
        Hf          0.0;
        emissivity  0.9;

        // Tabulated density rho(T) [kg/m^3]
        rho
        table
        (
            (300  2500)
            (600  2450)
            (1200 2380)
        );

        // Tabulated specific heat Cp(T) [J/(kg K)]
        Cp
        table
        (
            (300  800)
            (500  950)
            (800  1100)
            (1200 1250)
        );

        // Tabulated thermal conductivity [W/(m K)]
        kappa
        table
        (
            (300  1.5)
            (500  1.5)
            (800  1.5)
            (1200 1.5)
        );

    }
}
```

---

## Running Tests

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
