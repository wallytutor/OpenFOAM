# Custom Temperature-Dependent Solid Properties Library for OpenFOAM v13

This library extends OpenFOAM's `solidProperties` base class to support temperature-dependent specific heat capacity ($C_p(T)$) and density ($\rho(T)$) via 8th-order polynomials (`Polynomial<8>`).

---

## Directory Structure

```
polynomialSolidProperties/
├── polynomialSolid.H
├── polynomialSolid.C
├── polynomialSolidI.H
└── Make/
    ├── files
    └── options
```

---

## Compilation Instructions

1. Navigate to the `polynomialSolidProperties` directory:

   ```bash
   cd polynomialSolidProperties
   ```

2. Compile the dynamic library using `wmake`:

   ```bash
   wmake libso
   ```

   This builds `libpolynomialSolidProperties.so` inside `$FOAM_USER_LIBBIN`.

---

## Usage in Case Files

### `system/controlDict`

Load the compiled library at startup:

```openfoam
libs
(
    "libpolynomialSolidProperties.so"
);
```

### `constant/physicalProperties`

Specify `polynomialSolid` under `solids {}`:

```openfoam
solids
{
    polynomialSolid
    {
        molWeight   60.084;
        Hf          0.0;
        emissivity  0.9;

        // Polynomial coefficients for rho(T) = r0 + r1*T + r2*T^2 + ...
        rhoCoeffs    (2500 -0.1 0 0 0 0 0 0);

        // Polynomial coefficients for Cp(T) = c0 + c1*T + c2*T^2 + ...
        CpCoeffs     (800 0.5 -1e-4 0 0 0 0 0);

        // Polynomial coefficients for kappa(T) = c0 + c1*T + c2*T^2 + ...
        kappaCoeffs  (1.5 0 0 0 0 0 0 0);
    }
}
```
