# Thermophysical Models

## Thermophysical Properties

### Polynomial Solid Properties

*Library:* `thermophysicalProperties/solidProperties/polynomialSolidProperties`

Extends OpenFOAM's `solidProperties` base class using `Polynomial<8>` to support temperature-dependent 8th-order polynomials density ($\rho(T)$), specific heat capacity ($C_p(T)$), and themal conductivity ($\kappa$).

It is structured as follows:

```
├── Make
│   ├── files
│   └── options
├── polynomialSolid.C
├── polynomialSolid.H
├── polynomialSolidI.H
└── test
    ├── Allclean
    ├── Allrun
    ├── Make
    │   ├── files
    │   └── options
    └── testPolynomialSolid.C
```

In `system/controlDict` one must load the compiled library at startup:

```openfoam
libs
(
    "libpolynomialSolidProperties.so"
);
```

In `constant/physicalProperties` you can now specify `polynomialSolid` under `solids {}`:

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

### Tabulated Solid Properties

Extends OpenFOAM's `solidProperties` base class using `Function1<scalar>` to support tabulated temperature-dependent density ($\rho(T)$), specific heat capacity ($C_p(T)$), and themal conductivity ($\kappa$).

It is structured as follows:

```
├── Make
│   ├── files
│   └── options
├── tabulatedSolid.C
├── tabulatedSolid.H
├── tabulatedSolidI.H
└── test
    ├── Allclean
    ├── Allrun
    ├── Make
    │   ├── files
    │   └── options
    └── testTabulatedSolid.C
```

In `system/controlDict` one must load the compiled library at startup:

```openfoam
libs
(
    "libtabulatedSolidProperties.so"
);
```

In `constant/physicalProperties` you can now specify `tabulatedSolid` under `solids {}`:

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
