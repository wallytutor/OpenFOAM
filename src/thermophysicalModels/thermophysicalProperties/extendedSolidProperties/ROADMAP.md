# Dynamic Multi-Component Solid Thermophysics for `extendedThermoCloud`

## Overview

This document specifies the architectural implementation plan to enable multi-compound, temperature-dependent solid thermophysics within `extendedThermoCloud` for OpenFOAM 13.

The design:

1. Replaces the legacy, constant `solidProperties` container with an extensible, run-time selectable `extendedSolidProperties` hierarchy.
2. Introduces an abstract `solidMixturePropertiesModel` supporting arbitrary mixing laws (e.g., mass-weighted specific heat $C_p$, ideal volumetric density $\rho$).
3. Integrates with OpenFOAM's native `compositionModel singlePhaseMixture` and `singlePhaseMixtureCoeffs` inside `cloudProperties` to parse and track multi-solid mass fraction vectors $\mathbf{Y}$.
4. Replaces `parcelThermo` with `extendedParcelThermo` on the cloud level.
5. Purges redundant constant material entries (`Cp0`, `epsilon0`, `f0`) from `constantProperties` inside `cloudProperties`.

---

## Architectural Data Flow

```
constant/physicalProperties
  └── solids
        ├── solidMixtureModel  massWeighted;
        ├── sand    { type tabulatedSolid; Cp { type table; ... } }
        ├── alumina { type tabulatedSolid; Cp { type table; ... } }
        └── carbon  { type tabulatedSolid; Cp { type table; ... } }
              │
              ▼
extendedSolidMixtureProperties
  ├── PtrList<extendedSolidProperties>
  └── autoPtr<solidMixturePropertiesModel>
              │
              ▼
extendedParcelThermo (held by extendedThermoCloud)
              │
              ▼
constant/cloudProperties
  └── subModels
        └── compositionModel singlePhaseMixture;
        └── singlePhaseMixtureCoeffs
              └── phases ( solid { sand 0.5; alumina 0.3; carbon 0.2; } )
              │
              ▼
ExtendedThermoParcel
  └── calcHeatTransfer(...) / calc(...)
        └── evaluates Cp_mix(Y, p, T) via cloud.thermo().solids().Cp(Y, p, T)

```

---

## Phase 1: `extendedSolidProperties` Base & Derived Models

* [ ] **1.1 Create `extendedSolidProperties` Base Class**
* Inherit from `thermophysicalProperties`.
* Expose pure virtual, state-dependent getters:
```cpp
virtual scalar rho(const scalar p, const scalar T) const = 0;
virtual scalar Cp(const scalar p, const scalar T) const = 0;
virtual scalar kappa(const scalar p, const scalar T) const = 0;
virtual scalar emissivity(const scalar T) const = 0;
virtual scalar hf() const = 0;

```


* Declare standard run-time selection table infrastructure:
```cpp
TypeName("extendedSolidProperties");
declareRunTimeSelectionTable
(
    autoPtr,
    extendedSolidProperties,
    dictionary,
    (const dictionary& dict),
    (dict)
);

```




* [ ] **1.2 Implement `tabulatedSolidProperties**`
* Derive from `extendedSolidProperties`.
* Own `autoPtr<Function1<scalar>>` members for `Cp_`, `kappa_`, and `emissivity_`.
* Read `rho` and `hf` as scalar constants or optional `Function1<scalar>` entries.
* Register in the run-time selection table via `addToRunTimeSelectionTable`.


* [ ] **1.3 Implement `constantExtendedSolidProperties**`
* Derive from `extendedSolidProperties`.
* Read scalar values (`rho`, `Cp`, `kappa`, `Hf`, `emissivity`) directly from the sub-dictionary.
* Return invariant scalar values matching legacy OpenFOAM `solidProperties` behavior.



---

## Phase 2: Polymorphic `solidMixturePropertiesModel`

* [ ] **2.1 Implement `solidMixturePropertiesModel` Abstract Base Class**
* Declare run-time selection table keyed by word:
```cpp
TypeName("solidMixturePropertiesModel");
declareRunTimeSelectionTable
(
    autoPtr,
    solidMixturePropertiesModel,
    dictionary,
    (const dictionary& dict, const PtrList<extendedSolidProperties>& solids),
    (dict, solids)
);

```


* Store `const PtrList<extendedSolidProperties>& solids_`.
* Expose pure virtual mixing operations taking composition vector $\mathbf{Y}$:
```cpp
virtual scalar Cp(const scalarField& Y, const scalar p, const scalar T) const = 0;
virtual scalar rho(const scalarField& Y, const scalar p, const scalar T) const = 0;
virtual scalar kappa(const scalarField& Y, const scalar p, const scalar T) const = 0;
virtual scalar emissivity(const scalarField& Y, const scalar T) const = 0;
virtual scalar hf(const scalarField& Y) const = 0;

```




* [ ] **2.2 Implement `massWeightedSolidMixturePropertiesModel**`
* Derive from `solidMixturePropertiesModel`.
* Implement linear mass weighting for specific heat and formation enthalpy:

$$C_p(\mathbf{Y}, p, T) = \sum_{i} Y_i \, C_{p, i}(p, T), \quad h_f(\mathbf{Y}) = \sum_{i} Y_i \, h_{f, i}$$


* Implement ideal volumetric inverse mixing for density:

$$\frac{1}{\rho(\mathbf{Y}, p, T)} = \sum_{i} \frac{Y_i}{\rho_i(p, T)}$$


* Implement mass-weighted emissivity and thermal conductivity:

$$\varepsilon(\mathbf{Y}, T) = \sum_{i} Y_i \, \varepsilon_i(T), \quad \kappa(\mathbf{Y}, p, T) = \sum_{i} Y_i \, \kappa_i(p, T)$$


* Register in the run-time selection table.



---

## Phase 3: Container & Cloud Thermo Upgrade

* [ ] **3.1 Upgrade `extendedSolidMixtureProperties**`
* Maintain the solid component inventory and mixing engine:
```cpp
List<word> components_;
PtrList<extendedSolidProperties> properties_;
autoPtr<solidMixturePropertiesModel> mixtureModel_;

```


* In constructor, instantiate each component using `extendedSolidProperties::New(...)`.
* Parse the `solidMixtureModel` keyword (default: `"massWeighted"`) and construct `mixtureModel_`.
* Delegate `Cp(Y, p, T)`, `rho(Y, p, T)`, `kappa(Y, p, T)`, `emissivity(Y, T)`, and `hf(Y)` to `mixtureModel_`.


* [ ] **3.2 Implement `extendedParcelThermo**`
* Create `extendedParcelThermo` replacing standard `parcelThermo`.
* Replace:
```cpp
autoPtr<solidMixtureProperties> solids_;

```


with:
```cpp
autoPtr<extendedSolidMixtureProperties> solids_;

```


* In constructor, parse `solids_` from `carrierThermo.properties().subDict("solids")`.
* Expose `const extendedSolidMixtureProperties& solids() const`.


* [ ] **3.3 Update `extendedThermoCloud**`
* In `lagrangian/parcel/clouds/derived/extendedThermoCloud/`:
* Ensure the cloud holds `extendedParcelThermo thermo_`.
* Expose `const extendedParcelThermo& thermo() const`.





---

## Phase 4: Composition Integration & Parcel Dynamics

* [ ] **4.1 Connect Composition Mass Fractions $\mathbf{Y}$**
* In `ExtendedThermoParcel`, query the solid phase composition vector from the cloud's composition model:
```cpp
const label solidPhaseId = cloud.composition().idSolid();
const scalarField& Y = cloud.composition().Y0(solidPhaseId);

```


* For parcels with variable or independent composition, store and serialize `scalarField Y_` inside `ExtendedThermoParcelIO.C`.


* [ ] **4.2 Refactor `ExtendedThermoParcel::calcHeatTransfer**`
* Retrieve the mixture properties engine:
```cpp
const auto& solidMixture = cloud.thermo().solids();
const scalarField& Y = cloud.composition().Y0(cloud.composition().idSolid());

```


* **Step 1: Start-of-step state**
```cpp
Cp_ = solidMixture.Cp(Y, td.pc(), T_);
const scalar epsilon = solidMixture.emissivity(Y, T_);

```


* **Step 2: Predictor step**
```cpp
scalar CpEff = Cp_;
scalar bcp = htc*As/(m*CpEff);
scalar acp = bcp*td.Tc();
scalar ancp = Sh;
if (cloud.radiation())
{
    const tetIndices tetIs = this->currentTetIndices(td.mesh);
    const scalar Gc = td.GInterp().interpolate(this->coordinates(), tetIs);
    const scalar sigma = physicoChemical::sigma.value();
    ancp += As*epsilon*(Gc/4.0 - sigma*pow4(T_));
}
ancp /= (m*CpEff);

scalar deltaT = cloud.TIntegrator().delta(T_, dt, acp + ancp, bcp);

```


* **Step 3: Corrector step using midpoint evaluation**
```cpp
const scalar Tmid = max(cloud.constProps().TMin(), T_ + 0.5*deltaT);
CpEff = solidMixture.Cp(Y, td.pc(), Tmid);

bcp = htc*As/(m*CpEff);
acp = bcp*td.Tc();
ancp = (Sh + (cloud.radiation() ? As*epsilon*(td.GInterp().interpolate(this->coordinates(), this->currentTetIndices(td.mesh))/4.0 - physicoChemical::sigma.value()*pow4(T_)) : 0.0)) / (m*CpEff);

deltaT = cloud.TIntegrator().delta(T_, dt, acp + ancp, bcp);
const scalar deltaTncp = ancp*dt;
const scalar deltaTcp = deltaT - deltaTncp;

scalar Tnew = min(max(T_ + deltaT, cloud.constProps().TMin()), cloud.constProps().TMax());

```


* **Step 4: Source coupling and state update**
```cpp
dhsTrans -= m*CpEff*deltaTcp;
Sph = dt*m*CpEff*bcp;
Cp_ = solidMixture.Cp(Y, td.pc(), Tnew);
return Tnew;

```




* [ ] **4.3 Update `ExtendedThermoParcel::calc**`
* Bypass standard `cloud.composition().Cp(...)` override.
* Maintain the newly integrated $C_p$:
```cpp
this->Cp_ = cloud.thermo().solids().Cp(Y, td.pc(), this->T_);

```




* [ ] **4.4 Purge Redundant Values from `constantProperties**`
* In `ExtendedThermoParcel.H` and `ExtendedThermoParcelI.H`:
* Remove `Cp0_`, `epsilon0_`, and `f0_`.
* Retain operational bounds: `T0_`, `TMin_`, `TMax_`.
* Remove accessor functions `Cp0()`, `epsilon0()`, `f0()`.


* In parcel injection and initialization routines:
* Initialize initial parcel $C_p$ and $\rho$ via mixture evaluation at `T0`:
```cpp
Cp_ = cloud.thermo().solids().Cp(Y, pc, constProps.T0());
rho_ = cloud.thermo().solids().rho(Y, pc, constProps.T0());

```







---

## Phase 5: Build Configuration & Case Setup

* [ ] **5.1 Update `Make/files**`
Add source targets:
```make
extendedSolidProperties/extendedSolidProperties.C
extendedSolidProperties/tabulatedSolidProperties.C
extendedSolidProperties/constantExtendedSolidProperties.C
solidMixturePropertiesModel/solidMixturePropertiesModel.C
solidMixturePropertiesModel/massWeightedSolidMixturePropertiesModel.C
extendedSolidMixtureProperties/extendedSolidMixtureProperties.C
extendedParcelThermo/extendedParcelThermo.C

parcels/Templates/ExtendedThermoParcel/ExtendedThermoParcelName.C
parcels/derived/extendedThermoParcel/makeExtendedThermoParcelSubmodels.C
clouds/derived/extendedThermoCloud/extendedThermoCloud.C

LIB = $(FOAM_USER_LIBBIN)/libmyLagrangianParcel

```


* [ ] **5.2 Update `Make/options**`
Ensure include headers and link targets are set:
```make
EXE_INC = \
    -I$(LIB_SRC)/lagrangian/basic/lnInclude \
    -I$(LIB_SRC)/lagrangian/parcel/lnInclude \
    -I./parcels/Templates/ExtendedThermoParcel \
    -I./parcels/derived/extendedThermoParcel \
    -I./clouds/derived/extendedThermoCloud \
    -I./extendedSolidProperties \
    -I./solidMixturePropertiesModel \
    -I./extendedSolidMixtureProperties \
    -I./extendedParcelThermo \
    -I$(LIB_SRC)/finiteVolume/lnInclude \
    -I$(LIB_SRC)/meshTools/lnInclude \
    -I$(LIB_SRC)/thermophysicalModels/basic/lnInclude \
    -I$(LIB_SRC)/physicalProperties/lnInclude

LIB_LIBS = \
    -llagrangian \
    -llagrangianParcel \
    -lfiniteVolume \
    -lmeshTools \
    -lfluidThermo

```


* [ ] **5.3 Case Configuration: `constant/physicalProperties**`
```foam
solids
{
    solidMixtureModel   massWeighted;

    sand
    {
        type            tabulatedSolid;
        rho             2600;
        hf              -910700;
        emissivity      0.90;
        kappa           1.4;
        Cp
        {
            type        table;
            file        "$FOAM_CASE/constant/Cp_sand.csv";
            outOfBounds clamp;
        }
    }

    alumina
    {
        type            tabulatedSolid;
        rho             3950;
        hf              -1675700;
        emissivity      0.82;
        kappa           30.0;
        Cp
        {
            type        table;
            file        "$FOAM_CASE/constant/Cp_alumina.csv";
            outOfBounds clamp;
        }
    }
}

```


* [ ] **5.4 Case Configuration: `constant/cloudProperties**`
```foam
type    extendedThermoCloud;

constantProperties
{
    T0              300;
    TMin            200;
    TMax            3500;
}

subModels
{
    compositionModel singlePhaseMixture;

    singlePhaseMixtureCoeffs
    {
        phases
        (
            solid
            {
                sand        0.7;
                alumina     0.3;
            }
        );
    }

    heatTransferModel   RanzMarshall;
    radiation           off;
}

```