# CasADi Migration Directives

## Version 3.8

CasADi 3.8 introduces NEP-18 (`__array_function__`) integration for CasADi types (`SX`, `MX`, `DM`), changing how CasADi objects interact with NumPy. Follow these directives when migrating code to CasADi 3.8:

### Avoid `numpy.poly1d` and `numpy.polyval` with symbolic types

- **Issue**: `numpy.poly1d(coeffs)(x)` internally delegates to `numpy.polyval(coeffs, x)`. In CasADi 3.7.2, symbolic inputs worked through Python duck-typing fallback. In CasADi 3.8, `SX` and `MX` implement `__array_function__`. Because `numpy.polyval` is not in CasADi's NEP-18 dispatch table, CasADi returns `NotImplemented`, causing NumPy to raise:

  ```
  TypeError: no implementation found for 'numpy.polyval' on types that implement __array_function__
  ```

- **Directive**: do not use `numpy.poly1d` or `numpy.polyval` on expressions that may be symbolic. Use a custom polynomial evaluator implementing Horner's method with standard arithmetic operators. Implement it such as that it works transparently with Python numeric scalars (`float`, `int`), NumPy `ndarray`, and CasADi `SX`/`MX` symbolic variables.

### Use CasADi Math Functions for Dual Numeric and Symbolic Evaluation

- **Issue**: calling NumPy mathematical functions (`np.exp`, `np.log10`, `np.sqrt`, `np.sin`, etc.) on CasADi symbolic values routes through NEP-18 dispatch and emits a `FutureWarning` in default NumPy mode (`mode = 0`).

- **Directive**: import and use native CasADi mathematical functions (`from casadi import exp, log, log10, sqrt, sin, cos`) in material properties and constitutive relations.

- **Numerical Compatibility**: all CasADi standard math functions natively accept Python numbers (`float`, `int`, `np.float64`) and return standard Python `float` values. This ensures that analytical expressions can be evaluated both numerically (e.g. calculating initial guesses or checking magnitudes) and symbolically (building DAE/ODE computational graphs) without branching.

### Global NumPy Mode Configuration (`GlobalOptions.setNumpyMode`)

CasADi 3.8 provides three operating modes via `casadi.GlobalOptions.setNumpyMode(mode)`:

  - `0` (default): Legacy 3.7.2 behavior, but emits a `FutureWarning` on the first NumPy operation called on a CasADi object.

  - `1`: Type-preserving CasADi-aware NumPy mode. NumPy functions return `NumpyArray` wrappers around CasADi expressions following NumPy shape and axis contracts.

  - `-1`: Legacy 3.7.2 behavior silently without emitting `FutureWarning`.

- **Directive**: Prefer writing code that does not depend on global mutable state by using native CasADi operations and operator-overloaded evaluators. If legacy NumPy behavior is strictly required across external packages, set `setNumpyMode(-1)` silently.
