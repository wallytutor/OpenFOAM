---
name: data-lua
description: Provides directives for generating and modifying Lua code as data
---

# Coding databases in Lua

## When to use this skill

- Use this when writing or modifying Lua code as data.

- Do not use this skill for other programming languages.

## How to use this skill

- Databases must be readable to be maintainable and auditable.

- For all coefficients, use scientific format with at least 6 decimal digits and signed exponent, such as `1.000000e+01` and `-5.500000e-02`. Align values based on the dot, so that negative values are indented 4 spaces and positive ones 5 spaces.

- For all temperature range limits, use normal floatting point notation with at least 2 decimal digits, such as `298.15` and `2000.00`.

- Indent coefficients inside data constructors, breaking like for each value. Also indent elemental compositions, aligning equal signs (because some elements are 2, others 1-character long).

- Use 1 blank lines between local variables.
