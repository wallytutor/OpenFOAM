---
name: review-rust
description: Provides directives for reviewing Rust code.
---

# Coding in Rust

## When to use this skill

- Use this skill whenever a pull request, commit, or code snippet modifies `.rs` source files.

- Use this skill when reviewing the boundary layer where Rust structures and functions are exposed to Python via PyO3 (e.g., modules, functions, or structs marked with `#[pymodule]`, `#[pyfunction]`, or `#[pyclass]`).

- Use this skill to evaluate performance, memory safety, idiomatic patterns, and compilation efficiency of the Rust codebase.

- Do not use this skill for other programming languages, except when evaluating how changes to Rust code impact the generated Python bindings.

## How to use this skill

- All `use` directives come on top of a module.

- Avoid deeply nested code. Keep the code flat.

- **Zero-Warning Tolerance:** Analyze the code for compiler warnings, `clippy` lints, and documentation warnings. Do not suggest or accept solutions that suppress warnings using attributes like `#[allow(...)]` or `#![allow(...)]` unless explicitly told to do so by the user. Instead, identify the root cause of the warning (e.g., unused variables, dead code, non-idiomatic types) and provide the exact refactoring steps required to solve the underlying issue.

- **PyO3 Boundary Validation:** Ensure that all types crossing the Rust/Python boundary are safe and correctly map to Python types. Check that:

  - Errors are properly converted into `PyResult<T>` using meaningful Python exception types.

  - Memory management is sound when Python objects (`PyObject`, `Bound<'py, ...>`) are held or manipulated on the Rust side.

  - Python's Global Interpreter Lock (GIL) is respected, ensuring `Python::with_gil` is used appropriately or avoided in heavy parallel computation loops where `py.allow_threads` should be invoked.

- **Cargo & Compilation Speed Auditing:** Pay close attention to dependency usage and local build configurations. Ensure that any newly introduced code or features do not inadvertently cause global dependency cache invalidation or break parallel compilation.

- **Idiomatic Rust Code Review:** Verify adherence to the Rust API Guidelines, ensuring proper ownership management, correct usage of references, lifetimes, and efficient pattern matching.
