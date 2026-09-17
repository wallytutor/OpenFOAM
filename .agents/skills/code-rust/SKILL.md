---
name: code-rust
description: Provides directives for generating and modifying Rust code.
---

# Coding in Rust

## When to use this skill

- Use this when writing or modifying Rust code

- Use this skill when writting PyO3 extensios for Python modules.

- Do not use this skill for other programming languages.

## How to use this skill

### Crates structure and management

- Rust crates live under the `src/` or `crates/` directory depending on project structure and each crate has its own `Cargo.toml` file.

- Crates living under `crates/` have the `lib.rs` at the root of the crate directory (same as `Cargo.toml`) to avoid unnecessary nesting of directories in the tree.

- Public code exported from submodule and intended for use in Rust must be always explicitly listed in `prelude.rs` files so that user code does not import from deep within the tree. For instance:

    ```rust
    // Expecting the user to write the following is not allowed:
    use majordome_numerical::linear_algebra::tridiagonal::TridiagonalSolver;

    // Instead, one of the following should be possible to use:
    use majordome_numerical::prelude::TridiagonalSolver;
    use majordome_numerical::prelude::*
    ```

- The path inside the crate must be as short as possible. So, if a submodule is not meant to be used directly, it should not be exposed in `prelude.rs`.

- Not all crates are meant to be used directly. Some are internal implementation details.

- Always place new Rust code in the most appropriate file and create new modules when needed.

- Include integration tests in `test/` when appropriate (*e.g.* for high-level functionality) or in a `test` submodule (singular, not `tests`).

- Each module must be fully unit tested with tests under the `#[cfg(test)]` attribute for its own functionalities. Crates may use more complex integration tests as discussed above.

- Ensure all tests pass with `cargo test` before proposing a change.

- Ensure modules are picky with code style and use `cargo clippy` for linting. All new code MUST pass `cargo clippy` and include the following on top:

    ```rust
    #![deny(clippy::unwrap_used)]
    #![deny(clippy::expect_used)]
    #![deny(clippy::panic)]
    #![deny(unused_must_use)]
    #![deny(warnings)]
    ```

### Code style

- Add documentation comments (`///`) and examples for all public APIs. Documentation lines are limited to 78 characters, including the leading triple slashes. Example:

    ```rust
    /// This is a short doc comment.
    /// It is very short indeed.
    ///
    /// The lines are not supposed to exceed 78 characters, including the leading slashes (so this line is not acceptable).
    ///
    /// This is a new paragraph.
    ```

- Although you are coding in Rust, the author has Irlen syndrome and needs code to be extremely readable.

- Instead of using Rust typical guidelines for style, stick to PEP8 adapted for Rust whenever possible.

- Use a blank line before if, match, for, while, and loop keywords. Use a blank line between functions and blank lines between methods in a struct/enum. Inside `impl` and `[#pymethods] impl` blocks use a single blank line between functions and methods.

- Whenever reasonable, allow aligning `=>` in matches so that it remains *tabular*.

- For any non-trivial match, use braces and indent the code block (no-inline).

- Lines are preferrably less than 80 characters long, with a maximum of 100 characters, except for documentation as discussed above.

### Python bindings with PyO3

- Do not use direct *in-struct* attribute PyO3 annotations, i.e., keep the python API logic completely disconnected from Rust structure. All elements to be exposed shall be given a dedicated interface in the `[#pymethods] impl` block.

- Do not use attributes to bypass any warning, such as `skip_from_py_object`. Instead, implement the new recommended approach so that the code base is already compliant with new standards.

- When creating Python API's, do not use `py_` naming schemes. Any created API must look native as they are the primary goal. Find alternatives in the code structure.
