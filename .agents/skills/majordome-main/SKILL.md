---
name: majordome-main
description: Provides directives for building and maintaining the Majordome package.
---

# Majordome Main

## Goal

Provide *repository-aware* implementation and maintenance instructions.

---

## When to use this skill

- Modifications to the Majordome project (python-majordome).

- Use of Majordome within other projects that depend on its structure.

---

## Project Map

- Packaging metadata: `pyproject.toml`, `Cargo.toml`
- Python package: `majordome/`
- Rust extension core: `src/`
- Quarto docs: `docs/`

---

## Repository Guardrails

- Prefer minimal patches over broad rewrites.
- Preserve existing style and layout.
- Add comments only where logic is non-obvious
- Prefer documenting over commenting code.
- Never edit generated docs by default (`docs/_book/`).
- Run smallest relevant verification command after edits.
- If a tool/command is missing, report that and provide an equivalent path.

---

## Task Routing

- Python utility/library changes -> Python Package workflow
- CPython extension or import/runtime failures -> Rust/PyO3 workflow
- Documentation/tutorial/site edits -> Quarto Documentation workflow
- Versioning/changelog/release changes -> Release workflow

---

## Python Package Workflow

### Use For

- Edits under `majordome/*.py`
- Entrypoints and script wiring
- Python dependency or metadata updates

### Key Files

- `pyproject.toml`
- `majordome/__init__.py`
- `majordome/entrypoints.py`
- `majordome/**/test.py`

### Steps

1. Read package metadata and scripts in `pyproject.toml`.
2. Keep changes targeted and API-compatible unless API change is requested.
3. Update tests when behavior changes.
4. Verify script references if entrypoints move.

### Validate

```bash
pytest -q
```

Entrypoints currently expected, see `[project.scripts]` in `pyproject.toml`.

---

## Rust/PyO3 Workflow

### Use For

- Changes in `src/core/`
- PyO3 bindings and build behavior
- Extension import/load issues

### Key Files

- `Cargo.toml`
- `build.rs`
- `src/core/lib.rs`
- `pyproject.toml` (`[tool.maturin]`)

### Steps

1. Keep module mapping consistent:
   - Rust crate/lib: `_core`
   - Python import target: `majordome._core`
2. Keep ABI target coherent unless compatibility change is requested (`abi3-py312`).
3. Preserve `build.rs` runtime staging behavior on Windows.

### Validate

```bash
cargo check
cargo test
```

Optional local extension check:

```bash
maturin develop
python -c "import majordome, majordome._core; print('ok')"
```

### Windows Pitfall

If Rust tests fail with `STATUS_DLL_NOT_FOUND (0xc0000135)`, verify `build.rs` still copies `pythonXY.dll` into `target/<profile>/` and `target/<profile>/deps/` based on the active interpreter.

---

## Quarto Documentation Workflow

### Use For

- Source documentation edits in `.qmd` files

### Key Files

- `docs/_quarto.yaml`
- `docs/*.qmd`
- `devel/*.qmd`
- `reactors/*.qmd`
- `diffusion/*.qmd`

### Steps

1. Edit source `.qmd` files.
2. Do not manually edit generated `docs/_book/*` unless explicitly requested.
3. Keep equations/links/code snippets coherent.

### Validate

```bash
uv run majordome-build-qmd
```

---

## Release and Dependency Workflow

### Use For

- Version bumps
- Dependency updates
- Changelog/release prep

### Key Files

- `pyproject.toml`
- `Cargo.toml`
- `CHANGELOG.md`
- `release.py`, `release.ps1`, `release.sh`
- `Containerfile`

### Steps

1. Keep Python and Rust version metadata aligned.
2. Update changelog with user-visible changes.
3. Prefer minimal dependency expansion and justify heavy additions.

### Validate

```bash
cargo check
pytest -q
```
