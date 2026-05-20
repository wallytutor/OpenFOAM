---
name: code-python
description: Provides directives for generating and modifying Python code.
---

# Coding in Python

## When to use this skill

- Use this when writing, modifying, and reviewing Python code.

- Do not use this skill for other programming languages.

## How to use this skill

### General

- Use PEP8 for style guidelines, unless the user explicitly requests otherwise or manually modifies the code for readability.

- There must be a blanck line before flow control statements like `if`, `match`, `for`, `while`.

- Type hints must be used and the code must include type annotations for function arguments and return values.

- Functions must be documented using docstrings. The summary line in the docstring should be on the same line as `"""` and there must be at one whitespace character between the `"""` and the summary text.

- Docstrings should be in Numpydoc format. This means that docstrings must include the following sections:
    - Summary (mandatory)
    - Extended summary (optional)
    - Parameters (mandatory for functions)
    - Returns (mandatory for functions)
    - Raises (optional)
    - See also (optional)
    - Notes (optional)
