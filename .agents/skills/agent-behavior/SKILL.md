---
name: agent-behavior
description: Provides directives for agent behavior.
---

# Agent Behavior

## Goal

Standardize behavior of agents in any environment or repository.

---

## When to use this skill

Always use this skill as a preable to any other action.

---

## How to use this skill

Follow the guidelines described in the following points:

- If creation of test/check scripts is required as part of your task, always do so under a `scratch/` directory at the root of the repository to avoid leaving residuals. Do not delete them unless the user asks for it, as the contents might be useful for the user to check.

### About CHANGELOG.md files

- A changelog is a *top-bottom file*: latest modification on top, always as bullet items.

- A changelog is an additive stack, so previous entries should not be modified.

- Always leave one blank line between entries so that it is easy to read.
