---
name: tool-gmsh
description: Provides directives for using Gmsh scripting.
---

# Tool Gmsh

## Goal

Provide guidelines for using Gmsh to generate CAD geometry and mesh files, as applicable.

---

## When to use this skill

Use this skill when using package `gmsh` to generate CAD geometry and mesh files from Python scripts.

---

## How to use this skill

- Whenever in doubt, consult the provided references in the end of the document.

### Common misconceptions to avoid

- Agents tend to believe that STL files are *meshes* in the sense of a simulation mesh. They are not. STL files are just collections of triangular faces. Do not insist on using `gmsh.model.mesh.generate` before writing an STL file.

## References

- [Gmsh Docs](https://gmsh.info/doc/texinfo/gmsh.txt)
- [Gmsh API](https://gitlab.onelab.info/gmsh/gmsh/-/raw/gmsh_4_15_2/api/gmsh.py)
- [Gmsh Examples](https://gitlab.onelab.info/gmsh/gmsh/-/tree/master/examples?ref_type=heads)
