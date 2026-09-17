# -*- coding: utf-8 -*-
import os
import shutil

from pathlib import Path


def get_blender():
    """ Resolve the Blender executable path from PATH or BLENDER_HOME.

    Checks if blender is already available on the system PATH. If not,
    falls back to the BLENDER_HOME environment variable, checking for
    typical executable names. If resolved, configures the local process
    environment and trimesh. If not found, raises a RuntimeError.

    Returns
    -------
    str
        Resolved path to the Blender executable.
    """
    if not (blender_exe := shutil.which("blender")):
        if (blender_home := os.environ.get("BLENDER_HOME")):
            for name in ["blender", "blender.exe"]:
                candidate = os.path.join(blender_home, name)
                if os.path.exists(candidate):
                    blender_exe = candidate
                    break

    if not blender_exe:
        raise RuntimeError(
            "Blender is required for robust CSG boolean operations but "
            "was not found. Please ensure 'blender' is in your PATH or "
            "set the 'BLENDER_HOME' environment variable pointing to "
            "your Blender installation directory."
        )

    blender_dir = str(Path(blender_exe).parent)

    if blender_dir not in os.environ["PATH"]:
        os.environ["PATH"] = blender_dir + os.pathsep + os.environ["PATH"]

    return blender_exe