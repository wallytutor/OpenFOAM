# -*- coding: utf-8 -*-

import os

from pathlib import Path
from subprocess import run, STDOUT, PIPE
from typing import Callable


def banner(message: str) -> None:
    """ Minimal banner for printing messages in workflow. """
    print(f"{78 * '='}\n{message}\n{78 * '='}\n")


def source_openfoam_env(
        foam_root: str | Path = "/opt/openfoam13",
        shell: str = "bash"
    ) -> None:
    """ Sources OpenFOAM environment variables. """
    banner(f"Sourcing OpenFOAM environment for {shell}")
    rc = Path(foam_root) / f"etc/{shell}rc"

    if not rc.exists():
        raise FileNotFoundError(
            f"OpenFOAM environment file not found: {rc}"
        )

    # TODO check if in csh it is also source.
    args = [shell, "-c", f"source {rc} && env"]
    proc = run(args, stdout=PIPE, text=True, check=True)

    for line in proc.stdout.splitlines():
        key, _, val = line.partition("=")

        if key and val:
            os.environ[key] = val


class Runner:
    @staticmethod
    def log_file(
            log_name: str | None,
            app_name: str
        ) -> Path:
        """ Standardized name of a log for a given application. """
        return Path(log_name if log_name else f"log.{app_name}")

    @classmethod
    def serial(
            cls,
            args: list[str],
            log_name: str | None = None,
            force: bool = False,
        ) -> None:
        """ Wraps running an application with simultaneous logging. """
        log_file = cls.log_file(log_name, args[0])

        if log_file.exists() and not force:
            raise FileExistsError(log_file)

        with log_file.open("w") as f:
            run(args, stdout=f, stderr=STDOUT, check=True)

    @classmethod
    def parallel(
            cls,
            args: list[str],
            log_name: str | None = None,
            cores: int = 1,
            force: bool = False,
        ) -> None:
        """ Wraps running an application with simultaneous logging. """
        if cores > 1:
            cmd = ["mpirun", "-np", str(cores), args[0], "-parallel"]
            args = cmd + args[1:]

        cls.serial(args, log_name, force)


class Meshing:
    @classmethod
    def gmsh_to_foam_single_region(
            cls,
            mesh_file: str | Path,
            patching: Callable[[], None] | None = None
        ) -> None:
        """ Convert a Gmsh .msh into an OpenFOAM mesh. """
        banner("Workflow gmshToFoam for a single region")

        if not isinstance(mesh_file, Path):
            geometry = Path(mesh_file)

        if not geometry.exists():
            raise FileNotFoundError(geometry)

        Runner.serial(["gmshToFoam", str(geometry)])

        if callable(patching):
            patching()

        Runner.serial(["renumberMesh"])
        Runner.serial(["checkMesh"])

        banner("Finished!")

    @classmethod
    def snappyhexmesh(cls, *args, **kwargs) -> None:
        raise NotImplemented
