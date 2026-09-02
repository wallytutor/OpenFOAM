# -*- coding: utf-8 -*-

import os
import re
import shutil
import sys

from argparse import ArgumentParser
from pathlib import Path
from subprocess import run, STDOUT, PIPE
from typing import Callable


if sys.platform != "linux":
    raise OSError("Only Linux is supported!")


def common_arguments(source_env: bool = True) -> ArgumentParser:
    """ Reusable argument parsing for common OpenFOAM cases.

    By default it will source environment variables (controlled by
    `source_env` argument). If the defaults of `source_openfoam_env`
    do not suit the local installation, set that value to false and
    configure manually.
    """
    parser = ArgumentParser(description="OpenFOAM case workflow")

    parser.add_argument(
        "--cores",
        type    = int,
        default = 1,
        help    = "Number of cores to use"
    )
    parser.add_argument(
        "--clean",
        action = "store_true",
        help   = "Clean case files and logs"
    )

    #  TODO make --latest be accepted only if --reconstruct
    parser.add_argument(
        "--reconstruct",
        action = "store_true",
        help   = "Reconstruct parallel mesh/data"
    )
    parser.add_argument(
        "--latest",
        action = "store_true",
        help   = "Reconstruct only the latest time step"
    )

    action_group = parser.add_mutually_exclusive_group()
    action_group.add_argument(
        "--mesh",
        action = "store_true",
        help   = "Run meshing workflow"
    )
    action_group.add_argument(
        "--run",
        action = "store_true",
        help   = "Run solver workflow"
    )

    if source_env:
        source_openfoam_env()

    return parser


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


def get_latest_time(
        path: Path,
        regex: str = r"^[0-9]+(\.[0-9]+)?$"
    ) -> float | None:
    """ Find the latest execution time directory of a simulation. """
    times = []

    if not path.is_dir():
        return None

    for item in path.iterdir():
        if item.is_dir() and re.match(regex, item.name):
            try:
                times.append(float(item.name))
            except ValueError:
                pass

    return max(times) if times else None


def get_processor_dirs(root_dir) -> list[Path]:
    return [p for p in root_dir.glob("processor*") if p.is_dir()]


def is_openfoam_case(root_dir: Path | None = None) -> bool:
    if root_dir and not root_dir.exists():
        return False

    here = root_dir if root_dir else Path.cwd()
    return (here / "system/controlDict").exists()


def is_restart(
        cores: int,
        root_dir: Path | None = None,
    ) -> bool:
    """ Check if there is past data compatible with a restart. """
    here = root_dir if root_dir else Path.cwd()

    if not is_openfoam_case(here):
        raise FileNotFoundError(f"Not an OpenFOAM case: {here}")

    procs = get_processor_dirs(here)

    if len(procs) != cores:
        return False

    if (top_latest := get_latest_time(here)) is None:
        return False

    proc_times = [get_latest_time(p) for p in procs]

    return all(t is not None and t == top_latest for t in proc_times)


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

        cls.serial(args, log_name=log_name, force=force)

    @classmethod
    def decompose(
            cls,
            *,
            log_name: str | None = None,
            force: bool = False,
            patching: Callable[[], None] | None = None
        ) -> None:
        """ Manages the domain decomposition. """
        dict_file = Path("system/decomposeParDict")

        if callable(patching):
            patching()

        if not dict_file:
            raise FileNotFoundError(dict_file)

        cls.serial(["decomposePar"], log_name=log_name, force=force)

    @classmethod
    def reconstruct(
            cls,
            *,
            latest: bool = False,
            force_clean: bool = False,
            root_dir: str | Path | None = None,
        ) -> None:
        here = root_dir if root_dir else Path.cwd()
        procs = get_processor_dirs(here)

        if not procs:
            return

        if latest:
            cls.serial(["reconstructPar", "-latestTime"])
        else:
            cls.serial(["reconstructPar", "-constant"])

        if force_clean:
            for proc in procs:
                shutil.rmtree(proc, ignore_errors=True)

    @classmethod
    def foam_run(
            cls,
            *,
            log_name: str | None = None,
            cores: int = 1,
            force: bool = False,
            reconstruct: bool = False,
            preprocess: Callable[[], None] | None = None,
            decomposing: Callable[[], None] | None = None,
            **kwargs
        ) -> None:
        restart = is_restart(cores)

        if not restart and callable(preprocess):
            preprocess()

        if not restart and cores > 1:
            cls.decompose(patching=decomposing)

        cls.parallel(
            args     = ["foamRun"],
            log_name = log_name,
            cores    = cores,
            force    = force
        )

        if reconstruct:
            pass
            # _reconstruct(latest=latest)

    @classmethod
    def dict_set_entry(
            cls,
            file: str | Path,
            entry: str,
            value: str,
            *,
            log_name: str | None = None,
            force: bool = False
        ) -> None:
        # If not log name is provided, assume the user is ok with
        # overwritting (as setting entries is generally a batch).
        if not log_name:
            force = True

        Runner.serial(
            args = [
                "foamDictionary", file,
                "-entry", entry,
                "-set", value
            ],
            force = force
        )


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
