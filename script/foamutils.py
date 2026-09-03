# -*- coding: utf-8 -*-

import functools
import os
import re
import shlex
import shutil
import sys

from argparse import ArgumentParser
from pathlib import Path
from subprocess import run, STDOUT, PIPE
from time import perf_counter, time_ns
from typing import Any, Callable, Sequence


if sys.platform != "linux":
    raise OSError("Only Linux is supported!")


TIME_DIR_REGEX = r"^[0-9]+(\.[0-9]+)?([eE][+-]?[0-9]+)?$"


class ArgumentsPresets:
    """ Reusable ArgumentParser preset builders for OpenFOAM workflows. """

    __slots__ = ()

    @classmethod
    def common(cls, source_env: bool = True) -> ArgumentParser:
        """ Construct reusable argument parser for OpenFOAM cases.

        Parameters
        ----------
        source_env : bool = True
            Whether to source OpenFOAM environment setup variables.

        Returns
        -------
        ArgumentParser
            Configured argument parser instance for OpenFOAM CLI workflows.
        """
        parser = ArgumentParser(description="OpenFOAM case workflow")

        parser.add_argument(
            "--cores",
            type    = int,
            default = 1,
            help    = "Number of cores to use"
        )

        ###
        # Case cleaning
        ###

        parser.add_argument(
            "--clean",
            action = "store_true",
            help   = "Clean case files and logs"
        )
        parser.add_argument(
            "--clean-logs",
            action = "store_true",
            help   = "Clean log files only"
        )
        parser.add_argument(
            "--clean-processors",
            action = "store_true",
            help   = "Clean processor directories only"
        )

        # TODO make --latest be accepted only if --reconstruct
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

        ###
        # Main workflows
        ###

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

    @classmethod
    def meshing(cls, source_env: bool = True) -> ArgumentParser:
        """ Construct reusable argument parser for OpenFOAM cases.

        Parameters
        ----------
        source_env : bool = True
            Whether to source OpenFOAM environment setup variables.

        Returns
        -------
        ArgumentParser
            Configured argument parser instance for OpenFOAM CLI workflows.
        """
        parser = cls.common(source_env=source_env)

        parser.add_argument(
            "--renumber-mesh",
            action = "store_true",
            help   = "Renumber mesh"
        )
        parser.add_argument(
            "--check-mesh",
            action = "store_true",
            help   = "Check mesh quality"
        )
        # parser.add_argument(
        #     "--check-mesh-options",
        #     nargs = "*",
        #     default = [],
        #     help   = "Options for checkMesh tool"
        # )

        return parser

    @classmethod
    def snappyhexmesh(cls, source_env: bool = True) -> ArgumentParser:
        """ Construct reusable argument parser for OpenFOAM cases.

        Parameters
        ----------
        source_env : bool = True
            Whether to source OpenFOAM environment setup variables.

        Returns
        -------
        ArgumentParser
            Configured argument parser instance for OpenFOAM CLI workflows.
        """
        parser = cls.meshing(source_env=source_env)
        return parser


def banner(message: str) -> None:
    """ Print a formatted banner message to standard output.

    Parameters
    ----------
    message : str
        Message text to display inside standard workflow banner.

    Returns
    -------
    None
        Output is written directly to standard stdout stream.
    """
    print(f"{78 * '='}\n> {message}\n")


def validate_input(msg: str) -> bool:
    while True:
        match (ans := input(f"{msg} (y/N): ").lower().strip()):
            case "y":
                return True
            case "n":
                return False
            case _:
                print(f"Invalid input {ans}; please, answer (y/N)")
                continue


def source_openfoam_env(
        foam_root: str | Path = "/opt/openfoam13",
        shell: str = "bash"
    ) -> None:
    """ Source OpenFOAM environment variables into current environment.

    Parameters
    ----------
    foam_root : str | Path = "/opt/openfoam13"
        Root installation path of OpenFOAM distribution.
    shell : str = "bash"
        Shell executable used to run environment configuration script.

    Returns
    -------
    None
        Updates active environment variables in os.environ directly.
    """
    # Do not source if already sourced, it's slow...
    wm_project_dir = os.environ.get("WM_PROJECT_DIR", None)
    if wm_project_dir and Path(wm_project_dir).resolve() == Path(foam_root).resolve():
        return

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
        regex: str = TIME_DIR_REGEX
    ) -> float | None:
    """ Find numerical value of latest simulation time directory.

    Parameters
    ----------
    path : Path
        Target directory containing execution time step folders.
    regex : str = TIME_DIR_REGEX
        Regular expression matching valid time directory names.

    Returns
    -------
    float | None
        Highest numerical time step found, or None if empty.
    """
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


def get_processor_dirs(
        root_dir: Path | None = None
    ) -> list[Path]:
    """ Return processor directories.

    Parameters
    ----------
    root_dir : Path | None = None
        Case directory path to check. Defaults to current working dir.

    Returns
    -------
    list[Path]
        List of processor directory paths.
    """
    here = root_dir if root_dir else Path.cwd()
    return [p for p in here.glob("processor*") if p.is_dir()]


def is_openfoam_case(root_dir: Path | None = None) -> bool:
    """ Check if target path contains a valid OpenFOAM case.

    Parameters
    ----------
    root_dir : Path | None = None
        Case directory path to verify. Defaults to current working dir.

    Returns
    -------
    bool
        True if system/controlDict file is present.
    """
    if root_dir and not root_dir.exists():
        return False

    here = root_dir if root_dir else Path.cwd()
    return (here / "system/controlDict").exists()


def ensure_openfoam_case(func: Callable) -> Callable:
    """ Decorate workflow function to ensure OpenFOAM case context.

    Parameters
    ----------
    func : Callable
        Workflow entrypoint function to wrap.

    Returns
    -------
    Callable
        Wrapped workflow execution function.
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        t0 = perf_counter()
        here = kwargs.get("root_dir", Path.cwd())

        if not is_openfoam_case(here):
            raise NotADirectoryError(f"Not an OpenFOAM case: {here}")

        Path("case.foam").touch()

        values = func(*args, **kwargs)

        print(f"Workflow took {perf_counter() - t0} seconds")
        return values

    return wrapper


def is_restart(
        cores: int,
        root_dir: Path | None = None
    ) -> bool:
    """ Verify if directory contains valid past data for simulation restart.

    Parameters
    ----------
    cores : int
        Number of processor subdomains configured for parallel run.
    root_dir : Path | None = None
        Case directory path to check. Defaults to current working dir.

    Returns
    -------
    bool
        True if consistent time step outputs exist across all domains.
    """
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


def clean_times(
        root_dir: Path | None = None,
        *,
        remove_zero: bool = True
    ) -> None:
    """ Clean simulation time step output directories.

    Parameters
    ----------
    root_dir : Path | None = None
        Case directory path to clean. Defaults to current working dir.
    remove_zero : bool = True
        Whether to remove the initial conditions zero directory.

    Returns
    -------
    None
        Time step directories are removed in-place.
    """
    cwd = Path(root_dir) if root_dir else Path.cwd()

    for item in cwd.iterdir():
        if item.is_dir() and re.match(TIME_DIR_REGEX, item.name):
            try:
                t = float(item.name)
            except ValueError:
                continue

            if t > 0.0 or (t == 0.0 and remove_zero):
                shutil.rmtree(item, ignore_errors=True)


def clean_logs(
        root_dir: Path | None = None,
    ) -> None:
    """ Clean execution log files matching log.* pattern.

    Parameters
    ----------
    root_dir : Path | None = None
        Case directory path to clean. Defaults to current working dir.

    Returns
    -------
    None
        Log files matching log.* are removed in-place.
    """
    cwd = Path(root_dir) if root_dir else Path.cwd()

    for log_file in cwd.glob("log.*"):
        if log_file.is_file():
            log_file.unlink(missing_ok=True)


def clean_processors_dirs(
        root_dir: Path | None = None,
    ) -> None:
    """ Remove all processor directories from a case.

    Parameters
    ----------
    root_dir : Path | None = None
        Case directory path to clean. Defaults to current working dir.

    Returns
    -------
    None
        Processor directories are removed in-place.
    """
    for p in get_processor_dirs(root_dir):
        if p.is_dir():
            shutil.rmtree(p, ignore_errors=True)


def clean_case(
        root_dir: str | Path | None = None,
        *,
        remove_zero: bool = True,
        extra_dirs: Sequence[str | Path] | None = None,
        extra_files: Sequence[str | Path] | None = None,
        extra_patterns: Sequence[str] | None = None,
        **kwargs: Any,
    ) -> None:
    """ Clean execution outputs, mesh directories, and log files in case.

    Parameters
    ----------
    root_dir : str | Path | None = None
        Target case root directory. Defaults to current working dir.
    remove_zero : bool = True
        Whether to remove the 0 initial conditions directory.
    extra_dirs : Sequence[str | Path] | None = None
        Additional directory paths to remove.
    extra_files : Sequence[str | Path] | None = None
        Additional file paths to remove.
    extra_patterns : Sequence[str] | None = None
        Glob patterns for matching files/directories to delete.
    **kwargs : Any
        Additional case cleaning options and arguments.

    Returns
    -------
    None
        Case workspace is cleaned in place.
    """
    cwd = Path(root_dir) if root_dir else Path.cwd()

    if not is_openfoam_case(cwd):
        raise FileNotFoundError(f"Not an OpenFOAM case: {cwd}")

    clean_times(cwd, remove_zero=remove_zero)
    clean_processors_dirs(cwd)
    clean_logs(cwd)

    to_remove = [
        cwd / "constant" / "extendedFeatureEdgeMesh",
        cwd / "constant" / "polyMesh",
        cwd / "postProcessing",
    ]

    if extra_dirs:
        for ed in extra_dirs:
            to_remove.append(Path(ed) if not isinstance(ed, Path) else ed)

    for p in to_remove:
        if p.is_dir():
            shutil.rmtree(p, ignore_errors=True)

    if (geom_dir := cwd / "constant" / "geometry").is_dir():
        for f in geom_dir.glob("*.eMesh"):
            if f.is_file():
                f.unlink(missing_ok=True)

    if (case_foam := cwd / "case.foam").exists():
        case_foam.unlink(missing_ok=True)

    if extra_files:
        for ef in extra_files:
            p = Path(ef) if not isinstance(ef, Path) else ef

            if p.is_file():
                p.unlink(missing_ok=True)

    if extra_patterns:
        for pattern in extra_patterns:
            for match in cwd.glob(pattern):
                if match.is_file():
                    match.unlink(missing_ok=True)
                elif match.is_dir():
                    shutil.rmtree(match, ignore_errors=True)


class Runner:
    """ Manage serial and parallel execution of OpenFOAM applications. """

    __slots__ = ()

    @staticmethod
    def log_file(log_name: str | None, app_name: str) -> Path:
        """ Return standardized log file Path for an application.

        Parameters
        ----------
        log_name : str | None
            Explicit log filename override, or None for standard format.
        app_name : str
            Executable application name used to generate default filename.

        Returns
        -------
        Path
            Target log file path location.
        """
        return Path(log_name if log_name else f"log.{app_name}")

    @classmethod
    def serial(
            cls,
            args: str | list[str],
            *,
            log_name: str | None = None,
            force: bool = False
        ) -> None:
        """ Execute application in serial mode while logging stdout/stderr.

        Parameters
        ----------
        args : str | list[str]
            Command arguments array starting with application binary name.
        log_name : str | None = None
            Custom log filename override.
        force : bool = False
            Whether to overwrite an existing log file of the same name.

        Returns
        -------
        None
            Subprocess executes and directs output into log file.
        """
        if isinstance(args, str):
            args = shlex.split(args)

        log_file = cls.log_file(log_name, args[0])

        if log_file.exists() and not force:
            raise FileExistsError(log_file)

        with log_file.open("w") as f:
            run(args, stdout=f, stderr=STDOUT, check=True)

    @classmethod
    def parallel(
            cls,
            args: str | list[str],
            *,
            log_name: str | None = None,
            cores: int = 1,
            force: bool = False
        ) -> None:
        """ Execute application in parallel mode using mpirun launcher.

        Parameters
        ----------
        args : str | list[str]
            Command arguments array for OpenFOAM solver/tool.
        log_name : str | None = None
            Custom log filename override.
        cores : int = 1
            Number of MPI process slots to allocate.
        force : bool = False
            Whether to overwrite existing log file.

        Returns
        -------
        None
            Constructs MPI execution command and delegates to serial log.
        """
        if isinstance(args, str):
            args = shlex.split(args)

        app_bin = str(args[0])

        if cores > 1:
            cmd = ["mpirun", "-np", str(cores), app_bin, "-parallel"]
            args = cmd + args[1:]

        log_name = cls.log_file(
            log_name,
            f"{app_bin}_parallel" if cores > 1 else app_bin
        ).name

        cls.serial(args, log_name=log_name, force=force)

    @classmethod
    def batch(
            cls,
            args: list[list[str]],
        ) -> None:
        """ Execute application in serial mode while logging stdout/stderr.

        Parameters
        ----------
        args : list[list[str]]
            List of command arguments arrays, each starting with application binary name. Commands are run in sequence and logs use the time
            stamp as unique identifiers.

        Returns
        -------
        None
            Subprocess executes and directs output into log file.
        """
        for arg in args:
            log_name = f"log.batch_{time_ns()}"
            cls.serial(arg, log_name=log_name, force=True)

    @classmethod
    def decompose(
            cls,
            *,
            log_name: str | None = None,
            cores: int = 1,
            force: bool = False,
            patching: Callable[[int], None] | None = None
        ) -> None:
        """ Run domain decomposition using decomposePar utility.

        Parameters
        ----------
        log_name : str | None = None
            Custom log filename for decomposePar command.
        cores : int = 1
            Number of processor slots to allocate.
        force : bool = False
            Whether to overwrite existing decomposePar log.
        patching : Callable[[int], None] | None = None
            Optional setup function executed before domain decomposition.
            Takes the number of cores as an argument.

        Returns
        -------
        None
            Invokes decomposePar tool after optional patching step.
        """
        if cores < 2:
            return

        if is_restart(cores) and not force:
            return

        dict_file = Path("system/decomposeParDict")

        if callable(patching):
            patching(cores)

        if not dict_file.exists():
            raise FileNotFoundError(dict_file)

        cls.serial(["decomposePar"], log_name=log_name, force=force)

    @classmethod
    def reconstruct(
            cls,
            *,
            log_name: str | None = None,
            force: bool = False,
            latest: bool = False,
            constant: bool = False,
            options: list[str] | None = None,
        ) -> None:
        """ Reconstruct parallel simulation data using reconstructPar.

        Parameters
        ----------
        log_name : str | None = None
            Custom log filename for reconstructPar output.
        force : bool = False
            Whether to overwrite existing reconstructPar log file.
        latest : bool = False
            Whether to reconstruct only the latest time step.
        constant : bool = False
            Whether to reconstruct constant directory data.
        options : list[str] | None = None
            Additional reconstructPar command line options.

        Returns
        -------
        None
            Runs reconstructPar in serial mode if processor dirs exist.
        """
        procs = [p for p in Path.cwd().glob("processor*") if p.is_dir()]

        if not procs:
            return


        options = options or []

        if latest:
            options.append("-latestTime")

        if constant:
            options.append("-constant")

        cmd = ["reconstructPar"] + list(set(options))
        cls.serial(cmd, log_name=log_name, force=force)

    @classmethod
    def surface_features(
            cls,
            log_name: str | None = None,
            force: bool = False,
        ) -> None:
        """ Run surfaceFeatures utility on STL surfaces.

        Parameters
        ----------
        log_name : str | None = None
            Custom log filename for surfaceFeatures command.
        force : bool = False
            Whether to force execution even if eMesh files exist.

        Returns
        -------
        None
            Executes surfaceFeatures tool when needed.
        """
        stl_files = Path("constant/geometry").glob("*.stl")

        if (
            force or
            any(not f.with_suffix(".eMesh").exists() for f in stl_files)
        ):
            cls.serial("surfaceFeatures", log_name=log_name, force=force)
            return

        print("Skipping surfaceFeatures - eMesh files already exist.")

    @classmethod
    def foam_run(
            cls,
            app_args: list[str] | None = None,
            *,
            log_name: str | None = None,
            cores: int = 1,
            force: bool = False,
            reconstruct: bool = False,
            preprocess: Callable[[], None] | None = None,
            decomposing: Callable[[], None] | None = None,
            latest: bool = False
        ) -> None:
        """ Manage end-to-end OpenFOAM solver workflow execution.

        Parameters
        ----------
        app_args : list[str] | None = None
            Application command list. Defaults to ['foamRun'].
        log_name : str | None = None
            Log filename for solver execution.
        cores : int = 1
            Number of processor slots to use for computation.
        force : bool = False
            Whether to force overwriting existing log outputs.
        reconstruct : bool = False
            Whether to run reconstruction after solver finishes.
        preprocess : Callable[[], None] | None = None
            Callback executed before running initial setup.
        decomposing : Callable[[], None] | None = None
            Callback executed during domain decomposition phase.
        latest : bool = False
            Whether to reconstruct only the latest time step.

        Returns
        -------
        None
            Executes preprocessing, decomposition, solver, and optional rec.
        """
        if not is_restart(cores) and callable(preprocess):
            preprocess()

        cls.decompose(cores=cores, patching=decomposing)

        cmd_args = app_args if app_args is not None else ["foamRun"]

        cls.parallel(
            args     = cmd_args,
            log_name = log_name,
            cores    = cores,
            force    = force
        )

        if reconstruct:
            cls.reconstruct(latest=latest)

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
        """ Set a dictionary entry value via foamDictionary tool.

        Parameters
        ----------
        file : str | Path
            Target OpenFOAM dictionary file path.
        entry : str
            Target key/entry path within dictionary.
        value : str
            New string value to write to entry.
        log_name : str | None = None
            Custom log filename for foamDictionary call.
        force : bool = False
            Whether to force log file replacement.

        Returns
        -------
        None
            Executes foamDictionary -entry -set command in serial mode.
        """
        if not log_name:
            force = True

        cls.serial(
            args = [
                "foamDictionary", str(file),
                "-entry", entry,
                "-set", value
            ],
            log_name = log_name,
            force    = force
        )


class Meshing:
    """ Utility routines for mesh generation and conversion operations. """

    __slots__ = ()

    @staticmethod
    def post_meshing(
            renumber_mesh: bool = True,
            check_mesh: bool = True,
            **kwargs: Any
        ) -> None:
        """ Run post-meshing routines (renumberMesh, checkMesh).

        Parameters
        ----------
        renumber_mesh : bool = True
            Whether to renumber mesh to reduce bandwidth.
        check_mesh : bool = True
            Whether to run checkMesh for quality evaluation.
        **kwargs : Any
            Additional execution options (e.g., checkMesh_options).

        Returns
        -------
        None
            Executes post-meshing utilities.
        """
        if renumber_mesh:
            Runner.serial("renumberMesh")

        if check_mesh:
            opts = kwargs.get("checkMesh_options", [])
            Runner.serial(["checkMesh"] + opts)

    @classmethod
    def gmsh_to_foam_single_region(
            cls,
            mesh_file: str | Path,
            *,
            renumber_mesh: bool = True,
            check_mesh: bool = True,
            patching: Callable[[], None] | None = None,
            **kwargs: Any
        ) -> None:
        """ Convert Gmsh mesh into single-region OpenFOAM polyMesh format.

        Parameters
        ----------
        mesh_file : str | Path
            Path to input Gmsh .msh mesh file.
        renumber_mesh : bool = True
            Whether to renumber mesh to reduce matrix bandwidth.
        check_mesh : bool = True
            Whether to run checkMesh for mesh quality evaluation.
        patching : Callable[[], None] | None = None
            Optional callback executed after conversion to update boundaries.
        **kwargs : Any
            Additional execution options (e.g., checkMesh_options).

        Returns
        -------
        None
            Runs gmshToFoam, renumberMesh, and checkMesh sequentially.
        """
        banner("Workflow gmshToFoam for a single region")

        geometry = Path(mesh_file)

        if not geometry.exists():
            raise FileNotFoundError(geometry)

        Runner.serial(["gmshToFoam", str(geometry)])

        if callable(patching):
            patching()

        cls.post_meshing(renumber_mesh, check_mesh, **kwargs)

    @classmethod
    def snappyhexmesh(
            cls,
            *,
            cores: int = 1,
            force: bool = False,
            extract_surface_features: bool = True,
            clean_extended_features: bool = True,
            create_background_mesh: bool = True,
            reconstruct: bool = True,
            clean_parallel_dirs: bool = True,
            renumber_mesh: bool = True,
            check_mesh: bool = True,
            geometry: Callable[[int], None] | None = None,
            preprocess: Callable[[], None] | None = None,
            decomposing: Callable[[int], None] | None = None,
            postprocess: Callable[[], None] | None = None,
            **kwargs: Any,
        ) -> None:
        """ Generate mesh using snappyHexMesh utility workflow.

        Parameters
        ----------
        cores : int = 1
            Number of processor slots to use for computation.
        force : bool = False
            Whether to force overwriting existing log outputs.
        extract_surface_features : bool = True
            Whether to run surfaceFeatures prior to meshing.
        clean_extended_features : bool = True
            Whether to clean constant/extendedFeatureEdgeMesh.
        create_background_mesh : bool = True
            Whether to create background mesh via blockMesh. To force it
            running you must manually clean the case before running; if
            constant/polyMesh already exists, blockMesh will be skipped.
            This is intended to allow the different snappyHexMesh steps
            to be run independently.
        reconstruct : bool = True
            Whether to reconstruct mesh after parallel execution.
        clean_parallel_dirs : bool = True
            Whether to remove processor dirs post-reconstruct.
        renumber_mesh : bool = True
            Whether to renumber mesh to reduce bandwidth.
        check_mesh : bool = True
            Whether to run checkMesh for quality evaluation.
        geometry : Callable[[int], None] | None = None
            Callback for retrieving or generating geometry files.
            Takes the number of cores as an argument.
        preprocess : Callable[[], None] | None = None
            Callback executed before initial mesh generation.
        decomposing : Callable[[int], None] | None = None
            Callback executed during domain decomposition phase.
            Takes the number of cores as an argument.
        postprocess : Callable[[], None] | None = None
            Callback executed after mesh reconstruction.
        **kwargs : Any
            Additional execution options and arguments.

        Returns
        -------
        None
            Executes complete snappyHexMesh workflow.
        """
        banner("Workflow snappyHexMesh")

        if callable(geometry):
            geometry(cores)

        if extract_surface_features:
            Runner.surface_features(force=False)

        if clean_extended_features:
            if (feat_dir := Path("constant/extendedFeatureEdgeMesh")).is_dir():
                shutil.rmtree(feat_dir, ignore_errors=True)

        if create_background_mesh:
            if Path("constant/polyMesh").exists():
                print("Skipping blockMesh — polyMesh already exists.")
            else:
                opts = kwargs.get("blockMesh_options", [])
                Runner.serial(["blockMesh"] + opts)

        if callable(preprocess):
            preprocess()

        Runner.decompose(cores=cores, patching=decomposing)

        cmd = ["snappyHexMesh"] + kwargs.get("snappyHexMesh_options", [])

        Runner.parallel(
            args     = cmd,
            cores    = cores,
            force    = force
        )

        if reconstruct:
            Runner.reconstruct(constant=True)

        if clean_parallel_dirs:
            clean_processors_dirs()

        if callable(postprocess):
            postprocess()

        cls.post_meshing(renumber_mesh, check_mesh, **kwargs)


class CommonProjectManager:
    """ Manage OpenFOAM case project execution workflows.

    Parameters
    ----------
    root_dir : Path
        Base directory of the OpenFOAM project.
    how_to_mesh : Callable
        Workflow callback handling meshing tasks.
    how_to_run : Callable
        Workflow callback handling solver execution tasks.
    how_to_clean : Callable = clean_case
        Workflow callback handling case cleanup tasks.
    get_args : str | Callable = ArgumentsPresets.common
        Callback for parsing command line arguments.
    """

    __slots__ = (
        "_root_dir",
        "_mesher",
        "_runner",
        "_cleaner",
        "_get_args",
    )

    def __init__(
            self,
            root_dir: Path,
            how_to_mesh: Callable,
            how_to_run: Callable,
            how_to_clean: Callable = clean_case,
            get_args: str | Callable = ArgumentsPresets.common,
        ) -> None:
        if not root_dir.exists():
            raise FileNotFoundError(root_dir)

        if isinstance(get_args, str):
            get_args = getattr(ArgumentsPresets, get_args)

        self._root_dir = root_dir
        self._mesher   = ensure_openfoam_case(how_to_mesh)
        self._runner   = ensure_openfoam_case(how_to_run)
        self._cleaner  = ensure_openfoam_case(how_to_clean)
        self._get_args = get_args

    def valid_options(self, args: Any) -> bool:
        """ Check if parsed arguments contain valid execution options.

        Parameters
        ----------
        args : Any
            Parsed command line arguments object.

        Returns
        -------
        bool
            True if at least one workflow execution option is set.
        """
        return args.mesh or args.run or args.clean or args.reconstruct

    def __call__(self, *args: Any, **kwargs: Any) -> None:
        """ Execute CLI workflow based on parsed command line arguments.

        Parameters
        ----------
        *args : Any
            Positional CLI arguments passed from entrypoint.
        **kwargs : Any
            Keyword CLI arguments passed from entrypoint.

        Returns
        -------
        None
            Executes target workflow phase directly.
        """
        os.chdir(self._root_dir)

        parser = self._get_args()
        args = parser.parse_args()

        process = args.mesh or args.run

        if not self.valid_options(args):
            parser.print_help()
            return

        if not process and args.clean:
            self._cleaner(self._root_dir)
            return

        ###
        # Partial cleans to enable process steps
        ###

        if args.clean_logs:
            clean_logs(self._root_dir)

        if args.clean_processors:
            clean_processors_dirs(self._root_dir)

        ###
        # Reconstruction / etc
        ###

        if args.reconstruct and not process:
            Runner.reconstruct(latest=args.latest)
            return

        ###
        # Main workflows
        ###

        if args.mesh:
            if args.clean or validate_input(
                "Do you want to clean the case before meshing?"
            ):
                self._cleaner(root_dir=self._root_dir)

            self._mesher(args, root_dir=self._root_dir)
            return

        if args.run:
            self._runner(args, root_dir=self._root_dir)
            return
