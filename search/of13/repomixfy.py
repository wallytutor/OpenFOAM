# -*- coding: utf-8 -*-

import logging
import sys

from pathlib import Path
from subprocess import run


logging.basicConfig(
    stream = sys.stdout,
    level  = logging.INFO,
    format = '%(asctime)s - %(levelname)s - %(message)s'
)


class RepoMixfy:
    __slots__ = (
        "_url",
        "_branch",
        "_repo_dir",
        "_output_dir",
        "_ignore_files",
        "_ignore_dirs",
        "_ignore_ext",
    )

    def __init__(
            self,
            url: str,
            branch: str = "main",
            repo_dir: str | Path | None = None,
            output_dir: str | Path | None = None,
            ignore_files: list[str] | None = None,
            ignore_dirs: list[str] | None = None,
            ignore_ext: list[str] | None = None
        ) -> None:
        self._url: str = url
        self._branch: str = branch
        self._repo_dir: Path = self._resolve_dir(repo_dir)
        self._output_dir: Path = self._resolve_dir(output_dir)

        self._init_clone()

    def _resolve_dir(self, tmp_dir) -> Path:
        if not tmp_dir:
            tmp_dir = repository_name(self._url)

        if isinstance(tmp_dir, str):
            tmp_dir = Path(tmp_dir)

        if not tmp_dir.is_absolute():
            tmp_dir = Path.cwd().joinpath(tmp_dir)

        logging.info(f"Repository at {tmp_dir}")
        return tmp_dir

    def _init_clone(self) -> None:
        """ Initialize clone of the repository. """
        if self._repo_dir.exists() and self._repo_dir.is_dir():
            logging.info(f"Repository already cloned at {self._repo_dir}")
            return

        run(
            [
                "git",
                "clone",
                "--branch",
                self._branch,
                self._url,
                self._repo_dir
            ],
            check=True,
            capture_output=True
        )


def repository_name(url: str) -> str:
    """ Return the repository name from the URL. """
    if not url.endswith(".git"):
        raise ValueError(f"URL must end with .git: {url}")

    if len(blocks := url.split("/")) < 2:
        raise ValueError(f"Invalid URL: {url}")

    return blocks[-1].split(".")[0]


#region: tests
def test_repository_name():
    url = "https://github.com/OpenFOAM/OpenFOAM-13.git"
    assert repository_name(url) == "OpenFOAM-13"


def tests():
    test_repository_name()
#endregion

if __name__ == "__main__":
    tests()

    ignore_files = [
        "Allwmake",
        "Allwclean",
        "Allmake",
        "Allclean",
        "Alltest",
    ]

    ignore_dirs = [
        ".git",
        "wmake",
        "Make",
        "Doxygen",
    ]

    ignore_ext = [
        ".pdf",
        ".png",
        ".jpg",
        ".jpeg",
        ".gif",
        ".svg",
    ]

    mix = RepoMixfy(
        url = "https://github.com/OpenFOAM/OpenFOAM-13.git",
        branch = "master"
    )

# file_path = Path('example.txt')
# file_size_bytes = file_path.stat().st_size
# print(file_size_bytes)