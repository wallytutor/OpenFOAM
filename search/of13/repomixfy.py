# -*- coding: utf-8 -*-

import fnmatch
import logging
import os
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
        self._repo_dir: Path = self._resolve_repo_dir(repo_dir)
        self._output_dir: Path = self._resolve_output_dir(output_dir)
        self._ignore_files: list[str] = ignore_files or []
        self._ignore_dirs: list[str] = ignore_dirs or []
        self._ignore_ext: list[str] = ignore_ext or []

        self._init_clone()
        self._init_files()

    def _resolve_repo_dir(self, repo_dir) -> Path:
        if not repo_dir:
            repo_dir = repository_name(self._url)

        return self._resolve_dir(repo_dir)

    def _resolve_output_dir(self, output_dir) -> Path:
        if not output_dir:
            output_dir = self._repo_dir.name + "_mix"

        return self._resolve_dir(output_dir)

    def _resolve_dir(self, tmp_dir) -> Path:
        if isinstance(tmp_dir, str):
            tmp_dir = Path(tmp_dir)

        if not tmp_dir.is_absolute():
            tmp_dir = Path.cwd().joinpath(tmp_dir)

        logging.info(f"Repository at {tmp_dir}")
        return tmp_dir

    def _is_dir_ignored(self, dir_parts: tuple[str, ...]) -> bool:
        """ Check if directory parts match any rule in _ignore_dirs. """
        if not dir_parts or not self._ignore_dirs:
            return False

        for raw_pattern in self._ignore_dirs:
            pattern = raw_pattern.strip("/").replace("\\", "/")

            if not pattern:
                continue

            pattern_parts = tuple(pattern.split("/"))

            if len(pattern_parts) == 1:
                pat = pattern_parts[0]

                for part in dir_parts:
                    if part == pat or fnmatch.fnmatch(part, pat):
                        return True

            else:
                n = len(pattern_parts)

                for i in range(len(dir_parts) - n + 1):
                    slice_parts = dir_parts[i:i + n]

                    if all(
                        s == p or fnmatch.fnmatch(s, p)
                        for s, p in zip(slice_parts, pattern_parts)
                    ):
                        return True

        return False

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

    def _init_files(self) -> None:
        """ Initialize .repomixfy file with repository files list. """
        if not self._repo_dir.exists():
            return

        self._output_dir.mkdir(parents=True, exist_ok=True)

        repomixfy_path = self._output_dir / ".repomixfy"

        ignore_ext = {
            ext if ext.startswith(".") else f".{ext}"
            for ext in self._ignore_ext
        }

        file_paths: list[str] = []

        for root, dirs, files in os.walk(self._repo_dir):
            root_path = Path(root)
            rel_root = root_path.relative_to(self._repo_dir)

            dirs[:] = [
                d for d in dirs
                if not self._is_dir_ignored((rel_root / d).parts)
            ]

            for file in files:
                # if file == ".repomixfy":
                #     continue

                file_path = root_path / file
                rel_path = file_path.relative_to(self._repo_dir)

                if self._is_dir_ignored(rel_path.parts[:-1]):
                    continue

                if any(
                    file == pat or fnmatch.fnmatch(file, pat)
                    for pat in self._ignore_files
                ):
                    continue

                if file_path.suffix in ignore_ext:
                    continue

                file_paths.append(rel_path.as_posix())

        file_paths.sort()

        with repomixfy_path.open("w", encoding="utf-8") as f:
            for rel_path in file_paths:
                f.write(f"{rel_path}\n")

        logging.info(
            f"Created .repomixfy with {len(file_paths)} files at "
            f"{repomixfy_path}"
        )


def repository_name(url: str) -> str:
    """ Return the repository name from the URL. """
    if not url.endswith(".git"):
        raise ValueError(f"URL must end with .git: {url}")

    if len(blocks := url.split("/")) < 2:
        raise ValueError(f"Invalid URL: {url}")

    return blocks[-1].split(".")[0]


def is_text_file(file_path: str | Path, chunk_size: int = 1024) -> bool:
    """ Return True if the file is a text file. """
    with open(file_path, "rb") as f:
        chunk = f.read(chunk_size)

    if b"\x00" in chunk:
        return False

    try:
        chunk.decode("utf-8")
        return True
    except UnicodeDecodeError:
        return False


#region: tests
def test_repository_name():
    url = "https://github.com/OpenFOAM/OpenFOAM-13.git"
    assert repository_name(url) == "OpenFOAM-13"


def test_init_files():
    import tempfile

    with tempfile.TemporaryDirectory() as tmpdir:
        repo_dir = Path(tmpdir) / "test_repo"
        repo_dir.mkdir()

        (repo_dir / "file1.txt").write_text("hello", encoding="utf-8")
        (repo_dir / "Allwmake").write_text("clean", encoding="utf-8")
        (repo_dir / "image.png").write_text("img", encoding="utf-8")

        git_dir = repo_dir / ".git"
        git_dir.mkdir()
        (git_dir / "config").write_text("git config", encoding="utf-8")

        sub_dir = repo_dir / "src"
        sub_dir.mkdir()
        (sub_dir / "main.py").write_text("print('hi')", encoding="utf-8")

        # Test deep nested directory matching
        make_dir = sub_dir / "OpenFOAM" / "Make"
        make_dir.mkdir(parents=True)
        (make_dir / "options").write_text("options", encoding="utf-8")

        # Test multi-level directory path matching
        rules_dir = repo_dir / "wmake" / "rules"
        rules_dir.mkdir(parents=True)
        (rules_dir / "linux").write_text("rule", encoding="utf-8")

        # Test glob pattern directory matching
        build_dir = repo_dir / "build_temp"
        build_dir.mkdir()
        (build_dir / "temp.o").write_text("obj", encoding="utf-8")

        ignore_files = ["Allwmake"]
        ignore_dirs = [".git", "Make", "wmake/rules", "build_*"]
        ignore_ext = [".png"]

        mix = RepoMixfy.__new__(RepoMixfy)
        mix._repo_dir = repo_dir
        mix._ignore_files = ignore_files
        mix._ignore_dirs = ignore_dirs
        mix._ignore_ext = ignore_ext

        mix._init_files()

        repomixfy_file = repo_dir / ".repomixfy"
        assert repomixfy_file.exists()

        lines = repomixfy_file.read_text(encoding="utf-8").splitlines()
        assert lines == ["file1.txt", "src/main.py"]


def tests():
    test_repository_name()
    test_init_files()
#endregion

if __name__ == "__main__":
    # tests()

    ignore_files = [
        ".gitattributes",
        ".gitignore",
        "COPYING",
        "Allwmake",
        "Allwclean",
        "Allmake",
        "Allclean",
        "Allrun",
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
        branch = "master",
        ignore_files = ignore_files,
        ignore_dirs = ignore_dirs,
        ignore_ext = ignore_ext
    )

# file_path = Path('example.txt')
# file_size_bytes = file_path.stat().st_size
# print(file_size_bytes)