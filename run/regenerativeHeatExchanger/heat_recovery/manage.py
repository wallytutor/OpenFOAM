# -*- coding: utf-8 -*-
from pathlib import Path
import shutil
import subprocess


class CaseManager:
    @staticmethod
    def run(log, cmd, blocking=True, force=False):
        """ Runs a command and logs the output to a file.

        If the log file already exists, it will not run the command
        unless force is True. By default the command is run in blocking
        mode, which means that the function will wait for the command to
        finish before returning.
        """
        if Path(log).exists() and not force:
            print(f"Already run, check log {log} for details")
            return

        with open(log, "w") as f:
            print(f"Logging to {log}: {' '.join(cmd)}")
            if blocking:
                subprocess.run(cmd, stdout=f, stderr=f, check=True)
            else:
                subprocess.Popen(cmd, stdout=f, stderr=f)

    @staticmethod
    def handle_field(region, field, init_dir="0.000000e+00", copy=True):
        """ Copy or expand a field from the original case to the runtime. """
        log = f"{init_dir}/{region}/{field}"
        src = f"0.orig/{region}/{field}"

        if copy:
            shutil.copy(src, log)
            return

        CaseManager.run(log, ["foamDictionary", src, "-expand"])