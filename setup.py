from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

from setuptools import Command, setup
from setuptools.command.build_py import build_py as _build_py
from setuptools.command.develop import develop as _develop
from setuptools.command.sdist import sdist as _sdist


ROOT = Path(__file__).resolve().parent


def _skip_native_build() -> bool:
    return os.environ.get("ASGARD_SKIP_NATIVE_BUILD") == "1"


class build_native(Command):
    description = "build ASGARD native Fortran extensions in-place"
    user_options: list[tuple[str, str | None, str]] = []

    def initialize_options(self) -> None:
        self.force = None

    def finalize_options(self) -> None:
        self.force = False

    def run(self) -> None:
        if _skip_native_build():
            return
        command = [sys.executable, "build_extensions.py"]
        subprocess.run(command, cwd=ROOT, check=True)


class build_py(_build_py):
    def run(self) -> None:
        if not _skip_native_build():
            self.run_command("build_native")
        super().run()


class develop(_develop):
    def run(self) -> None:
        if not _skip_native_build():
            self.run_command("build_native")
        super().run()


class sdist(_sdist):
    def run(self) -> None:
        super().run()


setup(
    py_modules=["lc_spec_demo", "build_extensions"],
    cmdclass={
        "build_native": build_native,
        "build_py": build_py,
        "develop": develop,
        "sdist": sdist,
    },
    package_data={
        "asgard_core": ["data_light_curve/*", "data_spectrum/*"],
        "src": ["Constants*.so", "Constants*.pyd", "*.f90"],
        "src.Dynamics": ["Dynamics_forward*.so", "Dynamics_forward*.pyd", "Dynamics_reverse*.so", "Dynamics_reverse*.pyd", "*.f90"],
        "src.Electron": [
            "FS_electron_weno5_1d*.so",
            "FS_electron_weno5_1d*.pyd",
            "FS_electron_slc1_1d*.so",
            "FS_electron_slc1_1d*.pyd",
            "FS_electron_charint_1d*.so",
            "FS_electron_charint_1d*.pyd",
            "FS_electron_charint_2d*.so",
            "FS_electron_charint_2d*.pyd",
            "FS_electron_fullhide_1d*.so",
            "FS_electron_fullhide_1d*.pyd",
            "FS_electron_fullhide_2d*.so",
            "FS_electron_fullhide_2d*.pyd",
            "FS_electron_t2g1_1d*.so",
            "FS_electron_t2g1_1d*.pyd",
            "electron_radiation*.so",
            "electron_radiation*.pyd",
            "electron_reverse_kernel*.so",
            "electron_reverse_kernel*.pyd",
            "*.f90",
        ],
        "src.Interpolation": ["SED_interpolation*.so", "SED_interpolation*.pyd", "SED_interpolation_structured*.so", "SED_interpolation_structured*.pyd", "*.f90"],
        "src.Radiation": ["Annihilation*.so", "Annihilation*.pyd", "Seed_reverse*.so", "Seed_reverse*.pyd", "SSC_spec*.so", "SSC_spec*.pyd", "*.f90", "*.py"],
    },
)
