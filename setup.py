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


def _utf8_env() -> dict[str, str]:
    env = os.environ.copy()
    env["PYTHONUTF8"] = "1"
    env["PYTHONIOENCODING"] = "utf-8"
    env["LANG"] = "C.UTF-8"
    env["LC_ALL"] = "C.UTF-8"
    return env


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
        subprocess.run(command, cwd=ROOT, check=True, env=_utf8_env(), text=True, encoding="utf-8")


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
            "electron_forward_weno5_1d*.so",
            "electron_forward_weno5_1d*.pyd",
            "electron_forward_slc1_1d*.so",
            "electron_forward_slc1_1d*.pyd",
            "electron_forward_charint_1d*.so",
            "electron_forward_charint_1d*.pyd",
            "electron_forward_charint_2d*.so",
            "electron_forward_charint_2d*.pyd",
            "electron_forward_fullhide_1d*.so",
            "electron_forward_fullhide_1d*.pyd",
            "electron_forward_transport_2d*.so",
            "electron_forward_transport_2d*.pyd",
            "electron_forward_t2g1_1d*.so",
            "electron_forward_t2g1_1d*.pyd",
            "electron_radiation*.so",
            "electron_radiation*.pyd",
            "electron_reverse_kernel*.so",
            "electron_reverse_kernel*.pyd",
            "*.f90",
        ],
        "src.Interpolation": ["SED_interpolation*.so", "SED_interpolation*.pyd", "SED_interpolation_structured*.so", "SED_interpolation_structured*.pyd", "*.f90"],
        "src.Radiation": ["radiation_gamma_gamma_absorption*.so", "radiation_gamma_gamma_absorption*.pyd", "radiation_reverse_seed*.so", "radiation_reverse_seed*.pyd", "radiation_ssc_spectrum*.so", "radiation_ssc_spectrum*.pyd", "*.f90", "*.py"],
        "src.Structured": ["structured_jet_1d*.so", "structured_jet_1d*.pyd", "*.f90", "*.py"],
    },
)
