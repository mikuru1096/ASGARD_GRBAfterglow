from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

from setuptools import Command, setup
from setuptools.command.build_py import build_py as _build_py
from setuptools.command.develop import develop as _develop


ROOT = Path(__file__).resolve().parent


def _skip_native_build() -> bool:
    return os.environ.get("ASGARD_SKIP_NATIVE_BUILD") == "1"


def _native_patterns(*stems: str) -> list[str]:
    return [pattern for stem in stems for pattern in (f"{stem}*.so", f"{stem}*.pyd")]


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
        subprocess.run(
            [sys.executable, "build_extensions.py"],
            cwd=ROOT,
            check=True,
            env=dict(os.environ, PYTHONUTF8="1", PYTHONIOENCODING="utf-8", LANG="C.UTF-8", LC_ALL="C.UTF-8"),
            text=True,
            encoding="utf-8",
        )


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


setup(
    py_modules=["lc_spec_demo", "build_extensions"],
    cmdclass={
        "build_native": build_native,
        "build_py": build_py,
        "develop": develop,
    },
    package_data={
        "asgard_core": ["data_light_curve/*", "data_spectrum/*"],
        "src": [*_native_patterns("Constants"), "*.f90"],
        "src.Dynamics": [*_native_patterns("Dynamics_forward", "Dynamics_reverse"), "*.f90"],
        "src.Electron": [
            *_native_patterns(
                "electron_forward_weno5_1d",
                "electron_forward_slc1_1d",
                "electron_forward_charint_1d",
                "electron_forward_charint_2d",
                "electron_forward_fullhide_1d",
                "electron_forward_transport_2d",
                "electron_forward_t2g1_1d",
                "electron_radiation",
                "electron_reverse_kernel",
            ),
            "*.f90",
        ],
        "src.Interpolation": [
            *_native_patterns("SED_interpolation", "SED_interpolation_structured"),
            "*.f90",
        ],
        "src.Radiation": [
            *_native_patterns(
                "radiation_gamma_gamma_absorption",
                "radiation_reverse_seed",
                "radiation_ssc_spectrum",
            ),
            "*.f90",
            "*.py",
        ],
        "src.Structured": [*_native_patterns("structured_jet_1d"), "*.f90", "*.py"],
    },
)
