from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def _detect_platform() -> str:
    if os.name == "nt":
        return "windows"
    if sys.platform.startswith("linux"):
        return "ubuntu" if Path("/proc/version").exists() else "linux"
    return sys.platform


def _require_fortran_compiler(platform_name: str) -> None:
    if platform_name == "windows":
        mingw_bin = Path(os.environ.get("ASGARD_MINGW_BIN", r"C:\msys64\mingw64\bin"))
        candidates = [
            shutil.which("gfortran"),
            str(mingw_bin / "gfortran.exe") if mingw_bin.is_dir() else None,
        ]
    else:
        candidates = [shutil.which("gfortran"), "/usr/bin/gfortran"]
    for candidate in candidates:
        if candidate and Path(candidate).exists():
            return
    raise SystemExit("Missing gfortran toolchain. Install GNU Fortran first, then rerun install.py.")


def _run(command: list[str]) -> None:
    subprocess.run(
        command,
        cwd=ROOT,
        check=True,
        env=dict(
            os.environ,
            PYTHONUTF8="1",
            PYTHONIOENCODING="utf-8",
            LANG="C.UTF-8",
            LC_ALL="C.UTF-8",
        ),
        text=True,
        encoding="utf-8",
    )


def _require_uv() -> str:
    for candidate in (shutil.which("uv"), str(Path.home() / ".local" / "bin" / "uv")):
        if candidate and Path(candidate).exists():
            return candidate
    _run([sys.executable, "-m", "pip", "install", "--user", "uv"])
    for candidate in (shutil.which("uv"), str(Path.home() / ".local" / "bin" / "uv")):
        if candidate and Path(candidate).exists():
            return candidate
    raise SystemExit("uv installation failed. Install uv and rerun install.py.")


def main() -> None:
    platform_name = _detect_platform()
    if platform_name not in {"windows", "linux", "ubuntu"}:
        raise SystemExit(f"Unsupported platform: {platform_name}")

    _require_fortran_compiler(platform_name)
    uv = _require_uv()

    py = ROOT / ".venv" / ("Scripts/python.exe" if platform_name == "windows" else "bin/python")

    if not py.exists():
        _run([uv, "venv", "--python", sys.executable, ".venv"])
    _run([uv, "pip", "install", "--python", str(py), "--upgrade", "pip", "setuptools", "wheel", "numpy", "meson", "ninja"])
    install_command = [uv, "pip", "install", "--python", str(py), "--upgrade", "--no-build-isolation"]
    install_command.append(".")
    os.environ["ASGARD_SKIP_NATIVE_BUILD"] = "1"
    _run(install_command)
    _run([str(py), "build_extensions.py", "--force"])


if __name__ == "__main__":
    main()
