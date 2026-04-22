from __future__ import annotations

import argparse
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
    subprocess.run(command, cwd=ROOT, check=True)


def _venv_python(platform_name: str) -> Path:
    if platform_name == "windows":
        return ROOT / ".venv" / "Scripts" / "python.exe"
    return ROOT / ".venv" / "bin" / "python"


def _find_uv() -> str | None:
    candidates = [
        shutil.which("uv"),
        str(Path.home() / ".local" / "bin" / "uv"),
    ]
    for candidate in candidates:
        if candidate and Path(candidate).exists():
            return candidate
    return None


def _require_uv() -> str:
    uv = _find_uv()
    if uv is not None:
        return uv
    _run([sys.executable, "-m", "pip", "install", "--user", "uv"])
    uv = _find_uv()
    if uv is None:
        raise SystemExit("uv installation failed. Install uv and rerun install.py.")
    return uv


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--editable", action="store_true", help="install in editable mode")
    parser.add_argument("--with-fit", action="store_true", help="install fitting extras")
    parser.add_argument("--with-compare", action="store_true", help="install comparison/benchmark dependencies")
    parser.add_argument("--skip-build", action="store_true", help="install Python package without rebuilding Fortran extensions")
    args = parser.parse_args()

    platform_name = _detect_platform()
    if platform_name not in {"windows", "linux", "ubuntu"}:
        raise SystemExit(f"Unsupported platform: {platform_name}")

    _require_fortran_compiler(platform_name)
    uv = _require_uv()

    extras_requested = []
    if args.with_fit:
        extras_requested.append("fit")
    if args.with_compare:
        extras_requested.append("compare")
    extras = f"[{','.join(extras_requested)}]" if extras_requested else ""
    target = f".{extras}"
    py = _venv_python(platform_name)

    if not py.exists():
        _run([uv, "venv", "--python", sys.executable, ".venv"])
    _run([uv, "pip", "install", "--python", str(py), "--upgrade", "pip", "setuptools", "wheel", "numpy", "meson", "ninja"])
    install_command = [uv, "pip", "install", "--python", str(py), "--upgrade", "--no-build-isolation"]
    if args.editable:
        install_command.extend(["-e", target])
    else:
        install_command.append(target)
    if args.skip_build:
        os.environ["ASGARD_SKIP_NATIVE_BUILD"] = "1"
    else:
        os.environ.pop("ASGARD_SKIP_NATIVE_BUILD", None)
    _run(install_command)
    if not args.skip_build:
        _run([str(py), "build_extensions.py", "--force"])


if __name__ == "__main__":
    main()
