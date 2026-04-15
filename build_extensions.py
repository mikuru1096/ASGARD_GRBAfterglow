from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path


COMMON_FLAGS = "-Ofast -march=native -funroll-loops -ffast-math -fno-signed-zeros -fno-trapping-math"
OMP_FLAGS = f"-fopenmp {COMMON_FLAGS} -flto"
OPENMP_LIBS = ["-lgomp"]
BUILD_LOGIC_VERSION = 2


def _configure_windows_toolchain_env() -> None:
    if os.name != "nt":
        return

    mingw_bin = Path(os.environ.get("ASGARD_MINGW_BIN", r"C:\msys64\mingw64\bin"))
    if not mingw_bin.is_dir():
        return

    path_entries = os.environ.get("PATH", "").split(os.pathsep)
    mingw_str = str(mingw_bin)
    if mingw_str not in path_entries:
        os.environ["PATH"] = mingw_str + os.pathsep + os.environ.get("PATH", "")

    gcc = mingw_bin / "gcc.exe"
    gxx = mingw_bin / "g++.exe"
    gfortran = mingw_bin / "gfortran.exe"
    ar = mingw_bin / "gcc-ar.exe"
    if gcc.is_file():
        os.environ["CC"] = str(gcc)
    if gxx.is_file():
        os.environ["CXX"] = str(gxx)
    if gfortran.is_file():
        os.environ["FC"] = str(gfortran)
    if ar.is_file():
        os.environ["AR"] = str(ar)
    if hasattr(os, "add_dll_directory"):
        os.add_dll_directory(mingw_str)


def _clean_build_outputs(directory: Path) -> None:
    for pattern in ("*.so", "*.pyd", "*.o", "*.mod"):
        for path in directory.glob(pattern):
            path.unlink()


def _module_output_paths(directory: Path, module_name: str) -> list[Path]:
    outputs: list[Path] = []
    for pattern in (f"{module_name}*.pyd", f"{module_name}*.so"):
        outputs.extend(directory.glob(pattern))
    return outputs


def _hash_module_inputs(root: Path, cwd: Path, sources: list[str], fflags: str | None, extra_args: list[str] | None) -> str:
    digest = hashlib.sha256()
    digest.update(str(BUILD_LOGIC_VERSION).encode())
    digest.update(sys.version.encode())
    digest.update(os.name.encode())
    digest.update(str(fflags).encode())
    digest.update(json.dumps(extra_args or []).encode())
    for source in sources:
        source_path = (cwd / source).resolve()
        digest.update(str(source_path.relative_to(root)).encode())
        digest.update(source_path.read_bytes())
    return digest.hexdigest()


def _latest_timestamp(paths: list[Path]) -> float:
    return max(path.stat().st_mtime for path in paths)


def _is_timestamp_fresh(directory: Path, module_name: str, sources: list[str], script_path: Path) -> bool:
    outputs = _module_output_paths(directory, module_name)
    if not outputs:
        return False
    newest_output = _latest_timestamp(outputs)
    dependency_paths = [(directory / source).resolve() for source in sources]
    dependency_paths.append(script_path)
    newest_dependency = _latest_timestamp(dependency_paths)
    return newest_output >= newest_dependency


def _load_cache(cache_path: Path) -> dict[str, str]:
    if not cache_path.is_file():
        return {}
    try:
        data = json.loads(cache_path.read_text())
    except json.JSONDecodeError:
        return {}
    if not isinstance(data, dict):
        return {}
    return {str(key): str(value) for key, value in data.items()}


def _save_cache(cache_path: Path, cache: dict[str, str]) -> None:
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(cache, indent=2, sort_keys=True))


def _write_build_log(log_path: Path, command: list[str], cwd: Path, result: subprocess.CompletedProcess[str]) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        f"cwd: {cwd}",
        f"command: {' '.join(command)}",
        f"returncode: {result.returncode}",
        "",
        "--- stdout ---",
        result.stdout or "",
        "",
        "--- stderr ---",
        result.stderr or "",
    ]
    log_path.write_text("\n".join(lines), encoding="utf-8")


def _tail_output(text: str, max_lines: int = 120) -> str:
    lines = text.strip().splitlines()
    if len(lines) <= max_lines:
        return "\n".join(lines)
    return "\n".join(lines[-max_lines:])


def _build_module(
    module_name: str,
    cwd: Path,
    sources: list[str],
    log_dir: Path,
    verbose: bool,
    fflags: str | None = None,
    extra_args: list[str] | None = None,
) -> float:
    env = os.environ.copy()
    py_dir = Path(sys.executable).resolve().parent
    py_scripts = py_dir / "Scripts"
    extra_path_entries = [str(py_dir)]
    if py_scripts.is_dir():
        extra_path_entries.append(str(py_scripts))
    if os.name == "nt":
        mingw_bin = env.get("ASGARD_MINGW_BIN", r"C:\msys64\mingw64\bin")
        mingw_dir = Path(mingw_bin)
        if mingw_dir.is_dir():
            extra_path_entries.insert(0, mingw_bin)
            gcc = mingw_dir / "gcc.exe"
            gxx = mingw_dir / "g++.exe"
            gfortran = mingw_dir / "gfortran.exe"
            ar = mingw_dir / "gcc-ar.exe"
            if gcc.is_file():
                env["CC"] = str(gcc)
            if gxx.is_file():
                env["CXX"] = str(gxx)
            if gfortran.is_file():
                env["FC"] = str(gfortran)
            if ar.is_file():
                env["AR"] = str(ar)
    env["PATH"] = os.pathsep.join(extra_path_entries) + os.pathsep + env["PATH"]
    if fflags is not None:
        env["FFLAGS"] = fflags
        env["F90FLAGS"] = fflags

    command = [sys.executable, "-m", "numpy.f2py", "-m", module_name, "-c", *sources]
    if extra_args:
        command.extend(extra_args)
    start = time.perf_counter()
    if verbose:
        subprocess.run(command, cwd=cwd, check=True, env=env)
        return time.perf_counter() - start

    result = subprocess.run(
        command,
        cwd=cwd,
        env=env,
        check=False,
        capture_output=True,
        text=True,
        errors="replace",
    )
    elapsed = time.perf_counter() - start
    log_path = log_dir / f"{module_name}.log"
    _write_build_log(log_path, command, cwd, result)
    if result.returncode != 0:
        output = "\n".join(part for part in (result.stdout, result.stderr) if part)
        tail = _tail_output(output)
        print(f"Build failed for {module_name}. Full log: {log_path}")
        if tail:
            print(tail)
        raise subprocess.CalledProcessError(result.returncode, command)
    return elapsed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true", help="rebuild all selected modules")
    parser.add_argument("--clean", action="store_true", help="remove built artifacts and cache before building")
    parser.add_argument("--module", action="append", dest="modules", help="only build the named module; can be repeated")
    parser.add_argument("--verbose", action="store_true", help="stream raw f2py/meson output")
    args = parser.parse_args()

    _configure_windows_toolchain_env()
    root = Path(__file__).resolve().parent
    cache_path = root / ".buildcache" / "extensions.json"
    log_dir = root / ".buildcache" / "logs"
    src = root / "src"
    dyn = src / "Dynamics"
    ele = src / "Electron"
    itp = src / "Interpolation"
    rad = src / "Radiation"
    modules = [
        ("Constants", src, ["Constants.f90"], None, None),
        ("Dynamics_reverse", dyn, ["../Constants.f90", "Dynamics_reverse.f90"], COMMON_FLAGS, None),
        ("Dynamics_forward", dyn, ["../Constants.f90", "dynamics_common.f90", "Dynamics_forward.f90"], COMMON_FLAGS, None),
        ("FS_electron_weno5", ele, ["../Constants.f90", "../utils/adaptive_2nd_resampling_mod.f90", "electron_common.f90", "calling_modules.f90", "FS_electron_weno5.f90"], OMP_FLAGS, OPENMP_LIBS),
        ("FS_electron_slc1", ele, ["../Constants.f90", "../utils/adaptive_2nd_resampling_mod.f90", "electron_common.f90", "calling_modules.f90", "FS_electron_slc1.f90"], OMP_FLAGS, OPENMP_LIBS),
        ("FS_electron_fullhide", ele, ["../Constants.f90", "../utils/adaptive_2nd_resampling_mod.f90", "electron_common.f90", "calling_modules.f90", "FS_electron_fullhide.f90"], OMP_FLAGS, OPENMP_LIBS),
        ("FS_electron_t2g1", ele, ["../Constants.f90", "../utils/adaptive_2nd_resampling_mod.f90", "electron_common.f90", "calling_modules.f90", "FS_electron_t2g1.f90"], OMP_FLAGS, OPENMP_LIBS),
        ("SED_interpolation", itp, ["../Constants.f90", "interpolation_common.f90", "SED_interpolation.f90"], OMP_FLAGS, OPENMP_LIBS),
        ("SED_interpolation_structured", itp, ["../Constants.f90", "interpolation_common.f90", "SED_interpolation_structured.f90"], OMP_FLAGS, OPENMP_LIBS),
        ("Annihilation", rad, ["../Constants.f90", "radiation_common.f90", "Annihilation.f90"], OMP_FLAGS, OPENMP_LIBS),
        ("Seed_reverse", rad, ["../Constants.f90", "radiation_common.f90", "Seed_reverse.f90"], OMP_FLAGS, OPENMP_LIBS),
        ("SSC_spec", rad, ["../Constants.f90", "radiation_common.f90", "SSC_spec.f90"], OMP_FLAGS, OPENMP_LIBS),
    ]
    selected = set(args.modules or [])
    if selected:
        modules = [spec for spec in modules if spec[0] in selected]
        missing = selected.difference({spec[0] for spec in modules})
        if missing:
            raise SystemExit(f"Unknown module(s): {', '.join(sorted(missing))}")
    if args.clean:
        for directory in (src, dyn, ele, itp, rad):
            _clean_build_outputs(directory)
        if cache_path.exists():
            cache_path.unlink()

    cache = _load_cache(cache_path)

    print("Compile start")
    for module_name, cwd, sources, fflags, extra_args in modules:
        signature = _hash_module_inputs(root, cwd, sources, fflags, extra_args)
        if not args.force:
            cached_signature = cache.get(module_name)
            if cached_signature == signature and _module_output_paths(cwd, module_name):
                print(f"Skip {module_name}: unchanged")
                continue
            if cached_signature is None and _is_timestamp_fresh(cwd, module_name, sources, Path(__file__).resolve()):
                cache[module_name] = signature
                print(f"Skip {module_name}: outputs are newer than sources")
                continue
        print(f"Build {module_name}")
        elapsed = _build_module(
            module_name,
            cwd,
            sources,
            log_dir,
            args.verbose,
            fflags,
            extra_args,
        )
        cache[module_name] = signature
        print(f"Done {module_name}: {elapsed:.2f}s")
    _save_cache(cache_path, cache)
    print("Compile complete!")


if __name__ == "__main__":
    main()
