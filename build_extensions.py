from __future__ import annotations

import argparse
import os
import shlex
import shutil
import subprocess
import sys
import sysconfig
import time
from pathlib import Path


COMMON_FLAGS = "-Ofast -march=native -funroll-loops -ffast-math -fno-signed-zeros -fno-trapping-math -ffree-line-length-none"
OMP_FLAGS = f"-fopenmp {COMMON_FLAGS} -flto"
OPENMP_LIBS = ["-lgomp"]


def _detect_build_platform() -> str:
    if os.name == "nt":
        return "windows-mingw"
    if os.name == "posix":
        return "linux-gfortran"
    return f"unsupported:{os.name}"


def _prepare_build_env() -> dict[str, str]:
    env = os.environ.copy()
    py_dir = Path(sys.executable).resolve().parent
    extra_path_entries = [str(py_dir)]
    py_scripts = py_dir / "Scripts"
    if py_scripts.is_dir():
        extra_path_entries.append(str(py_scripts))
    if os.name == "nt":
        mingw_bin = env.get("ASGARD_MINGW_BIN", r"C:\msys64\mingw64\bin")
        mingw_dir = Path(mingw_bin)
        if mingw_dir.is_dir():
            extra_path_entries.insert(0, mingw_bin)
            compiler_map = {
                "CC": mingw_dir / "gcc.exe",
                "CXX": mingw_dir / "g++.exe",
                "FC": mingw_dir / "gfortran.exe",
                "AR": mingw_dir / "gcc-ar.exe",
            }
            for name, path in compiler_map.items():
                if path.is_file():
                    env[name] = str(path)
            if hasattr(os, "add_dll_directory"):
                os.add_dll_directory(mingw_bin)
    env["PATH"] = os.pathsep.join(extra_path_entries) + os.pathsep + env["PATH"]
    return env


def _clean_build_outputs(directory: Path) -> None:
    for pattern in ("*.so", "*.pyd", "*.o", "*.mod"):
        for path in directory.glob(pattern):
            path.unlink()


def _module_output_paths(directory: Path, module_name: str) -> list[Path]:
    outputs: list[Path] = []
    for pattern in (f"{module_name}*.pyd", f"{module_name}*.so"):
        outputs.extend(directory.glob(pattern))
    return outputs

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


def _run_command(
    command: list[str],
    cwd: Path,
    env: dict[str, str],
    log_path: Path,
    verbose: bool,
) -> subprocess.CompletedProcess[str]:
    if verbose:
        return subprocess.run(command, cwd=cwd, check=True, env=env, text=True)
    result = subprocess.run(
        command,
        cwd=cwd,
        env=env,
        check=False,
        capture_output=True,
        text=True,
        errors="replace",
    )
    _write_build_log(log_path, command, cwd, result)
    if result.returncode != 0:
        output = "\n".join(part for part in (result.stdout, result.stderr) if part)
        tail = _tail_output(output)
        print(f"Build failed for {' '.join(command[:3])}. Full log: {log_path}")
        if tail:
            print(tail)
        raise subprocess.CalledProcessError(result.returncode, command)
    return result


def _build_fs_electron_ordered_fallback(
    root: Path,
    module_name: str,
    cwd: Path,
    sources: list[str],
    log_dir: Path,
    verbose: bool,
    fflags: str | None,
    extra_args: list[str] | None,
) -> float:
    env = _prepare_build_env()
    fflags_override = env.get("ASGARD_FFLAGS_OVERRIDE")
    if fflags_override:
        fflags = fflags_override
    if fflags is not None:
        env["FFLAGS"] = fflags
        env["F90FLAGS"] = fflags

    fc = env.get("FC") or shutil.which("gfortran", path=env["PATH"])
    if not fc:
        raise RuntimeError(f"{module_name} fallback build requires gfortran in PATH or FC.")

    build_dir = root / ".buildcache" / "fs_electron_ordered_fallback" / module_name
    build_dir.mkdir(parents=True, exist_ok=True)
    for path in build_dir.glob("*"):
        if path.is_file():
            path.unlink()

    compile_flags = shlex.split(fflags or "")
    compile_flags.append("-fPIC")
    compile_flags.extend(["-I", str(build_dir), "-J", str(build_dir)])
    object_paths: list[Path] = []
    start = time.perf_counter()
    for source in sources:
        source_path = (cwd / source).resolve()
        object_path = build_dir / f"{source_path.stem}.o"
        command = [fc, "-c", *compile_flags, str(source_path), "-o", str(object_path)]
        _run_command(command, cwd, env, log_dir / f"{module_name}_fallback_compile_{source_path.stem}.log", verbose)
        object_paths.append(object_path)

    pyf_path = build_dir / f"{module_name}.pyf"
    main_source_name = Path(sources[-1]).name
    source_rel = os.path.relpath((cwd / main_source_name).resolve(), build_dir)
    entry_names = [module_name.lower()]
    if module_name == "FS_electron_charint_1d":
        entry_names.append("fs_electron_charint_1d_affine_step_test")
    if module_name == "electron_get_y":
        entry_names = ["get_nu_a", "get_syn_selected"]
    if module_name == "SSC_spec":
        entry_names = ["ssc_spec", "ssc_spec_nonuniform"]
    signature_command = [
        sys.executable,
        "-m",
        "numpy.f2py",
        "-m",
        module_name,
        "-h",
        pyf_path.name,
        "--overwrite-signature",
        source_rel,
        "only:",
        *entry_names,
        ":",
    ]
    _run_command(signature_command, build_dir, env, log_dir / f"{module_name}_fallback_signature.log", verbose)

    wrapper_command = [
        sys.executable,
        "-m",
        "numpy.f2py",
        pyf_path.name,
    ]
    _run_command(wrapper_command, build_dir, env, log_dir / f"{module_name}_fallback_wrapper.log", verbose)

    import numpy as np

    cc = env.get("CC") or shutil.which("cc", path=env["PATH"])
    if not cc:
        raise RuntimeError(f"{module_name} fallback build requires a C compiler in PATH or CC.")
    py_include = sysconfig.get_paths()["include"]
    py_platinclude = sysconfig.get_paths().get("platinclude", py_include)
    f2py_src = Path(np.__file__).resolve().parent / "f2py" / "src"
    ext_suffix = sysconfig.get_config_var("EXT_SUFFIX") or ".so"

    c_objects: list[Path] = []
    for source in (build_dir / f"{module_name}module.c", f2py_src / "fortranobject.c"):
        object_path = build_dir / f"{source.stem}.o"
        compile_c = [
            cc,
            "-c",
            "-fPIC",
            f"-I{py_include}",
            f"-I{py_platinclude}",
            f"-I{np.get_include()}",
            f"-I{f2py_src}",
            str(source),
            "-o",
            str(object_path),
        ]
        _run_command(compile_c, build_dir, env, log_dir / f"{module_name}_fallback_compile_{source.stem}.log", verbose)
        c_objects.append(object_path)

    wrapper_source = build_dir / f"{module_name}-f2pywrappers.f"
    wrapper_object = build_dir / f"{module_name}-f2pywrappers.o"
    if wrapper_source.is_file():
        compile_wrapper = [fc, "-c", "-fPIC", str(wrapper_source), "-o", str(wrapper_object)]
        _run_command(compile_wrapper, build_dir, env, log_dir / f"{module_name}_fallback_compile_wrappers.log", verbose)
        object_paths.append(wrapper_object)

    output_path = build_dir / f"{module_name}{ext_suffix}"
    link_command = [
        fc,
        "-shared",
        *(str(path) for path in c_objects),
        *(str(path) for path in object_paths),
        "-o",
        str(output_path),
    ]
    if extra_args:
        link_command.extend(extra_args)
    _run_command(link_command, build_dir, env, log_dir / f"{module_name}_fallback_link.log", verbose)
    for built_path in _module_output_paths(build_dir, module_name):
        target_path = cwd / built_path.name
        if target_path.exists():
            target_path.unlink()
        built_path.replace(target_path)
    return time.perf_counter() - start


def _build_module(
    module_name: str,
    cwd: Path,
    sources: list[str],
    log_dir: Path,
    verbose: bool,
    fflags: str | None = None,
    extra_args: list[str] | None = None,
) -> float:
    env = _prepare_build_env()
    fflags_override = env.get("ASGARD_FFLAGS_OVERRIDE")
    if fflags_override:
        fflags = fflags_override
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
    parser.add_argument("--clean", action="store_true", help="remove built artifacts before building")
    parser.add_argument("--module", action="append", dest="modules", help="only build the named module; can be repeated")
    parser.add_argument("--verbose", action="store_true", help="stream raw f2py/meson output")
    args = parser.parse_args()

    root = Path(__file__).resolve().parent
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
        ("FS_electron_weno5_1d", ele, ["../Constants.f90", "../Radiation/radiation_common.f90", "adaptive_resampling_mod.f90", "electron_common.f90", "electron_radiation_kernel.f90", "electron_cooling_kernel.f90", "calling_modules.f90", "FS_electron_weno5_1d.f90"], OMP_FLAGS, OPENMP_LIBS),
        ("FS_electron_slc1_1d", ele, ["../Constants.f90", "../Radiation/radiation_common.f90", "adaptive_resampling_mod.f90", "electron_common.f90", "electron_radiation_kernel.f90", "electron_cooling_kernel.f90", "calling_modules.f90", "FS_electron_slc1_1d.f90"], OMP_FLAGS, OPENMP_LIBS),
        ("FS_electron_charint_1d", ele, ["../Constants.f90", "../Radiation/radiation_common.f90", "adaptive_resampling_mod.f90", "electron_common.f90", "electron_radiation_kernel.f90", "electron_cooling_kernel.f90", "calling_modules.f90", "FS_electron_charint_1d.f90"], OMP_FLAGS, OPENMP_LIBS),
        ("FS_electron_fullhide_1d", ele, ["../Constants.f90", "../Radiation/radiation_common.f90", "adaptive_resampling_mod.f90", "electron_common.f90", "electron_radiation_kernel.f90", "electron_cooling_kernel.f90", "electron_seed_history_kernel.f90", "calling_modules.f90", "FS_electron_fullhide_1d.f90"], OMP_FLAGS, OPENMP_LIBS),
        ("FS_electron_fullhide_2d", ele, ["../Constants.f90", "../Radiation/radiation_common.f90", "adaptive_resampling_mod.f90", "electron_common.f90", "electron_radiation_kernel.f90", "electron_cooling_kernel.f90", "electron_seed_history_kernel.f90", "electron_transport_2d_kernel.f90", "calling_modules.f90", "FS_electron_fullhide_2d.f90"], OMP_FLAGS, OPENMP_LIBS),
        ("FS_electron_charint_2d", ele, ["../Constants.f90", "../Radiation/radiation_common.f90", "adaptive_resampling_mod.f90", "electron_common.f90", "electron_radiation_kernel.f90", "electron_cooling_kernel.f90", "electron_seed_history_kernel.f90", "electron_transport_2d_kernel.f90", "calling_modules.f90", "FS_electron_charint_2d.f90"], OMP_FLAGS, OPENMP_LIBS),
        ("FS_electron_t2g1_1d", ele, ["../Constants.f90", "../Radiation/radiation_common.f90", "adaptive_resampling_mod.f90", "electron_common.f90", "electron_radiation_kernel.f90", "electron_cooling_kernel.f90", "calling_modules.f90", "FS_electron_t2g1_1d.f90"], OMP_FLAGS, OPENMP_LIBS),
        ("electron_get_y", ele, ["../Constants.f90", "../Radiation/radiation_common.f90", "adaptive_resampling_mod.f90", "electron_common.f90", "electron_radiation_kernel.f90", "electron_cooling_kernel.f90", "calling_modules.f90", "electron_get_y.f90"], OMP_FLAGS, OPENMP_LIBS),
        ("SED_interpolation", itp, ["../Constants.f90", "interpolation_common.f90", "SED_interpolation.f90"], OMP_FLAGS, OPENMP_LIBS),
        ("SED_interpolation_structured", itp, ["../Constants.f90", "interpolation_common.f90", "SED_interpolation_structured.f90"], OMP_FLAGS, OPENMP_LIBS),
        ("Annihilation", rad, ["../Constants.f90", "radiation_common.f90", "Annihilation.f90"], OMP_FLAGS, OPENMP_LIBS),
        ("Seed_reverse", rad, ["../Constants.f90", "radiation_common.f90", "Seed_reverse.f90"], OMP_FLAGS, OPENMP_LIBS),
        ("SSC_spec", rad, ["../Constants.f90", "radiation_common.f90", "SSC_spec.f90"], OMP_FLAGS, OPENMP_LIBS),
    ]
    module_aliases = {
        "FS_electron_weno5": "FS_electron_weno5_1d",
        "FS_electron_slc1": "FS_electron_slc1_1d",
        "FS_electron_charint": "FS_electron_charint_1d",
        "FS_electron_fullhide": "FS_electron_fullhide_1d",
        "FS_electron_t2g1": "FS_electron_t2g1_1d",
    }
    selected = set(args.modules or [])
    if selected:
        selected = {module_aliases.get(name, name) for name in selected}
    if selected:
        modules = [spec for spec in modules if spec[0] in selected]
        missing = selected.difference({spec[0] for spec in modules})
        if missing:
            raise SystemExit(f"Unknown module(s): {', '.join(sorted(missing))}")
    if args.clean:
        for directory in (src, dyn, ele, itp, rad):
            _clean_build_outputs(directory)

    print(f"Compile start ({_detect_build_platform()})")
    direct_ordered_build = {
        "electron_get_y",
        "FS_electron_weno5_1d",
        "FS_electron_slc1_1d",
        "FS_electron_charint_1d",
        "FS_electron_charint_2d",
        "FS_electron_fullhide_1d",
        "FS_electron_fullhide_2d",
        "FS_electron_t2g1_1d",
        "Seed_reverse",
        "SSC_spec",
    }
    for module_name, cwd, sources, fflags, extra_args in modules:
        print(f"Build {module_name}")
        if module_name in direct_ordered_build:
            elapsed = _build_fs_electron_ordered_fallback(
                root,
                module_name,
                cwd,
                sources,
                log_dir,
                args.verbose,
                fflags,
                extra_args,
            )
            print(f"Done {module_name}: {elapsed:.2f}s")
            continue
        try:
            elapsed = _build_module(
                module_name,
                cwd,
                sources,
                log_dir,
                args.verbose,
                fflags,
                extra_args,
            )
        except subprocess.CalledProcessError:
            if module_name not in {
                "FS_electron_weno5_1d",
                "FS_electron_slc1_1d",
                "FS_electron_charint_1d",
                "FS_electron_charint_2d",
                "FS_electron_fullhide_1d",
                "FS_electron_fullhide_2d",
                "FS_electron_t2g1_1d",
                "Seed_reverse",
                "SSC_spec",
            }:
                raise
            print(f"Retry {module_name} with ordered object build fallback")
            elapsed = _build_fs_electron_ordered_fallback(
                root,
                module_name,
                cwd,
                sources,
                log_dir,
                args.verbose,
                fflags,
                extra_args,
            )
        print(f"Done {module_name}: {elapsed:.2f}s")
    print("Compile complete!")


if __name__ == "__main__":
    main()
