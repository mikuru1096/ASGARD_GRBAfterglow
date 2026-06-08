from __future__ import annotations

import argparse
import hashlib
import os
import shlex
import shutil
import subprocess
import sys
import sysconfig
import time
from pathlib import Path


COMMON_FLAGS = "-Ofast -march=native -funroll-loops -ffast-math -fno-signed-zeros -fno-trapping-math -ffree-line-length-none"
OMP_FLAGS = f"-fopenmp {COMMON_FLAGS}"
OMP_LTO_FLAGS = f"{OMP_FLAGS} -flto"
OPENMP_LIBS = ["-lgomp"]
DIRECT_ORDERED_BUILD_MODULES = {
    "electron_forward_weno5_1d",
    "electron_forward_slc1_1d",
    "electron_forward_charint_1d",
    "electron_forward_charint_2d",
    "electron_forward_transport_2d_pic",
    "electron_forward_fullhide_1d",
    "electron_forward_fullhide_1d_hybrid",
    "electron_forward_transport_2d",
    "electron_forward_t2g1_1d",
    "hadronic_forward_1d",
    "hadronic_reverse_1d",
    "electron_radiation",
    "electron_reverse_kernel",
    "radiation_reverse_seed",
    "radiation_ssc_spectrum",
    "structured_jet_1d",
}
F2PY_ENTRYPOINTS = {
    "electron_forward_weno5_1d": ("fs_electron_weno5_1d",),
    "electron_forward_slc1_1d": ("fs_electron_slc1_1d",),
    "electron_forward_charint_1d": ("fs_electron_charint_1d",),
    "electron_forward_charint_2d": ("fs_electron_transport_2d_core",),
    "electron_forward_transport_2d_pic": ("fs_electron_transport_2d_pic_core",),
    "electron_forward_fullhide_1d": ("fs_electron_fullhide_1d",),
    "electron_forward_fullhide_1d_hybrid": ("fs_electron_fullhide_1d_hz",),
    "electron_forward_transport_2d": ("fs_electron_transport_2d_core",),
    "electron_forward_t2g1_1d": ("fs_electron_t2g1_1d",),
    "electron_radiation": ("get_nu_a", "get_syn_selected", "get_syn_transfer", "get_syn_polarization_selected"),
    "electron_reverse_kernel": ("electron_reverse_evolve",),
    "radiation_gamma_gamma_absorption": ("annihilation",),
    "radiation_reverse_seed": ("seed_reverse",),
    "radiation_ssc_spectrum": ("ssc_spec", "ssc_spec_nonuniform"),
    "hadronic_reverse_1d": ("fs_hadronic_reverse_1d",),
    "hadronic_forward_1d": (
        "fs_hadronic_1d",
        "fs_hadronic_proton_syn_shell",
        "fs_hadronic_syn_polarization_shell",
        "fs_hadronic_pgamma_operator_shell",
        "fs_hadronic_pair_production_shell",
        "fs_hadronic_pp_delta_shell",
        "fs_hadronic_bethe_heitler_shell",
        "fs_hadronic_hadronic_ic_shell",
        "fs_hadronic_species_transport_shell",
        "fs_hadronic_acceleration_shell",
        "fs_hadronic_secondary_radiation_shell",
        "fs_hadronic_decay_operator_shell",
        "fs_hadronic_pair_cascade_step",
        "fs_hadronic_pp_spectral_source",
        "fs_hadronic_quantum_syn_cooling_factor",
    ),
    "structured_jet_1d": ("structured_jet_flux_1d",),
}
DYNAMICS_COMMON_SOURCES = ("../Constants.f90", "dynamics_common.f90")
F2PY_SKIP_DYNAMICS_COMMON_INTERNALS = (
    "skip:",
    "dynamics_rk4_forward_ln_step",
    "dynamics_rk4_reverse",
    "dynamics_rk4_reverse_pre_m3",
    ":",
)
RADIATION_COMMON_SOURCES = (
    "../Constants.f90",
    "../Dynamics/dynamics_common.f90",
    "radiation_common.f90",
    "synchrotron_polarization_kernel.f90",
)
ELECTRON_COMMON_SOURCES = (
    "../Constants.f90",
    "../Dynamics/dynamics_common.f90",
    "../Radiation/radiation_common.f90",
    "../Radiation/synchrotron_polarization_kernel.f90",
    "electron_transport_common.f90",
    "adaptive_resampling_mod.f90",
    "electron_radiation_kernel.f90",
    "electron_injection_profiles.f90",
    "electron_common.f90",
    "electron_cooling_kernel.f90",
)
ELECTRON_RADIATION_SOURCES = (
    "../Constants.f90",
    "../Dynamics/dynamics_common.f90",
    "../Radiation/radiation_common.f90",
    "../Radiation/synchrotron_polarization_kernel.f90",
    "electron_transport_common.f90",
    "adaptive_resampling_mod.f90",
    "electron_radiation_kernel.f90",
)
ELECTRON_HISTORY_SOURCES = (
    *ELECTRON_COMMON_SOURCES,
    "electron_seed_history_kernel.f90",
)
ELECTRON_HISTORY_SOURCES_HZ = (
    *ELECTRON_HISTORY_SOURCES,
    "./slatec/dgamic.f",
    "./slatec/d1mach.f",
    "./slatec/d9gmic.f",
    "./slatec/d9gmit.f",
    "./slatec/d9lgic.f",
    "./slatec/d9lgit.f",
    "./slatec/dlgams.f",
    "./slatec/dlngam.f",
    "./slatec/d9lgmc.f",
    "./slatec/dgamma.f",
    "./slatec/dcsevl.f",
    "./slatec/initds.f",
    "./slatec/dgamlm.f",
    "./slatec/dlamch.f",
    "./slatec/lsame.f",
    "./slatec/xerclr.f",
    "./slatec/j4save.f",
    "./slatec/xermsg.f",
    "./slatec/fdump.f",
    "./slatec/xercnt.f",
    "./slatec/xerhlt.f",
    "./slatec/xerprn.f",
    "./slatec/xersve.f",
    "./slatec/i1mach.f",
    "./slatec/xgetua.f",
    "./slatec/dbesk.f",
    "./slatec/dasyik.f",
    "./slatec/dbesk0.f",
    "./slatec/dbesk1.f",
    "./slatec/dbsk0e.f",
    "./slatec/dbsk1e.f",
    "./slatec/dbsknu.f",
    "./slatec/dbesi0.f",
    "./slatec/dbesi1.f",
    "./slatec/dbsi0e.f",
    "./slatec/dbsi1e.f",
    "./slatec/specfun.f90",
    "hybrid_spectrum_kernel_fast.f90",
)
ELECTRON_2D_SOURCES = (
    *ELECTRON_COMMON_SOURCES,
    "electron_seed_history_kernel.f90",
    "electron_transport_2d_kernel.f90",
)
HADRONIC_1D_SOURCES = (
    "../Constants.f90",
    "../Dynamics/dynamics_common.f90",
    "../Radiation/radiation_common.f90",
    "../Radiation/quantum_synchrotron_kernel.f90",
    "../Radiation/synchrotron_polarization_kernel.f90",
    "hadronic_common.f90",
    "hadronic_transport_kernel.f90",
    "hadronic_transport_remap_kernel.f90",
    "hadronic_radiation_kernel.f90",
    "hadronic_interaction_kernel.f90",
    "hadronic_pair_production_kernel.f90",
    "hadronic_pp_kernel.f90",
    "hadronic_bethe_heitler_kernel.f90",
    "hadronic_hadronic_ic_kernel.f90",
    "hadronic_species_transport_kernel.f90",
    "hadronic_acceleration_kernel.f90",
    "hadronic_secondary_radiation_kernel.f90",
    "hadronic_decay_kernel.f90",
    "hadronic_pgamma_hummer_1d.f90",
    "hadronic_pair_cascade_kernel.f90",
    "hadronic_pp_models_kernel.f90",
)
STRUCTURED_JET_1D_SOURCES = (
    # Keep this list in module-topological order: every used module must be
    # compiled before files that import it in the ordered-object build path.
    "../Constants.f90",
    "../Dynamics/dynamics_common.f90",
    "../Radiation/radiation_common.f90",
    "../Radiation/quantum_synchrotron_kernel.f90",
    "../Radiation/synchrotron_polarization_kernel.f90",
    "../Electron/electron_transport_common.f90",
    "../Electron/adaptive_resampling_mod.f90",
    "../Electron/electron_radiation_kernel.f90",
    "../Electron/electron_injection_profiles.f90",
    "../Electron/electron_common.f90",
    "../Electron/electron_cooling_kernel.f90",
    "../Electron/electron_seed_history_kernel.f90",
    "../Hadronic/hadronic_common.f90",
    "../Hadronic/hadronic_transport_kernel.f90",
    "../Hadronic/hadronic_transport_remap_kernel.f90",
    "../Hadronic/hadronic_radiation_kernel.f90",
    "../Hadronic/hadronic_interaction_kernel.f90",
    "../Hadronic/hadronic_pair_production_kernel.f90",
    "../Hadronic/hadronic_pp_kernel.f90",
    "../Hadronic/hadronic_bethe_heitler_kernel.f90",
    "../Hadronic/hadronic_hadronic_ic_kernel.f90",
    "../Hadronic/hadronic_species_transport_kernel.f90",
    "../Hadronic/hadronic_acceleration_kernel.f90",
    "../Hadronic/hadronic_secondary_radiation_kernel.f90",
    "../Hadronic/hadronic_decay_kernel.f90",
    "../Hadronic/hadronic_pgamma_hummer_1d.f90",
    "../Hadronic/hadronic_pair_cascade_kernel.f90",
    "../Hadronic/hadronic_pp_models_kernel.f90",
    "../Dynamics/Dynamics_forward.f90",
    "../Dynamics/Dynamics_reverse.f90",
    "../Electron/electron_forward_fullhide_1d.f90",
    "../Electron/electron_reverse_kernel.f90",
    "../Radiation/radiation_ssc_spectrum.f90",
    "../Radiation/radiation_gamma_gamma_absorption.f90",
    "../Hadronic/hadronic_forward_1d.f90",
    "../Interpolation/interpolation_common.f90",
    "../Interpolation/SED_interpolation_structured.f90",
)

def _with_main(common_sources: tuple[str, ...], main_source: str) -> list[str]:
    return [*common_sources, main_source]


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


def _build_cache_root(root: Path) -> Path:
    override = os.environ.get("ASGARD_BUILD_CACHE_DIR")
    if override:
        return Path(override).expanduser().resolve()
    digest = hashlib.sha1(str(root).encode("utf-8")).hexdigest()[:16]
    return Path(os.environ.get("TMPDIR", "/tmp")) / "asgard_buildcache" / digest


def _module_output_paths(directory: Path, module_name: str) -> list[Path]:
    outputs: list[Path] = []
    for pattern in (f"{module_name}*.pyd", f"{module_name}*.so"):
        outputs.extend(directory.glob(pattern))
    return outputs


def _sources_newer_than(outputs: list[Path], cwd: Path, sources: list[str]) -> bool:
    if not outputs:
        return True
    output_mtime = min(path.stat().st_mtime for path in outputs)
    return any((cwd / source).resolve().stat().st_mtime > output_mtime for source in sources)


def _read_text_if_exists(path: Path) -> str:
    if not path.is_file():
        return ""
    return path.read_text(encoding="utf-8")


def _object_current(object_path: Path, source_path: Path, manifest_path: Path, manifest_text: str) -> bool:
    if not object_path.is_file():
        return False
    if _read_text_if_exists(manifest_path) != manifest_text:
        return False
    return object_path.stat().st_mtime >= source_path.stat().st_mtime


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


def _build_ordered_object_module(
    root: Path,
    module_name: str,
    cwd: Path,
    sources: list[str],
    log_dir: Path,
    verbose: bool,
    fflags: str | None,
    extra_args: list[str] | None,
    force: bool,
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

    outputs = _module_output_paths(cwd, module_name)
    if not force and not _sources_newer_than(outputs, cwd, sources):
        return 0.0

    build_dir = _build_cache_root(root) / "ordered_fallback" / module_name
    build_dir.mkdir(parents=True, exist_ok=True)

    compile_flags = shlex.split(fflags or "")
    compile_flags.append("-fPIC")
    compile_flags.extend(["-I", str(build_dir), "-J", str(build_dir)])
    object_paths: list[Path] = []
    start = time.perf_counter()
    dirty_seen = False
    for source in sources:
        source_path = (cwd / source).resolve()
        object_path = build_dir / f"{source_path.stem}.o"
        manifest_path = build_dir / f"{source_path.stem}.manifest"
        manifest_text = "\n".join(
            [
                f"source={source_path}",
                f"fc={fc}",
                f"flags={shlex.join(compile_flags)}",
            ]
        )
        if dirty_seen or not _object_current(object_path, source_path, manifest_path, manifest_text):
            command = [fc, "-c", *compile_flags, str(source_path), "-o", str(object_path)]
            _run_command(command, cwd, env, log_dir / f"{module_name}_fallback_compile_{source_path.stem}.log", verbose)
            manifest_path.write_text(manifest_text, encoding="utf-8")
            dirty_seen = True
        object_paths.append(object_path)

    pyf_path = build_dir / f"{module_name}.pyf"
    main_source_name = Path(sources[-1]).name
    main_source_path = (cwd / main_source_name).resolve()
    signature_source_path = build_dir / main_source_name
    if not signature_source_path.is_file() or signature_source_path.stat().st_mtime < main_source_path.stat().st_mtime:
        shutil.copy2(main_source_path, signature_source_path)
    entry_names = list(F2PY_ENTRYPOINTS.get(module_name, (module_name.lower(),)))
    wrapper_manifest_path = build_dir / "wrapper.manifest"
    wrapper_manifest_text = "\n".join(
        [
            f"main={main_source_path}",
            f"entries={','.join(entry_names)}",
            f"sources={'|'.join(str((cwd / source).resolve()) for source in sources)}",
            f"python={sys.executable}",
        ]
    )
    wrapper_outputs = [build_dir / f"{module_name}module.c"]
    wrapper_current = (
        _read_text_if_exists(wrapper_manifest_path) == wrapper_manifest_text
        and all(path.is_file() for path in wrapper_outputs)
        and all(path.stat().st_mtime >= (cwd / main_source_name).resolve().stat().st_mtime for path in wrapper_outputs)
    )
    if not wrapper_current:
        signature_command = [
            sys.executable,
            "-m",
            "numpy.f2py",
            "-m",
            module_name,
            "-h",
            pyf_path.name,
            "--overwrite-signature",
            signature_source_path.name,
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
        wrapper_manifest_path.write_text(wrapper_manifest_text, encoding="utf-8")

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
        c_flags = [
            "-fPIC",
            f"-I{py_include}",
            f"-I{py_platinclude}",
            f"-I{np.get_include()}",
            f"-I{f2py_src}",
        ]
        manifest_path = build_dir / f"{source.stem}.manifest"
        manifest_text = "\n".join([f"source={source}", f"cc={cc}", f"flags={shlex.join(c_flags)}"])
        if not _object_current(object_path, source, manifest_path, manifest_text):
            compile_c = [
                cc,
                "-c",
                *c_flags,
                str(source),
                "-o",
                str(object_path),
            ]
            _run_command(compile_c, build_dir, env, log_dir / f"{module_name}_fallback_compile_{source.stem}.log", verbose)
            manifest_path.write_text(manifest_text, encoding="utf-8")
        c_objects.append(object_path)

    wrapper_sources = [
        build_dir / f"{module_name}-f2pywrappers.f",
        build_dir / f"{module_name}-f2pywrappers2.f90",
    ]
    for wrapper_source in wrapper_sources:
        if not wrapper_source.is_file():
            continue
        wrapper_object = build_dir / f"{wrapper_source.stem}.o"
        manifest_path = build_dir / f"{wrapper_source.stem}.manifest"
        manifest_text = "\n".join([f"source={wrapper_source}", f"fc={fc}", "flags=-fPIC"])
        if not _object_current(wrapper_object, wrapper_source, manifest_path, manifest_text):
            compile_wrapper = [fc, "-c", "-fPIC", str(wrapper_source), "-o", str(wrapper_object)]
            _run_command(
                compile_wrapper,
                build_dir,
                env,
                log_dir / f"{module_name}_fallback_compile_{wrapper_source.stem}.log",
                verbose,
            )
            manifest_path.write_text(manifest_text, encoding="utf-8")
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
    built_outputs = _module_output_paths(build_dir, module_name)
    if not built_outputs:
        raise RuntimeError(f"{module_name} fallback build did not produce an extension in {build_dir}.")
    for built_path in built_outputs:
        target_path = cwd / built_path.name
        if target_path.exists():
            target_path.unlink()
        shutil.copy2(built_path, target_path)
        if not target_path.is_file():
            raise RuntimeError(f"{module_name} fallback build failed to copy {target_path}.")
        built_path.unlink()
    return time.perf_counter() - start


def _build_module(
    module_name: str,
    cwd: Path,
    sources: list[str],
    log_dir: Path,
    verbose: bool,
    fflags: str | None = None,
    extra_args: list[str] | None = None,
    force: bool = False,
) -> float:
    env = _prepare_build_env()
    fflags_override = env.get("ASGARD_FFLAGS_OVERRIDE")
    if fflags_override:
        fflags = fflags_override
    if fflags is not None:
        env["FFLAGS"] = fflags
        env["F90FLAGS"] = fflags

    outputs = _module_output_paths(cwd, module_name)
    if not force and not _sources_newer_than(outputs, cwd, sources):
        return 0.0

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
    parser.add_argument("--lto", action="store_true", help="enable link-time optimization for release/performance builds")
    args = parser.parse_args()

    root = Path(__file__).resolve().parent
    log_dir = root / ".buildcache" / "logs"
    src = root / "src"
    dyn = src / "Dynamics"
    ele = src / "Electron"
    itp = src / "Interpolation"
    rad = src / "Radiation"
    had = src / "Hadronic"
    structured = src / "Structured"
    omp_flags = OMP_LTO_FLAGS if args.lto else OMP_FLAGS
    modules = [
        ("Constants", src, ["Constants.f90"], None, None),
        ("Dynamics_reverse", dyn, _with_main(DYNAMICS_COMMON_SOURCES, "Dynamics_reverse.f90"), COMMON_FLAGS, None),
        ("Dynamics_forward", dyn, _with_main(DYNAMICS_COMMON_SOURCES, "Dynamics_forward.f90"), COMMON_FLAGS, None),
        ("electron_forward_weno5_1d", ele, _with_main(ELECTRON_COMMON_SOURCES, "electron_forward_weno5_1d.f90"), omp_flags, OPENMP_LIBS),
        ("electron_forward_slc1_1d", ele, _with_main(ELECTRON_COMMON_SOURCES, "electron_forward_slc1_1d.f90"), omp_flags, OPENMP_LIBS),
        ("electron_forward_charint_1d", ele, _with_main(ELECTRON_COMMON_SOURCES, "electron_forward_charint_1d.f90"), omp_flags, OPENMP_LIBS),
        ("electron_forward_fullhide_1d", ele, _with_main(ELECTRON_HISTORY_SOURCES, "electron_forward_fullhide_1d.f90"), omp_flags, OPENMP_LIBS),
        #
        ("electron_forward_fullhide_1d_hybrid", ele, _with_main(ELECTRON_HISTORY_SOURCES_HZ, "electron_forward_fullhide_1d_hybrid.f90"), omp_flags, OPENMP_LIBS),
        ("electron_forward_transport_2d", ele, _with_main(ELECTRON_2D_SOURCES, "electron_forward_transport_2d.f90"), omp_flags, OPENMP_LIBS),
        ("electron_forward_charint_2d", ele, _with_main(ELECTRON_2D_SOURCES, "electron_forward_transport_2d.f90"), omp_flags, OPENMP_LIBS),
        ("electron_forward_transport_2d_pic", ele, _with_main(ELECTRON_2D_SOURCES, "electron_forward_transport_2d_pic.f90"), omp_flags, OPENMP_LIBS),
        ("electron_forward_t2g1_1d", ele, _with_main(ELECTRON_COMMON_SOURCES, "electron_forward_t2g1_1d.f90"), omp_flags, OPENMP_LIBS),
        ("electron_radiation", ele, ELECTRON_RADIATION_SOURCES, omp_flags, OPENMP_LIBS),
        ("electron_reverse_kernel", ele, _with_main(ELECTRON_COMMON_SOURCES, "electron_reverse_kernel.f90"), omp_flags, OPENMP_LIBS),
        ("SED_interpolation", itp, ["../Constants.f90", "interpolation_common.f90", "SED_interpolation.f90"], omp_flags, OPENMP_LIBS),
        ("SED_interpolation_structured", itp, ["../Constants.f90", "interpolation_common.f90", "SED_interpolation_structured.f90"], omp_flags, OPENMP_LIBS),
        ("radiation_gamma_gamma_absorption", rad, _with_main(RADIATION_COMMON_SOURCES, "radiation_gamma_gamma_absorption.f90"), omp_flags, OPENMP_LIBS),
        ("radiation_reverse_seed", rad, _with_main(RADIATION_COMMON_SOURCES, "radiation_reverse_seed.f90"), omp_flags, OPENMP_LIBS),
        ("radiation_ssc_spectrum", rad, _with_main(RADIATION_COMMON_SOURCES, "radiation_ssc_spectrum.f90"), omp_flags, OPENMP_LIBS),
        ("hadronic_forward_1d", had, _with_main(HADRONIC_1D_SOURCES, "hadronic_forward_1d.f90"), omp_flags, OPENMP_LIBS),
        ("hadronic_reverse_1d", had, _with_main(HADRONIC_1D_SOURCES, "hadronic_reverse_1d.f90"), omp_flags, OPENMP_LIBS),
        ("structured_jet_1d", structured, _with_main(STRUCTURED_JET_1D_SOURCES, "structured_jet_1d.f90"), omp_flags, OPENMP_LIBS),
    ]
    module_aliases = {
        "electron_forward_weno5": "electron_forward_weno5_1d",
        "electron_forward_slc1": "electron_forward_slc1_1d",
        "electron_forward_charint": "electron_forward_charint_1d",
        "electron_forward_fullhide": "electron_forward_fullhide_1d",
        "electron_forward_t2g1": "electron_forward_t2g1_1d",
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
        for directory in (src, dyn, ele, had, itp, rad, structured):
            _clean_build_outputs(directory)

    print(f"Compile start ({_detect_build_platform()})")
    for module_name, cwd, sources, fflags, extra_args in modules:
        print(f"Build {module_name}")
        module_extra_args = list(extra_args or [])
        if module_name in DIRECT_ORDERED_BUILD_MODULES:
            elapsed = _build_ordered_object_module(
                root,
                module_name,
                cwd,
                sources,
                log_dir,
                args.verbose,
                fflags,
                module_extra_args,
                args.force,
            )
            print(f"Done {module_name}: {elapsed:.2f}s")
            continue
        if any(Path(source).name == "dynamics_common.f90" for source in sources):
            module_extra_args = [*F2PY_SKIP_DYNAMICS_COMMON_INTERNALS, *module_extra_args]
        elapsed = _build_module(
            module_name,
            cwd,
            sources,
            log_dir,
            args.verbose,
            fflags,
            module_extra_args,
            args.force,
        )
        print(f"Done {module_name}: {elapsed:.2f}s")
    print("Compile complete!")


if __name__ == "__main__":
    main()
