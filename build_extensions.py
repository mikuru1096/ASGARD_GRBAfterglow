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
from dataclasses import dataclass
from pathlib import Path


COMMON_FLAGS = "-Ofast -march=native -funroll-loops -ffast-math -fno-signed-zeros -fno-trapping-math -ffree-line-length-none"
OMP_FLAGS = f"-fopenmp {COMMON_FLAGS}"
OPENMP_LIBS = ["-lgomp"]


@dataclass(frozen=True)
class ModuleSpec:
    name: str
    cwd: Path
    sources: list[str]
    fflags: str | None = None
    extra_args: list[str] | None = None
    ordered: bool = False
    entrypoints: tuple[str, ...] = ()
    aliases: tuple[str, ...] = ()


DYNAMICS_DENSITY_SOURCES = ("../Constants.f90", "dynamics_density_profile.f90")
DYNAMICS_FORWARD_SOURCES = DYNAMICS_DENSITY_SOURCES
DYNAMICS_REVERSE_SOURCES = (
    *DYNAMICS_DENSITY_SOURCES,
    "reverse_shock_state.f90",
    "reverse_shock_mhd_jump.f90",
    "reverse_jump_conditions.f90",
    "reverse_rhs.f90",
    "reverse_shock.f90",
)
F2PY_SKIP_RADIATION_COMMON_INTERNALS = (
    "pair_delta",
    "pair_sigma",
    "syn_kernel",
)
RADIATION_COMMON_SOURCES = (
    "../Constants.f90",
    "rad_common.f90",
    "syn_polarization.f90",
)
ELECTRON_COMMON_SOURCES = (
    "../Constants.f90",
    "../Dynamics/dynamics_density_profile.f90",
    "../Radiation/rad_common.f90",
    "electron_energy_coordinate_common.f90",
    "electron_transport_common.f90",
    "electron_radiation_kernel.f90",
    "electron_injection_profiles.f90",
    "electron_shell_transport_common.f90",
    "electron_common.f90",
    "electron_cooling_ssa_kernel.f90",
    "electron_cooling_ic_kernel.f90",
    "electron_cooling_y_kernel.f90",
    "electron_cooling_kernel.f90",
)
ELECTRON_RADIATION_SOURCES = (
    "../Constants.f90",
    "../Radiation/rad_common.f90",
    "electron_transport_common.f90",
    "electron_radiation_kernel.f90",
)
ELECTRON_HISTORY_SOURCES = (
    *ELECTRON_COMMON_SOURCES,
    "electron_seed_history_kernel.f90",
)
ELECTRON_HYBRID_SOURCES = (
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
    "./slatec/dbesk0.f",
    "./slatec/dbesk1.f",
    "./slatec/dbsk0e.f",
    "./slatec/dbsk1e.f",
    "./slatec/dbesi0.f",
    "./slatec/dbesi1.f",
    "./slatec/dbsi0e.f",
    "./slatec/dbsi1e.f",
    "hybrid_special.f90",
    "hybrid_spectrum.f90",
)
ELECTRON_HISTORY_SOURCES_HZ = (
    *ELECTRON_HISTORY_SOURCES,
    *ELECTRON_HYBRID_SOURCES,
)
ELECTRON_2D_SOURCES = (
    *ELECTRON_COMMON_SOURCES,
    "electron_seed_history_kernel.f90",
    "electron_transport_2d_kernel.f90",
)
ELECTRON_DG_1D_SOURCES = (
    *ELECTRON_COMMON_SOURCES,
    *ELECTRON_HYBRID_SOURCES,
    "electron_dg_transport.f90",
)
ELECTRON_REVERSE_SOURCES = (
    *ELECTRON_COMMON_SOURCES,
    "electron_dg_transport.f90",
)
HADRONIC_COMMON_SOURCES = (
    "../Constants.f90",
    "../Radiation/rad_common.f90",
    "../Radiation/quantum_synch.f90",
    "../Radiation/syn_polarization.f90",
    "hadronic_base.f90",
    "hadronic_transport_kernel.f90",
    "hadronic_transport_remap_kernel.f90",
    "hadronic_rad.f90",
    "hadronic_pg.f90",
    "hadronic_pair.f90",
    "hadronic_pp.f90",
    "hadronic_bh.f90",
    "hadronic_ic.f90",
    "hadronic_species.f90",
    "hadronic_accel.f90",
    "hadronic_secondary.f90",
    "hadronic_decay.f90",
    "hadronic_hummer.f90",
    "../Electron/electron_radiation_kernel.f90",
    "hadronic_cascade.f90",
    "pp_models.f90",
    "hadronic_shell.f90",
)
HADRONIC_1D_SOURCES = (
    *HADRONIC_COMMON_SOURCES,
    "hadronic_formal.f90",
)
STRUCTURED_JET_1D_SOURCES = (
    # Keep this list in module-topological order: every used module must be
    # compiled before files that import it in the ordered-object build path.
    "../Constants.f90",
    "../Dynamics/dynamics_density_profile.f90",
    "../Dynamics/reverse_shock_state.f90",
    "../Dynamics/reverse_shock_mhd_jump.f90",
    "../Radiation/rad_common.f90",
    "../Radiation/quantum_synch.f90",
    "../Radiation/syn_polarization.f90",
    "../Electron/electron_energy_coordinate_common.f90",
    "../Electron/electron_transport_common.f90",
    "../Electron/electron_radiation_kernel.f90",
    "../Electron/electron_injection_profiles.f90",
    "../Electron/electron_shell_transport_common.f90",
    "../Electron/electron_common.f90",
    *(f"../Electron/{source}" for source in ELECTRON_HYBRID_SOURCES),
    "../Electron/electron_cooling_ssa_kernel.f90",
    "../Electron/electron_cooling_ic_kernel.f90",
    "../Electron/electron_cooling_y_kernel.f90",
    "../Electron/electron_cooling_kernel.f90",
    "../Electron/electron_dg_transport.f90",
    "../Electron/electron_seed_history_kernel.f90",
    "../Hadronic/hadronic_base.f90",
    "../Hadronic/hadronic_transport_kernel.f90",
    "../Hadronic/hadronic_transport_remap_kernel.f90",
    "../Hadronic/hadronic_rad.f90",
    "../Hadronic/hadronic_pg.f90",
    "../Hadronic/hadronic_pair.f90",
    "../Hadronic/hadronic_pp.f90",
    "../Hadronic/hadronic_bh.f90",
    "../Hadronic/hadronic_ic.f90",
    "../Hadronic/hadronic_species.f90",
    "../Hadronic/hadronic_accel.f90",
    "../Hadronic/hadronic_secondary.f90",
    "../Hadronic/hadronic_decay.f90",
    "../Hadronic/hadronic_hummer.f90",
    "../Hadronic/hadronic_cascade.f90",
    "../Hadronic/pp_models.f90",
    "../Dynamics/Dynamics_forward.f90",
    "../Dynamics/reverse_jump_conditions.f90",
    "../Dynamics/reverse_rhs.f90",
    "../Dynamics/reverse_shock.f90",
    "../Electron/electron_forward_dg_1d.f90",
    "../Electron/electron_forward_fullhide_1d.f90",
    "../Electron/electron_reverse_kernel.f90",
    "../Radiation/ssc_spectrum.f90",
    "../Radiation/pair_absorption.f90",
    "../Hadronic/hadronic_shell.f90",
    "../Hadronic/hadronic_formal.f90",
    "../Hadronic/hadronic_forward_1d.f90",
    "../Interpolation/interpolation_common.f90",
    "../Interpolation/SED_interpolation_structured.f90",
)

def _with_main(common_sources: tuple[str, ...], main_source: str) -> list[str]:
    return [*common_sources, main_source]


def _prepare_build_env() -> dict[str, str]:
    env = os.environ.copy()
    env["PYTHONUTF8"] = "1"
    env["PYTHONIOENCODING"] = "utf-8"
    env["LANG"] = "C.UTF-8"
    env["LC_ALL"] = "C.UTF-8"
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


def _effective_fflags(env: dict[str, str], fflags: str | None) -> str:
    override = env.get("ASGARD_FFLAGS_OVERRIDE")
    if override:
        return override
    return fflags or ""


def _module_output_paths(directory: Path, module_name: str) -> list[Path]:
    return [path for pattern in (f"{module_name}*.pyd", f"{module_name}*.so") for path in directory.glob(pattern)]


def _sources_newer_than(outputs: list[Path], cwd: Path, sources: list[str]) -> bool:
    if not outputs:
        return True
    output_mtime = min(path.stat().st_mtime for path in outputs)
    return any((cwd / source).resolve().stat().st_mtime > output_mtime for source in sources)


def _read_text_if_exists(path: Path) -> str:
    if not path.is_file():
        return ""
    return path.read_text(encoding="utf-8")


def _write_f2py_signature_source(source_path: Path, target_path: Path) -> None:
    lines = source_path.read_text(encoding="utf-8").splitlines()
    stripped_lines = [
        line for line in lines
        if not (
            line.split("!", 1)[0].strip().lower() in {"block", "end block", "end associate"}
            or line.split("!", 1)[0].strip().lower().startswith("associate(")
            or line.split("!", 1)[0].strip().lower().startswith("associate (")
        )
    ]
    if len(stripped_lines) == len(lines):
        shutil.copy2(source_path, target_path)
        return
    target_path.write_text("\n".join(stripped_lines) + "\n", encoding="utf-8")


def _object_current(object_path: Path, source_path: Path, manifest_path: Path, manifest_text: str) -> bool:
    if not object_path.is_file():
        return False
    if _read_text_if_exists(manifest_path) != manifest_text:
        return False
    return object_path.stat().st_mtime >= source_path.stat().st_mtime


def _run_command(
    command: list[str],
    cwd: Path,
    env: dict[str, str],
    log_path: Path,
    failure_label: str | None = None,
) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(
        command,
        cwd=cwd,
        env=env,
        check=False,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(
        f"cwd: {cwd}\ncommand: {' '.join(command)}\nreturncode: {result.returncode}\n\n"
        f"--- stdout ---\n{result.stdout or ''}\n\n--- stderr ---\n{result.stderr or ''}",
        encoding="utf-8",
    )
    if result.returncode != 0:
        output = "\n".join(part for part in (result.stdout, result.stderr) if part)
        tail = "\n".join(output.strip().splitlines()[-120:])
        print(f"Build failed for {failure_label or ' '.join(command[:3])}. Full log: {log_path}")
        if tail:
            print(tail)
        raise subprocess.CalledProcessError(result.returncode, command)
    return result


def _build_ordered_object_module(
    root: Path,
    module_name: str,
    cwd: Path,
    sources: list[str],
    entry_names: tuple[str, ...],
    log_dir: Path,
    fflags: str | None,
    extra_args: list[str] | None,
    force: bool,
) -> float:
    env = _prepare_build_env()
    fflags = _effective_fflags(env, fflags)
    if fflags:
        env["FFLAGS"] = fflags
        env["F90FLAGS"] = fflags

    fc = env.get("FC") or shutil.which("gfortran", path=env["PATH"])
    if not fc:
        raise RuntimeError(f"{module_name} ordered-object build requires gfortran in PATH or FC.")

    outputs = _module_output_paths(cwd, module_name)
    if not force and not _sources_newer_than(outputs, cwd, sources):
        return 0.0

    cache_override = os.environ.get("ASGARD_BUILD_CACHE_DIR")
    cache_digest = hashlib.sha1(str(root).encode("utf-8")).hexdigest()[:16]
    cache_root = (
        Path(cache_override).expanduser().resolve()
        if cache_override
        else Path(os.environ.get("TMPDIR", "/tmp")) / "asgard_buildcache" / cache_digest
    )
    build_dir = cache_root / "ordered_object" / module_name
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
        manifest_text = f"source={source_path}\nfc={fc}\nflags={shlex.join(compile_flags)}"
        if dirty_seen or not _object_current(object_path, source_path, manifest_path, manifest_text):
            command = [fc, "-c", *compile_flags, str(source_path), "-o", str(object_path)]
            _run_command(command, cwd, env, log_dir / f"{module_name}_ordered_compile_{source_path.stem}.log")
            manifest_path.write_text(manifest_text, encoding="utf-8")
            dirty_seen = True
        object_paths.append(object_path)

    pyf_path = build_dir / f"{module_name}.pyf"
    main_source_name = Path(sources[-1]).name
    main_source_path = (cwd / main_source_name).resolve()
    signature_source_path = build_dir / main_source_name
    if not signature_source_path.is_file() or signature_source_path.stat().st_mtime < main_source_path.stat().st_mtime:
        _write_f2py_signature_source(main_source_path, signature_source_path)
    wrapper_manifest_path = build_dir / "wrapper.manifest"
    wrapper_manifest_text = (
        f"main={main_source_path}\nentries={','.join(entry_names)}\n"
        f"sources={'|'.join(str((cwd / source).resolve()) for source in sources)}\npython={sys.executable}"
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
        _run_command(signature_command, build_dir, env, log_dir / f"{module_name}_ordered_signature.log")

        wrapper_command = [
            sys.executable,
            "-m",
            "numpy.f2py",
            pyf_path.name,
        ]
        _run_command(wrapper_command, build_dir, env, log_dir / f"{module_name}_ordered_wrapper.log")
        wrapper_manifest_path.write_text(wrapper_manifest_text, encoding="utf-8")

    import numpy as np

    cc = env.get("CC") or shutil.which("cc", path=env["PATH"])
    if not cc:
        raise RuntimeError(f"{module_name} ordered-object build requires a C compiler in PATH or CC.")
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
            _run_command(compile_c, build_dir, env, log_dir / f"{module_name}_ordered_compile_{source.stem}.log")
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
        wrapper_flags = ["-fPIC"]
        manifest_text = "\n".join([f"source={wrapper_source}", f"fc={fc}", f"flags={shlex.join(wrapper_flags)}"])
        if not _object_current(wrapper_object, wrapper_source, manifest_path, manifest_text):
            compile_wrapper = [fc, "-c", *wrapper_flags, str(wrapper_source), "-o", str(wrapper_object)]
            _run_command(
                compile_wrapper,
                build_dir,
                env,
                log_dir / f"{module_name}_ordered_compile_{wrapper_source.stem}.log",
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
    _run_command(link_command, build_dir, env, log_dir / f"{module_name}_ordered_link.log")
    built_outputs = _module_output_paths(build_dir, module_name)
    if not built_outputs:
        raise RuntimeError(f"{module_name} ordered-object build did not produce an extension in {build_dir}.")
    for built_path in built_outputs:
        target_path = cwd / built_path.name
        if target_path.exists():
            target_path.unlink()
        shutil.copy2(built_path, target_path)
        if not target_path.is_file():
            raise RuntimeError(f"{module_name} ordered-object build failed to copy {target_path}.")
        built_path.unlink()
    return time.perf_counter() - start


def _build_module(
    module_name: str,
    cwd: Path,
    sources: list[str],
    log_dir: Path,
    fflags: str | None = None,
    extra_args: list[str] | None = None,
    force: bool = False,
) -> float:
    env = _prepare_build_env()
    fflags = _effective_fflags(env, fflags)
    if fflags:
        env["FFLAGS"] = fflags
        env["F90FLAGS"] = fflags

    outputs = _module_output_paths(cwd, module_name)
    if not force and not _sources_newer_than(outputs, cwd, sources):
        return 0.0

    command = [sys.executable, "-m", "numpy.f2py", "-m", module_name, "-c", *sources]
    if fflags:
        command.extend([f"--f77flags={fflags}", f"--f90flags={fflags}"])
    if extra_args:
        command.extend(extra_args)
    start = time.perf_counter()
    _run_command(command, cwd, env, log_dir / f"{module_name}.log", failure_label=module_name)
    return time.perf_counter() - start


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true", help="rebuild all selected modules")
    parser.add_argument("--module", action="append", dest="modules", help="only build the named module; can be repeated")
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
    omp_flags = OMP_FLAGS
    module_specs = [
        ModuleSpec("Constants", src, ["Constants.f90"]),
        ModuleSpec("Dynamics_reverse", dyn, list(DYNAMICS_REVERSE_SOURCES), COMMON_FLAGS, ["only:", "rs_prompt_jump", "dynamics_reverse", ":"]),
        ModuleSpec("Dynamics_forward", dyn, _with_main(DYNAMICS_FORWARD_SOURCES, "Dynamics_forward.f90"), COMMON_FLAGS, ["only:", "dynamics_forward", ":"]),
        ModuleSpec("electron_forward_weno5_1d", ele, _with_main(ELECTRON_COMMON_SOURCES, "electron_forward_weno5_1d.f90"), omp_flags, OPENMP_LIBS, True, ("fs_weno5_1d",), ("electron_forward_weno5",)),
        ModuleSpec("electron_forward_slc1_1d", ele, _with_main(ELECTRON_COMMON_SOURCES, "electron_forward_slc1_1d.f90"), omp_flags, OPENMP_LIBS, True, ("fs_slc1_1d",), ("electron_forward_slc1",)),
        ModuleSpec("electron_forward_charint_1d", ele, _with_main(ELECTRON_COMMON_SOURCES, "electron_forward_charint_1d.f90"), omp_flags, OPENMP_LIBS, True, ("fs_charint_1d",), ("electron_forward_charint",)),
        ModuleSpec(
            "electron_forward_dg_1d",
            ele,
            _with_main(ELECTRON_DG_1D_SOURCES, "electron_forward_dg_1d.f90"),
            omp_flags,
            OPENMP_LIBS,
            True,
            ("fs_dg_1d",),
            ("electron_forward_dg",),
        ),
        ModuleSpec(
            "electron_forward_fullhide_1d",
            ele,
            _with_main(ELECTRON_HISTORY_SOURCES, "electron_forward_fullhide_1d.f90"),
            omp_flags,
            OPENMP_LIBS,
            True,
            ("fs_fullhide_1d", "fs_fullhide_bh_1d", "fs_fullhide_coupled", "fs_fullhide_coupled_bh"),
            ("electron_forward_fullhide",),
        ),
        ModuleSpec("electron_forward_fullhide_1d_hybrid", ele, _with_main(ELECTRON_HISTORY_SOURCES_HZ, "electron_forward_fullhide_1d_hybrid.f90"), omp_flags, OPENMP_LIBS, True, ("fs_fullhide_hz",)),
        ModuleSpec("electron_forward_transport_2d", ele, _with_main(ELECTRON_2D_SOURCES, "electron_forward_transport_2d.f90"), omp_flags, OPENMP_LIBS, True, ("fs_transport_2d",), ("electron_forward_charint_2d",)),
        ModuleSpec("electron_forward_t2g1_1d", ele, _with_main(ELECTRON_COMMON_SOURCES, "electron_forward_t2g1_1d.f90"), omp_flags, OPENMP_LIBS, True, ("fs_t2g1_1d",), ("electron_forward_t2g1",)),
        ModuleSpec("electron_radiation", ele, list(ELECTRON_RADIATION_SOURCES), omp_flags, OPENMP_LIBS, True, ("nua_solve", "syn_state", "syn_transfer")),
        ModuleSpec(
            "electron_reverse_kernel",
            ele,
            _with_main(ELECTRON_REVERSE_SOURCES, "electron_reverse_kernel.f90"),
            omp_flags,
            OPENMP_LIBS,
            True,
            (
                "electron_reverse_evolve",
                "multiple_synch",
                "branch_reaccel",
            ),
        ),
        ModuleSpec(
            "SED_interpolation",
            itp,
            ["../Constants.f90", "../Radiation/rad_common.f90", "interpolation_common.f90", "SED_interpolation.f90"],
            omp_flags,
            [*OPENMP_LIBS, "skip:", "accum_radial_batch", ":"],
        ),
        ModuleSpec(
            "SED_interpolation_structured",
            itp,
            ["../Constants.f90", "interpolation_common.f90", "SED_interpolation_structured.f90"],
            omp_flags,
            [*OPENMP_LIBS, "skip:", "accum_radial_batch", ":"],
        ),
        ModuleSpec("pair_absorption", rad, _with_main(RADIATION_COMMON_SOURCES, "pair_absorption.f90"), omp_flags, OPENMP_LIBS, False, ("annihilation",)),
        ModuleSpec("ssc_spectrum", rad, _with_main(RADIATION_COMMON_SOURCES, "ssc_spectrum.f90"), omp_flags, OPENMP_LIBS, True, ("ssc_spec", "ssc_spec_nonuniform")),
        ModuleSpec(
            "hadronic_forward_1d",
            had,
            _with_main(HADRONIC_1D_SOURCES, "hadronic_forward_1d.f90"),
            omp_flags,
            OPENMP_LIBS,
            True,
            (
                "hadronic_1d",
                "formal_transport_1d",
                "bh_support",
                "had_syn_pol",
                "pg_operator",
                "pair_production",
                "pp_shell",
                "bethe_heitler",
                "hadronic_ic",
                "decay_operator",
                "cascade_sequence",
            ),
        ),
        ModuleSpec("hadronic_reverse_1d", had, _with_main(HADRONIC_COMMON_SOURCES, "hadronic_reverse_1d.f90"), omp_flags, OPENMP_LIBS, True, ("reverse_hadronic_1d",)),
        ModuleSpec("structured_jet_1d", structured, _with_main(STRUCTURED_JET_1D_SOURCES, "structured_jet_1d.f90"), omp_flags, OPENMP_LIBS, True, ("jet_flux_1d",)),
    ]
    module_aliases = {alias: spec.name for spec in module_specs for alias in spec.aliases}
    selected = set(args.modules or [])
    if selected:
        selected = {module_aliases.get(name, name) for name in selected}
        module_specs = [spec for spec in module_specs if spec.name in selected]
        missing = selected.difference({spec.name for spec in module_specs})
        if missing:
            raise SystemExit(f"Unknown module(s): {', '.join(sorted(missing))}")
    platform = {"nt": "windows-mingw", "posix": "linux-gfortran"}.get(os.name, f"unsupported:{os.name}")
    print(f"Compile start ({platform})")
    for spec in module_specs:
        print(f"Build {spec.name}")
        module_extra_args = list(spec.extra_args or [])
        if spec.ordered:
            elapsed = _build_ordered_object_module(
                root,
                spec.name,
                spec.cwd,
                spec.sources,
                spec.entrypoints,
                log_dir,
                spec.fflags,
                module_extra_args,
                args.force,
            )
            print(f"Done {spec.name}: {elapsed:.2f}s")
            continue
        skip_names: list[str] = []
        if any(Path(source).name == "rad_common.f90" for source in spec.sources):
            skip_names.extend(F2PY_SKIP_RADIATION_COMMON_INTERNALS)
        if skip_names:
            module_extra_args = ["skip:", *skip_names, ":", *module_extra_args]
        elapsed = _build_module(
            spec.name,
            spec.cwd,
            spec.sources,
            log_dir,
            spec.fflags,
            module_extra_args,
            args.force,
        )
        print(f"Done {spec.name}: {elapsed:.2f}s")
    print("Compile complete!")


if __name__ == "__main__":
    main()
