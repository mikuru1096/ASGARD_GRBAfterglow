from __future__ import annotations

import json
import sys
from pathlib import Path
from time import perf_counter

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_models import default_num_threads
from asgard_presets import build_baseline_config
from asgard_runtime import solve_dynamics, solve_electron
from asgard_setup import build_simulation_setup
import src.Electron.FS_electron_fullhide as electron_module


OUTPUT_DIR = ROOT / "output" / "vegasafterglow_doc"
SERIAL_THREADS = 1
PARALLEL_THREADS = default_num_threads()
PROFILE_REPEATS = 12


def _build_shell_state():
    config = build_baseline_config(num_threads=SERIAL_THREADS)
    setup = build_simulation_setup(config)
    boundary = setup.boundary
    v_seed = np.asfortranarray(setup.seed_frequency_hz)

    dynamics = solve_dynamics(boundary, config)
    electron = solve_electron(boundary, dynamics, v_seed, config)

    i_shell = max(2, len(dynamics.radius) // 2)
    r_loc = float(dynamics.radius[i_shell - 1])
    gamma_loc = float(0.5 * (dynamics.r_gamma[i_shell - 1] + dynamics.r_gamma[i_shell]))
    if config.a_star > 0.0:
        d_ne = config.a_star * 3.0e35 / r_loc**2
    else:
        d_ne = config.d_ne
    db = 0.39 * np.sqrt(config.epsilon_b * d_ne * (gamma_loc * (gamma_loc - 1.0)))

    return {
        "r_loc": r_loc,
        "db": db,
        "gam_e": np.asfortranarray(electron.gam_e),
        "d_n_gam_e": np.asfortranarray(electron.d_n_gam_e[:, i_shell - 1]),
        "p_syn": np.asfortranarray(electron.l_syn_spec[:, i_shell]),
        "seed_syn": np.asfortranarray(electron.seed_syn[:, i_shell]),
        "v_seed": v_seed,
    }


def _time_call(func, repeats: int) -> float:
    start = perf_counter()
    for _ in range(repeats):
        func()
    return (perf_counter() - start) / repeats


def _profile_case(num_threads: int) -> dict[str, float]:
    shell = _build_shell_state()
    get_y = electron_module.get_y
    timings = {
        "get_syn_simpson": _time_call(
            lambda: get_y.get_syn_simpson(
                shell["r_loc"], shell["db"], num_threads, shell["gam_e"], shell["d_n_gam_e"], shell["v_seed"]
            ),
            PROFILE_REPEATS,
        ),
        "get_nu_a": _time_call(
            lambda: get_y.get_nu_a(shell["r_loc"], shell["db"], shell["gam_e"], shell["d_n_gam_e"]),
            PROFILE_REPEATS,
        ),
        "get_ssa_numerical": _time_call(
            lambda: get_y.get_ssa_numerical(
                shell["db"], num_threads, shell["gam_e"], shell["v_seed"], shell["seed_syn"]
            ),
            PROFILE_REPEATS,
        ),
        "get_y_nakar": _time_call(
            lambda: get_y.get_y_nakar(num_threads, shell["gam_e"], shell["v_seed"], shell["p_syn"]),
            PROFILE_REPEATS,
        ),
        "get_ic_numerical": _time_call(
            lambda: get_y.get_ic_numerical(num_threads, shell["gam_e"], shell["v_seed"], shell["seed_syn"]),
            PROFILE_REPEATS,
        ),
    }
    return {"num_threads": num_threads, "timings": timings}


def _plot(payload: dict[str, dict[str, object]], output_png: Path, output_pdf: Path) -> None:
    serial = payload["serial"]
    parallel = payload["parallel"]
    serial_timings = serial["timings"]
    parallel_timings = parallel["timings"]
    labels = [
        label
        for label in sorted(
            serial_timings,
            key=lambda name: max(serial_timings[name], parallel_timings[name]),
            reverse=True,
        )
    ]
    y = np.arange(len(labels), dtype=float)
    serial_values = np.array([serial_timings[label] for label in labels], dtype=float)
    parallel_values = np.array([parallel_timings[label] for label in labels], dtype=float)

    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    ax.barh(y - 0.2, serial_values, height=0.38, label=f"serial ({serial['num_threads']} thread)")
    ax.barh(y + 0.2, parallel_values, height=0.38, label=f"parallel ({parallel['num_threads']} threads)")
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("seconds per call")
    ax.set_title("ASGARD electron kernel micro-profile")
    ax.grid(axis="x", alpha=0.3)
    max_value = max(float(serial_values.max(initial=0.0)), float(parallel_values.max(initial=0.0)))
    text_dx = 0.01 * max(max_value, 1.0)
    for yi, s_val, p_val in zip(y, serial_values, parallel_values):
        if p_val > 0.0:
            ax.text(p_val + text_dx, yi + 0.2, f"x{s_val / p_val:.2f}", va="center", fontsize=8)
    ax.legend(loc="lower right")
    fig.savefig(output_png, dpi=200)
    fig.savefig(output_pdf)
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    payload = {
        "serial": _profile_case(SERIAL_THREADS),
        "parallel": _profile_case(PARALLEL_THREADS),
    }

    output_json = OUTPUT_DIR / "electron_kernel_profile.json"
    output_png = OUTPUT_DIR / "electron_kernel_profile.png"
    output_pdf = OUTPUT_DIR / "electron_kernel_profile.pdf"

    output_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    _plot(payload, output_png, output_pdf)

    print(output_json)
    print(output_png)
    print(output_pdf)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
