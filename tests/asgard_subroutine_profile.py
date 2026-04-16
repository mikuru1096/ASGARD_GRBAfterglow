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

from asgard_paths import ASGARD_DOC_DIR
from asgard_component_backend import build_query_setup, profile_observe_spectra_from_setup
from asgard_models import default_num_threads
from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet
from ASGARD.api import _build_fit_config_for_patch


OUTPUT_DIR = ASGARD_DOC_DIR
TIMES_S = np.logspace(2.0, 8.0, 100)
FREQS_HZ = np.array([1.0e9, 1.0e14, 1.0e18], dtype=float)
SERIAL_THREADS = 1
PARALLEL_THREADS = default_num_threads()
PROFILE_REPEATS = 2


def _build_model(*, ssc: bool, num_threads: int) -> Model:
    return Model(
        TophatJet(theta_c=0.1, E_iso=1.0e52, Gamma0=300.0),
        ISM(1.0),
        Observer(lumi_dist=1.0e26, z=0.1, theta_obs=0.0),
        Radiation(0.1, 1.0e-3, 2.3, ssc=ssc, kn=False),
        setups=Setups(
            num_threads=num_threads,
            num_gam_e=201,
            num_nu=201,
            num_r=300,
            num_theta=300,
            num_phi=1,
            num_tobs=200,
        ),
    )


def _build_direct_config(model: Model):
    return _build_fit_config_for_patch(
        model,
        phi_center=0.0,
        theta_v=model.observer.theta_obs,
        opening_angle_jet=model.jet.theta_j,
        e_iso=model.jet.E_iso,
        gamma0=model.jet.lf,
        theta_center=0.0,
    )


def _median_timings(records: list[dict[str, float]]) -> dict[str, float]:
    labels = sorted({label for record in records for label in record})
    result: dict[str, float] = {}
    for label in labels:
        values = [record.get(label, 0.0) for record in records]
        result[label] = float(np.median(values))
    return result


def _profile_case(*, ssc: bool, num_threads: int) -> dict[str, object]:
    model = _build_model(ssc=ssc, num_threads=num_threads)
    config = _build_direct_config(model)
    setup = build_query_setup(config, TIMES_S, FREQS_HZ)

    timings_runs: list[dict[str, float]] = []
    total_runs: list[float] = []

    for _ in range(PROFILE_REPEATS):
        start = perf_counter()
        _, observed, timings = profile_observe_spectra_from_setup(config, setup, FREQS_HZ)
        total_runs.append(perf_counter() - start)
        timings_runs.append(timings)
        assert observed["total"].shape == (FREQS_HZ.size, TIMES_S.size)
        assert np.all(np.isfinite(observed["total"]))

    timings = _median_timings(timings_runs)
    timings["Total pipeline"] = float(np.median(total_runs))
    return {
        "num_threads": num_threads,
        "timings": timings,
    }


def _plot_case(ax, title: str, serial: dict[str, float], parallel: dict[str, float]) -> None:
    serial_timings = serial["timings"]
    parallel_timings = parallel["timings"]
    labels = [
        label
        for label in sorted(
            set(serial_timings) | set(parallel_timings),
            key=lambda name: max(serial_timings.get(name, 0.0), parallel_timings.get(name, 0.0)),
            reverse=True,
        )
        if max(serial_timings.get(label, 0.0), parallel_timings.get(label, 0.0)) > 0.0
    ]
    y = np.arange(len(labels), dtype=float)
    serial_values = np.array([serial_timings.get(label, 0.0) for label in labels], dtype=float)
    parallel_values = np.array([parallel_timings.get(label, 0.0) for label in labels], dtype=float)

    ax.barh(y - 0.2, serial_values, height=0.38, label=f"serial ({serial['num_threads']} thread)")
    ax.barh(y + 0.2, parallel_values, height=0.38, label=f"parallel ({parallel['num_threads']} threads)")
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("seconds")
    ax.set_title(title)
    ax.grid(axis="x", alpha=0.3)

    max_value = max(float(serial_values.max(initial=0.0)), float(parallel_values.max(initial=0.0)))
    text_dx = 0.01 * max(max_value, 1.0)
    for yi, s_val, p_val in zip(y, serial_values, parallel_values):
        if p_val > 0.0:
            speedup = s_val / p_val if s_val > 0.0 else 0.0
            ax.text(p_val + text_dx, yi + 0.2, f"x{speedup:.2f}", va="center", fontsize=8)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    cases = {
        "Tophat sync": {
            "serial": _profile_case(ssc=False, num_threads=SERIAL_THREADS),
            "parallel": _profile_case(ssc=False, num_threads=PARALLEL_THREADS),
        },
        "Tophat full SSC": {
            "serial": _profile_case(ssc=True, num_threads=SERIAL_THREADS),
            "parallel": _profile_case(ssc=True, num_threads=PARALLEL_THREADS),
        },
    }

    output_json = OUTPUT_DIR / "subroutine_profile_tophat.json"
    output_png = OUTPUT_DIR / "subroutine_profile_tophat.png"
    output_pdf = OUTPUT_DIR / "subroutine_profile_tophat.pdf"

    with output_json.open("w", encoding="utf-8") as f:
        json.dump(cases, f, indent=2)

    fig, axes = plt.subplots(2, 1, figsize=(12, 10), constrained_layout=True)
    for ax, (title, payload) in zip(np.atleast_1d(axes), cases.items()):
        _plot_case(ax, title, payload["serial"], payload["parallel"])
    axes[0].legend(loc="lower right")
    fig.suptitle("ASGARD / ASGARD-style subroutine timings")
    fig.savefig(output_png, dpi=200)
    fig.savefig(output_pdf)
    plt.close(fig)

    print(output_json)
    print(output_png)
    print(output_pdf)
    print(json.dumps(cases, indent=2))


if __name__ == "__main__":
    main()
