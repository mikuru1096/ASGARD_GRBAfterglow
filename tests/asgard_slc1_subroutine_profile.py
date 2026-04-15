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

from asgard_component_backend import build_query_setup, profile_observe_spectra_from_setup
from asgard_models import default_num_threads
from VegasAfterglow import ISM, Model, Observer, Radiation, Setups, TophatJet
from VegasAfterglow.api import _build_fit_config_for_patch


OUTPUT_DIR = ROOT / "output" / "vegasafterglow_doc"
TIMES_S = np.logspace(2.0, 8.0, 100)
FREQS_HZ = np.array([1.0e9, 1.0e14, 1.0e18], dtype=float)
PROFILE_REPEATS = 2
NUM_GAM_BY_SOLVER = {"fullhide": 121, "slc1": 32, "slc1_mmg2": 32}


def _build_model(solver: str) -> Model:
    return Model(
        TophatJet(theta_c=0.1, E_iso=1.0e52, Gamma0=300.0),
        ISM(1.0),
        Observer(lumi_dist=1.0e26, z=0.1, theta_obs=0.0),
        Radiation(0.1, 1.0e-3, 2.3, ssc=False, kn=False),
        setups=Setups(
            electron_solver=solver,
            num_threads=default_num_threads(),
            num_gam_e=NUM_GAM_BY_SOLVER[solver],
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


def _profile_solver(solver: str) -> dict[str, object]:
    model = _build_model(solver)
    config = _build_direct_config(model)
    setup = build_query_setup(config, TIMES_S, FREQS_HZ)
    runs: list[dict[str, float]] = []
    totals: list[float] = []

    for _ in range(PROFILE_REPEATS):
        start = perf_counter()
        _, observed, timings = profile_observe_spectra_from_setup(config, setup, FREQS_HZ)
        totals.append(perf_counter() - start)
        runs.append(timings)
        assert observed["total"].shape == (FREQS_HZ.size, TIMES_S.size)
        assert np.all(np.isfinite(observed["total"]))

    timings = _median_timings(runs)
    timings["Total pipeline"] = float(np.median(totals))
    return {"timings": timings}


def _plot(compare: dict[str, dict[str, object]]) -> None:
    fullhide_timings = compare["fullhide"]["timings"]
    slc1_timings = compare["slc1"]["timings"]
    mmg2_timings = compare["slc1_mmg2"]["timings"]
    labels = [
        label
        for label in sorted(
            set(fullhide_timings) | set(slc1_timings) | set(mmg2_timings),
            key=lambda name: max(fullhide_timings.get(name, 0.0), slc1_timings.get(name, 0.0), mmg2_timings.get(name, 0.0)),
            reverse=True,
        )
        if max(fullhide_timings.get(label, 0.0), slc1_timings.get(label, 0.0), mmg2_timings.get(label, 0.0)) > 0.0
    ]
    y = np.arange(len(labels), dtype=float)
    v0 = np.array([fullhide_timings.get(label, 0.0) for label in labels], dtype=float)
    v1 = np.array([slc1_timings.get(label, 0.0) for label in labels], dtype=float)
    v2 = np.array([mmg2_timings.get(label, 0.0) for label in labels], dtype=float)

    fig, ax = plt.subplots(figsize=(12, 8), constrained_layout=True)
    ax.barh(y - 0.26, v0, height=0.24, label="fullhide")
    ax.barh(y, v1, height=0.24, label="slc1")
    ax.barh(y + 0.26, v2, height=0.24, label="slc1_mmg2")
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("seconds")
    ax.set_title("ASGARD subroutine timings: fullhide vs slc1 vs slc1_mmg2")
    ax.grid(axis="x", alpha=0.3)
    ax.legend(loc="lower right")

    max_value = max(float(v0.max(initial=0.0)), float(v1.max(initial=0.0)), float(v2.max(initial=0.0)))
    dx = 0.01 * max(max_value, 1.0)
    for yi, a, b in zip(y, v0, v1):
        if b > 0.0:
            ax.text(b + dx, yi, f"x{a / b:.2f}", va="center", fontsize=8)
    for yi, a, b in zip(y, v0, v2):
        if b > 0.0:
            ax.text(b + dx, yi + 0.26, f"x{a / b:.2f}", va="center", fontsize=8)

    fig.savefig(OUTPUT_DIR / "slc1_subroutine_profile.png", dpi=200)
    fig.savefig(OUTPUT_DIR / "slc1_subroutine_profile.pdf")
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    compare = {
        "fullhide": _profile_solver("fullhide"),
        "slc1": _profile_solver("slc1"),
        "slc1_mmg2": _profile_solver("slc1_mmg2"),
    }
    output_json = OUTPUT_DIR / "slc1_subroutine_profile.json"
    with output_json.open("w", encoding="utf-8") as fh:
        json.dump(compare, fh, indent=2)
    _plot(compare)
    print(output_json)
    print(OUTPUT_DIR / "slc1_subroutine_profile.png")
    print(OUTPUT_DIR / "slc1_subroutine_profile.pdf")
    print(json.dumps(compare, indent=2))


if __name__ == "__main__":
    main()
