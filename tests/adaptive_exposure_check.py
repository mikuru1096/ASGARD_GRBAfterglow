from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet


def _build_model() -> Model:
    return Model(
        medium=ISM(0.1),
        jet=TophatJet(E_iso=1.0e53, Gamma0=100.0, theta_c=0.1),
        observer=Observer(lumi_dist=1.0e28, z=0.4, theta_obs=0.0, phi_obs=0.0),
        fwd_rad=Radiation(eps_e=0.1, eps_B=1.0e-3, p=2.5, xi_N=0.1, ssc=True, kn=False),
        setups=Setups(
            num_threads=8,
            num_gam_e=121,
            num_nu=121,
            num_r=140,
            num_theta=120,
            num_phi=1,
            num_tobs=64,
            observer_time_min_s=1.0e2,
            observer_time_max_s=1.0e6,
        ),
    )


def main() -> None:
    model = _build_model()
    times_s = np.logspace(2.2, 4.8, 8)
    frequencies_hz = np.full(times_s.shape, 1.0e14)
    exposures_s = 0.3 * times_s

    adaptive = model.flux_density_exposures(times_s, frequencies_hz, exposures_s, num_subsamples=4).total

    dense_times = []
    dense_freqs = []
    spans = []
    start = 0
    for time_s, exposure_s, frequency_hz in zip(times_s, exposures_s, frequencies_hz):
        t1 = max(time_s - 0.5 * exposure_s, 1.0e-30)
        t2 = time_s + 0.5 * exposure_s
        nodes = np.linspace(t1, t2, 129, dtype=float)
        dense_times.append(nodes)
        dense_freqs.append(np.full(nodes.shape, frequency_hz, dtype=float))
        stop = start + nodes.shape[0]
        spans.append((start, stop, t1, t2))
        start = stop

    dense_sample = model.flux_density(np.concatenate(dense_times), np.concatenate(dense_freqs)).total
    reference = np.zeros_like(adaptive)
    for idx, (start, stop, t1, t2) in enumerate(spans):
        nodes = dense_times[idx]
        reference[idx] = float(np.trapezoid(dense_sample[start:stop], nodes) / (t2 - t1))

    rel = np.max(np.abs(adaptive - reference) / np.maximum(reference, 1.0e-99))
    assert rel <= 1.0e-2, rel
    print("PASS: adaptive exposure check succeeded.")


if __name__ == "__main__":
    main()
