from __future__ import annotations

from _repo_path import ensure_repo_root_on_path

import numpy as np


ensure_repo_root_on_path()

from asgard_core import Fitter, Param, Scale
from tests.public_api_builders import solver_options, top_hat_model


def _check_nu_callback() -> None:
    calls = []

    def collect(label, nu_m, nu_c, nu_a):
        calls.append((label, nu_m, nu_c, nu_a))

    model = top_hat_model(solver_options=solver_options(nu_callback=collect))
    model.flux_density_grid(np.array([1.0e3]), np.array([1.0e9]))
    assert calls
    label, nu_m, nu_c, nu_a = calls[0]
    assert label == "fullhide_1d"
    assert nu_m.shape == nu_c.shape == nu_a.shape
    assert np.all(np.isfinite(nu_m))
    assert not hasattr(model.details().fwd, "nu_m")


def main() -> None:
    _check_nu_callback()
    fitter = Fitter(model=top_hat_model())
    fixed = fitter.fit(
        param_defs=[Param("p", "fwd_rad.p", 2.3, 2.3, Scale.FIXED)],
        sampler="emcee",
    )
    assert fixed.best_params == {"p": 2.3}
    assert fixed.labels == ["p"]
    print("fitter_public_api_smoke: ok")


if __name__ == "__main__":
    main()
