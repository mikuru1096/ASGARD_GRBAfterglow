from __future__ import annotations

from _repo_path import ensure_repo_root_on_path


ensure_repo_root_on_path()

from asgard_core import Fitter, Param, Scale
from tests.public_api_builders import top_hat_model


def main() -> None:
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
