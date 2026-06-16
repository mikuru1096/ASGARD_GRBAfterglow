from __future__ import annotations

from _repo_path import ensure_repo_root_on_path
import importlib.util


ensure_repo_root_on_path()

from asgard_core import Fitter, Param, Scale
from tests.public_api_builders import top_hat_model


class DummyFitter(Fitter):
    def loglike(self, values: dict[str, float]) -> float:
        if "_compiled_problem" in self.__dict__:
            raise AssertionError("Fitter.fit left a stale cached _compiled_problem entry.")
        x = float(values["p"])
        return -(x - 2.3) ** 2


def main() -> None:
    if importlib.util.find_spec("emcee") is None:
        print("fitter_public_api_smoke: skipped, emcee is not installed")
        return

    fitter = DummyFitter(model=top_hat_model())
    fitter.__dict__["_compiled_problem"] = object()
    result = fitter.fit(
        param_defs=[Param("p", "fwd_rad.p", 2.1, 2.5, Scale.LINEAR)],
        sampler="emcee",
        total_steps=4,
        burn_frac=0.0,
        nwalkers=8,
    )
    assert result.labels == ["p"]
    assert "p" in result.best_params
    assert result.samples is not None

    fixed_fitter = DummyFitter(model=top_hat_model())
    fixed_fitter.__dict__["_compiled_problem"] = object()
    fixed = fixed_fitter.fit(
        param_defs=[Param("p", "fwd_rad.p", 2.3, 2.3, Scale.FIXED)],
        sampler="emcee",
    )
    assert fixed.best_params == {"p": 2.3}
    assert fixed.labels == ["p"]
    print("fitter_public_api_smoke: ok")


if __name__ == "__main__":
    main()
