from __future__ import annotations

from pathlib import Path
import sys
from unittest.mock import patch

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import asgard_component_backend as backend
from asgard_presets import build_baseline_config
from asgard_setup import build_simulation_setup


def main() -> None:
    config = build_baseline_config(
        num_threads=1,
        num_gam_e=32,
        num_nu=81,
        num_r=80,
        num_tobs=48,
        epsilon_e=0.2,
        epsilon_b=1.0e-5,
        d_ne=10.0,
        a_star=-1.0,
        electron_solver="slc1_mmg2",
        include_forward_ssc=True,
    )
    setup = build_simulation_setup(config)

    aux_calls = {"count": 0}
    original_aux = backend.compute_ssc_auxiliary_grid

    def _wrapped_aux(*args, **kwargs):
        aux_calls["count"] += 1
        return original_aux(*args, **kwargs)

    def _forbidden_nonuniform(*_args, **_kwargs):
        raise AssertionError("High-level slc1_mmg2 SSC path should not call Radiation.ssc_spec_nonuniform.")

    with patch.object(backend, "compute_ssc_auxiliary_grid", side_effect=_wrapped_aux):
        with patch.object(backend.Radiation, "ssc_spec_nonuniform", side_effect=_forbidden_nonuniform):
            component_spectra = backend.solve_component_spectra_from_setup(config, setup)

    assert aux_calls["count"] >= 1, aux_calls["count"]
    assert np.all(np.isfinite(component_spectra.fwd_ssc))
    assert np.any(component_spectra.fwd_ssc > 0.0)
    print("PASS: slc1_mmg2 SSC route uses auxiliary grid.")


if __name__ == "__main__":
    main()
