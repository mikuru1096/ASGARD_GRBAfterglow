from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
import json
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet
from tests.bug_audit_check import (
    check_component_finiteness_and_closure,
    check_monotonic_arrays,
    check_nu_a_bounds,
    check_thread_consistency,
)
from tests.physical_closure_check import (
    check_absorption_emission_closure,
    check_dynamics_radiation_closure,
    check_forward_reverse_closure,
    check_observer_frame_closure,
    check_regime_ordering_closure,
    check_spectral_slope_closure,
)


OUTPUT_DIR = ROOT / "output" / "asgard_doc"
JSON_OUT = ROOT / "tests" / "baseline_architecture.json"
NPZ_OUT = ROOT / "tests" / "baseline_architecture_observed.npz"


def _summary(items) -> dict[str, int]:
    total = len(items)
    failed = sum(0 if item.passed else 1 for item in items)
    return {"total": total, "failed": failed}


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
            num_tobs=48,
            observer_time_min_s=1.0e2,
            observer_time_max_s=1.0e6,
        ),
    )


def main() -> None:
    bug_checks = (
        check_thread_consistency()
        + check_nu_a_bounds()
        + check_monotonic_arrays()
        + check_component_finiteness_and_closure()
    )
    closure_checks = (
        check_regime_ordering_closure()
        + check_spectral_slope_closure()
        + check_dynamics_radiation_closure()
        + check_absorption_emission_closure()
        + check_forward_reverse_closure()
        + check_observer_frame_closure()
    )

    comprehensive_path = OUTPUT_DIR / "comprehensive_validation_asgard.json"
    subroutine_profile_path = OUTPUT_DIR / "subroutine_profile_tophat.json"
    electron_profile_path = OUTPUT_DIR / "electron_kernel_profile.json"

    payload = {
        "regression_baseline_npz": str((ROOT / "tests" / "baseline_lightcurves.npz").resolve()),
        "comprehensive_validation": json.loads(comprehensive_path.read_text(encoding="utf-8")),
        "bug_audit": {
            "summary": _summary(bug_checks),
            "checks": [asdict(item) for item in bug_checks],
        },
        "physical_closure": {
            "summary": _summary(closure_checks),
            "checks": [asdict(item) for item in closure_checks],
        },
        "subroutine_profile": json.loads(subroutine_profile_path.read_text(encoding="utf-8")),
        "electron_kernel_profile": json.loads(electron_profile_path.read_text(encoding="utf-8")),
    }

    times_s = np.logspace(2.0, 5.0, 12)
    frequencies_hz = np.array([9.0e9, 4.84e14, 1.0e18], dtype=float)
    model = _build_model()
    flux = model.flux_density_grid(times_s, frequencies_hz)
    details = model.details(1.0e2, 1.0e5)

    np.savez(
        NPZ_OUT,
        times_s=times_s,
        frequencies_hz=frequencies_hz,
        total=flux.total,
        fwd_sync=flux.fwd.sync,
        fwd_ssc=flux.fwd.ssc,
        rev_sync=flux.rev.sync,
        rev_ssc=flux.rev.ssc,
        cross_ic=np.zeros_like(flux.total) if flux.cross_ic is None else flux.cross_ic,
        detail_t_obs=details.fwd.t_obs,
        detail_nu_m=details.fwd.nu_m,
        detail_nu_c=details.fwd.nu_c,
        detail_nu_a=details.fwd.nu_a,
    )
    JSON_OUT.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print(JSON_OUT)
    print(NPZ_OUT)


if __name__ == "__main__":
    main()
