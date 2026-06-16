from __future__ import annotations

from _repo_path import ensure_repo_root_on_path


ensure_repo_root_on_path()

from asgard_core.asgard_config import RuntimeConfig
from asgard_core.asgard_state import _validate_joint_electron_photon_config


def test_default_coupling_is_separated() -> None:
    config = RuntimeConfig()
    assert config.electron_photon_coupling == "separated"


def test_joint_requires_bethe_heitler_boundary() -> None:
    config = RuntimeConfig(electron_photon_coupling="joint")
    try:
        _validate_joint_electron_photon_config(config)
    except ValueError as exc:
        assert "Bethe-Heitler" in str(exc)
    else:
        raise AssertionError("joint coupling without Bethe-Heitler must fail at the configuration boundary.")


def test_joint_accepts_gamma_gamma_feedback() -> None:
    config = RuntimeConfig(electron_photon_coupling="joint")
    config.hadronic.enabled = True
    config.hadronic.solver = "am3_1d"
    config.hadronic.epsilon_p = 0.1
    config.hadronic.include_bethe_heitler = True
    config.hadronic.include_pair_production = True
    config.hadronic.pgamma_scheme = "hummer_2010_response"
    config.index_y = 1
    _validate_joint_electron_photon_config(config)


def test_joint_requires_numeric_ic_cooling() -> None:
    config = RuntimeConfig(electron_photon_coupling="joint")
    config.hadronic.enabled = True
    config.hadronic.solver = "am3_1d"
    config.hadronic.epsilon_p = 0.1
    config.hadronic.include_bethe_heitler = True
    config.hadronic.pgamma_scheme = "hummer_2010_response"
    config.index_y = 2
    try:
        _validate_joint_electron_photon_config(config)
    except ValueError as exc:
        assert "ssc_cooling_mode='numeric_ic_kn'" in str(exc)
    else:
        raise AssertionError("joint coupling must reject non-numeric IC cooling.")


def main() -> None:
    test_default_coupling_is_separated()
    test_joint_requires_bethe_heitler_boundary()
    test_joint_accepts_gamma_gamma_feedback()
    test_joint_requires_numeric_ic_cooling()
    print("electron_photon_coupling_config_smoke: ok")


if __name__ == "__main__":
    main()
