from __future__ import annotations

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet


def main() -> None:
    model = Model(
        TophatJet(0.1, 1.0e52, 300.0),
        ISM(1.0),
        Observer(1.0e26, 0.1, 0.0),
        Radiation(0.1, 1.0e-3, 2.3, epsilon_p=0.2, proton_synch=True, pg=True, neutrino=True),
        setups=Setups(
            electron_solver="fullhide_1d",
            hadronic_enabled=True,
            hadronic_solver="am3_1d",
            pgamma_scheme="ka2008_reference",
            num_gam_e=24,
            num_gam_p=40,
            num_nu=40,
            num_nu_nu=24,
            num_r=24,
            num_theta=16,
            num_tobs=24,
        ),
    )
    try:
        model.details(1.0e3, 1.0e6)
    except ValueError as exc:
        message = str(exc)
        assert "tabulated" in message
        assert "eta/eta0" in message
        return
    raise AssertionError("Expected KA2008 domain guard for pgamma_scheme='ka2008_reference'.")


if __name__ == "__main__":
    main()
