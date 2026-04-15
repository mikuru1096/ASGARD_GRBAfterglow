from __future__ import annotations

from pathlib import Path
import sys
import warnings

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mergered import fit


def main() -> None:
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        redchi = fit(
            Num_threads=8,
            index_dyn=3,
            index_Y=2,
            index_syn_intger=2,
            Num_gam_e=101,
            Num_nu=101,
            Num_R=120,
            Num_theta=180,
            Num_phi=1,
            z=0.4,
            Eta_0=1.0e2,
            Epsilon_e=1.0e-1,
            Epsilon_b=1.0e-3,
            p=2.5,
            OpeningAngle_jet=1.0e-1,
            theta_v=0.0,
            f_e=1.0e-1,
            E_iso=1.0e53,
            dNe=1.0e-1,
            A_star=-1.0,
            R0=1.0e9,
            Ebv=0.0,
            Rv=2.93,
            Lyman_Ar=0.0,
            f_sys=-1.0,
            plot_LC=False,
        )

    assert np.isfinite(redchi)
    assert any(issubclass(item.category, DeprecationWarning) for item in caught)
    print("PASS: legacy fit(**kwargs) check succeeded.")


if __name__ == "__main__":
    main()
