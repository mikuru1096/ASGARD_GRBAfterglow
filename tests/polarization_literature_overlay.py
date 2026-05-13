from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet
from asgard_core.asgard_paths import ASGARD_DOC_DIR, ensure_output_dir


OUTPUT_DIR = ensure_output_dir(ASGARD_DOC_DIR)
OUTPUT_PNG = OUTPUT_DIR / "polarization_lan2023_overlay.png"
OUTPUT_CSV = OUTPUT_DIR / "polarization_lan2023_overlay_data.csv"


# Lan, Wu & Dai 2023, ApJ 952, 31, Fig. 6 lower-right PD panel.
# Case II is forward-shock dominated; this curve is q=2.0, k=0, optical R-band.
LAN2023_CASEII_Q20_K0 = np.array(
    [
        [2.190623e00, 8.163265e-01],
        [4.873186e00, 6.122449e-01],
        [1.117928e01, 6.122449e-01],
        [3.084243e01, 6.122449e-01],
        [6.756413e01, 6.122449e-01],
        [1.480075e02, 6.122449e-01],
        [4.342397e02, 8.163265e-01],
        [9.367406e02, 1.020408e00],
        [2.052045e03, 1.632653e00],
        [4.036547e03, 4.285714e00],
        [8.443923e03, 1.918367e01],
        [1.298741e04, 3.428571e01],
        [2.515752e04, 4.612245e01],
        [4.873186e04, 4.530612e01],
        [1.225973e05, 4.346939e01],
        [2.486906e05, 3.530612e01],
        [4.194735e05, 2.428571e01],
        [8.001510e05, 1.428571e01],
        [1.699746e06, 8.367347e00],
        [3.666687e06, 5.510204e00],
        [1.011599e07, 3.469388e00],
        [2.052045e07, 2.857143e00],
        [5.162437e07, 2.448980e00],
        [7.352653e07, 2.244898e00],
    ],
    dtype=float,
)


def _build_asgard_model() -> Model:
    """构建 Lan 2023 Case II k=0、q=2.0 的 forward-shock top-hat 偏振模型。"""
    theta_j = 0.1
    return Model(
        TophatJet(theta_j, 1.0e50, 100.0),
        ISM(1.0),
        Observer(1.0e28, 0.3, 2.0 * theta_j),
        Radiation(0.1, 0.1, 2.5),
        setups=Setups(
            electron_solver="fullhide_1d",
            num_threads=1,
            num_gam_e=32,
            num_nu=40,
            num_r=40,
            num_phi=48,
            num_tobs=40,
            patch_theta=8,
            patch_phi=48,
        ),
    )


def _asgard_polarization_curve() -> tuple[np.ndarray, np.ndarray]:
    """计算 ASGARD 对应模型的光学 R-band 偏振曲线。"""
    model = _build_asgard_model()
    times = np.logspace(0.0, 8.0, 32)
    pol = model.polarization(times, np.array([4.56e14]), magnetic_geometry="shock_random")
    return times, 100.0 * pol.linear_polarization[0]


def _write_overlay_csv(
    lan_t: np.ndarray,
    lan_p: np.ndarray,
    asgard_x: np.ndarray,
    asgard_p: np.ndarray,
) -> None:
    """保存 Lan 2023 数字化点和 ASGARD 曲线，便于后续复查。"""
    rows = ["series,time_s,polarization_percent"]
    rows.extend(f"Lan2023_Fig6_CaseII_q2_k0,{x:.8e},{p:.8e}" for x, p in zip(lan_t, lan_p, strict=True))
    rows.extend(f"ASGARD_CaseII_q2_k0,{x:.8e},{p:.8e}" for x, p in zip(asgard_x, asgard_p, strict=True))
    OUTPUT_CSV.write_text("\n".join(rows) + "\n", encoding="utf-8")


def main() -> None:
    """生成 Lan 2023 数字化曲线与 ASGARD 当前输出的叠图。"""
    lan_t = LAN2023_CASEII_Q20_K0[:, 0]
    lan_p = LAN2023_CASEII_Q20_K0[:, 1]
    asgard_x, asgard_p = _asgard_polarization_curve()
    _write_overlay_csv(lan_t, lan_p, asgard_x, asgard_p)

    fig, ax = plt.subplots(figsize=(6.8, 4.4), constrained_layout=True)
    ax.plot(lan_t, lan_p, "o-", color="black", lw=1.7, ms=4.0, label="Lan et al. 2023 Fig. 6, Case II q=2.0 k=0")
    ax.plot(asgard_x, asgard_p, "s-", color="#1f77b4", lw=1.7, ms=3.8, label="ASGARD FS shock_random, ISM q=2.0")
    ax.set_xscale("log")
    ax.set_xlim(1.0, 1.0e8)
    ax.set_ylim(0.0, 1.1 * max(float(np.max(lan_p)), float(np.max(asgard_p))))
    ax.set_xlabel("Observer time (s)")
    ax.set_ylabel("Linear polarization (%)")
    ax.set_title("Lan et al. 2023 off-axis forward-shock polarization")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="upper right", fontsize=8)
    fig.savefig(OUTPUT_PNG, dpi=180)
    print(f"wrote {OUTPUT_PNG}")
    print(f"wrote {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
