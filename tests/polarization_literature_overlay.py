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
OUTPUT_PNG = OUTPUT_DIR / "polarization_gl99_overlay.png"
OUTPUT_CSV = OUTPUT_DIR / "polarization_gl99_overlay_data.csv"


# Ghisellini & Lazzati 1999, MNRAS 309, L7, Fig. 4 lower panel.
# Curve: theta_o / theta_c = 0.9, P0 = 60 per cent; x-axis is Gamma^{-1}.
GL99_THETA_RATIO_09 = np.array(
    [
        [0.010, 2.0],
        [0.013, 3.0],
        [0.016, 3.7],
        [0.020, 4.0],
        [0.026, 3.3],
        [0.035, 2.0],
        [0.046, 0.4],
        [0.055, 0.0],
        [0.070, 1.4],
        [0.085, 3.6],
        [0.105, 6.3],
        [0.125, 8.2],
        [0.150, 9.4],
        [0.180, 8.8],
        [0.220, 7.0],
        [0.280, 4.7],
        [0.350, 3.0],
        [0.450, 1.8],
    ],
    dtype=float,
)


def _build_asgard_model() -> Model:
    """构建与 GL99 Fig. 4 曲线相同视角比的 top-hat 偏振模型。"""
    theta_j = 0.1
    return Model(
        TophatJet(theta_j, 1.0e52, 300.0),
        ISM(1.0),
        Observer(1.0e26, 0.1, 0.9 * theta_j),
        Radiation(0.1, 1.0e-3, 2.3),
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
    """计算 ASGARD 偏振曲线，并按文献 P0=60% 归一化本征偏振率。"""
    model = _build_asgard_model()
    times = np.logspace(2.0, 8.0, 28)
    pol = model.polarization(times, np.array([1.0e14]), magnetic_geometry="shock_random")
    details = model.details(float(times[0]), float(times[-1]))
    gamma = np.interp(np.log(times), np.log(details.fwd.t_obs), details.fwd.Gamma)
    intrinsic_pi = (2.3 + 1.0) / (2.3 + 7.0 / 3.0)
    scaled_percent = 100.0 * pol.linear_polarization[0] * (0.60 / intrinsic_pi)
    return 1.0 / gamma, scaled_percent


def _unit_direction(theta: float, phi: float) -> np.ndarray:
    """球坐标方向向量，供 GL99 几何积分使用。"""
    return np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)], dtype=float)


def _gl99_sky_basis(theta_obs: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """构造视线和天球基底，使 x 轴位于 jet-axis/observer 平面。"""
    sightline = _unit_direction(theta_obs, 0.0)
    sky_y = np.array([0.0, 1.0, 0.0], dtype=float)
    sky_x = np.cross(sky_y, sightline)
    sky_x = sky_x / np.linalg.norm(sky_x)
    sky_y = np.cross(sightline, sky_x)
    return sightline, sky_x, sky_y / np.linalg.norm(sky_y)


def _gl99_geometry_curve() -> tuple[np.ndarray, np.ndarray]:
    """按 GL99 的薄壳几何假设直接积分 top-hat 面元偏振。"""
    theta_j = 0.1
    theta_obs = 0.9 * theta_j
    alpha = 0.6
    sightline, sky_x, sky_y = _gl99_sky_basis(theta_obs)
    theta_edges = np.linspace(0.0, theta_j, 121)
    phi_edges = np.linspace(0.0, 2.0 * np.pi, 289)
    gamma_inverse = np.logspace(np.log10(0.01), np.log10(0.45), 48)
    polarization_percent = np.zeros_like(gamma_inverse)

    for i_gamma, gamma_inv in enumerate(gamma_inverse):
        gamma = 1.0 / gamma_inv
        beta = np.sqrt(1.0 - gamma ** -2)
        total_i = 0.0
        total_q = 0.0
        total_u = 0.0
        for i_theta in range(theta_edges.size - 1):
            theta = 0.5 * (theta_edges[i_theta] + theta_edges[i_theta + 1])
            domega_theta = np.cos(theta_edges[i_theta]) - np.cos(theta_edges[i_theta + 1])
            for i_phi in range(phi_edges.size - 1):
                phi = 0.5 * (phi_edges[i_phi] + phi_edges[i_phi + 1])
                patch_direction = _unit_direction(theta, phi)
                mu = float(np.dot(patch_direction, sightline))
                doppler = 1.0 / (gamma * (1.0 - beta * mu))
                intensity_weight = doppler ** (3.0 + alpha) * domega_theta * (phi_edges[i_phi + 1] - phi_edges[i_phi])
                e_vector = patch_direction - mu * sightline
                e_x = float(np.dot(e_vector, sky_x))
                e_y = float(np.dot(e_vector, sky_y))
                norm = np.hypot(e_x, e_y)
                e_x = e_x / norm
                e_y = e_y / norm
                mu_prime = (mu - beta) / (1.0 - beta * mu)
                local_pi = (1.0 - mu_prime * mu_prime) / (1.0 + mu_prime * mu_prime)
                total_i += intensity_weight
                total_q += intensity_weight * local_pi * (e_x * e_x - e_y * e_y)
                total_u += intensity_weight * local_pi * (2.0 * e_x * e_y)
        polarization_percent[i_gamma] = 60.0 * np.hypot(total_q, total_u) / total_i
    return gamma_inverse, polarization_percent


def _write_overlay_csv(
    gl_x: np.ndarray,
    gl_p: np.ndarray,
    geom_x: np.ndarray,
    geom_p: np.ndarray,
    asgard_x: np.ndarray,
    asgard_p: np.ndarray,
) -> None:
    """保存文献数字化点和 ASGARD 曲线，便于后续复查。"""
    rows = ["series,gamma_inverse,polarization_percent"]
    rows.extend(f"GL99_theta_ratio_0.9,{x:.8e},{p:.8e}" for x, p in zip(gl_x, gl_p, strict=True))
    rows.extend(f"GL99_geometry_reproduction,{x:.8e},{p:.8e}" for x, p in zip(geom_x, geom_p, strict=True))
    rows.extend(f"ASGARD_theta_ratio_0.9,{x:.8e},{p:.8e}" for x, p in zip(asgard_x, asgard_p, strict=True))
    OUTPUT_CSV.write_text("\n".join(rows) + "\n", encoding="utf-8")


def main() -> None:
    """生成 GL99 数字化曲线与 ASGARD 当前输出的叠图。"""
    gl_x = GL99_THETA_RATIO_09[:, 0]
    gl_p = GL99_THETA_RATIO_09[:, 1]
    geom_x, geom_p = _gl99_geometry_curve()
    asgard_x, asgard_p = _asgard_polarization_curve()
    order = np.argsort(asgard_x)
    asgard_x = asgard_x[order]
    asgard_p = asgard_p[order]
    _write_overlay_csv(gl_x, gl_p, geom_x, geom_p, asgard_x, asgard_p)

    fig, ax = plt.subplots(figsize=(6.8, 4.4), constrained_layout=True)
    ax.plot(gl_x, gl_p, "o-", color="black", lw=1.7, ms=4.0, label="GL99 Fig. 4 digitized, theta_o/theta_c=0.9")
    ax.plot(geom_x, geom_p, "-", color="#2ca02c", lw=1.8, label="GL99 geometry reproduced")
    ax.plot(asgard_x, asgard_p, "s-", color="#1f77b4", lw=1.7, ms=3.8, label="ASGARD shock_random, scaled to P0=60%")
    ax.set_xscale("log")
    ax.set_xlim(0.01, 0.55)
    ax.set_ylim(0.0, max(12.0, 1.1 * float(np.max(asgard_p))))
    ax.set_xlabel(r"$\Gamma^{-1}$")
    ax.set_ylabel("Linear polarization (%)")
    ax.set_title("Off-axis top-hat synchrotron polarization")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="upper right", fontsize=8)
    fig.savefig(OUTPUT_PNG, dpi=180)
    print(f"wrote {OUTPUT_PNG}")
    print(f"wrote {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
