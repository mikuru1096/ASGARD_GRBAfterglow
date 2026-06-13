from __future__ import annotations

from pathlib import Path
import argparse
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet
from asgard_core.asgard_paths import ASGARD_DOC_DIR


OUTPUT_DIR = ASGARD_DOC_DIR / "reverse_density_jump_tests"
OUTPUT_PATH = OUTPUT_DIR / "triple_density_jump_rs_fs_tophat.png"
JUMP_RADII_CM = (1.0e14, 1.0e15, 1.0e16)
JUMP_FACTOR = 1.0e3
JUMP_WIDTH_LOG10 = 1.0e-2
BANDS_HZ = np.array([1.0e9, 1.0e14, 1.0e18], dtype=float)
BAND_LABELS = ("1 GHz radio", "1e14 Hz optical", "1e18 Hz X-ray")
MODE_GRIDS = {
    "quick": {"num_r": 60, "num_theta": 32, "num_tobs": 32, "num_gam_e": 41, "num_nu": 31, "times": 80},
    "formal": {"num_r": 120, "num_theta": 80, "num_tobs": 48, "num_gam_e": 81, "num_nu": 49, "times": 160},
}


def _model(*, with_jumps: bool, grid: dict[str, int]) -> Model:
    jump_kwargs = {}
    if with_jumps:
        jump_kwargs = {
            "jump_r_cm": JUMP_RADII_CM,
            "jump_factor": (JUMP_FACTOR,) * len(JUMP_RADII_CM),
            "jump_width_log10": (JUMP_WIDTH_LOG10,) * len(JUMP_RADII_CM),
        }
    return Model(
        jet=TophatJet(theta_c=0.1, E_iso=1.0e52, Gamma0=300.0, duration=10.0),
        medium=ISM(n_ism=1.0),
        observer=Observer(lumi_dist=1.0e26, z=0.1, theta_obs=0.0),
        fwd_rad=Radiation(eps_e=3.0e-1, eps_B=1.0e-5, p=2.3, xi_N=1.0e-1, ssc=False),
        rvs_rad=Radiation(eps_e=3.0e-1, eps_B=1.0e-5, p=2.4, xi_N=1.0e-1, ssc=False),
        setups=Setups(
            rvs_shock=True,
            reverse_delta_t_s=10.0,
            fwd_ssc=False,
            rvs_ssc=False,
            ssc_cooling=False,
            num_threads=1,
            num_gam_e=grid["num_gam_e"],
            num_nu=grid["num_nu"],
            num_r=grid["num_r"],
            num_theta=grid["num_theta"],
            num_tobs=grid["num_tobs"],
            observer_time_min_s=1.0e0,
            observer_time_max_s=1.0e8,
            electron_solver="fullhide_1d",
            **jump_kwargs,
        ),
    )


def _density_enhancement(radius_cm: np.ndarray) -> np.ndarray:
    log_radius = np.log10(radius_cm)
    enhancement = np.ones_like(radius_cm, dtype=float)
    for radius_j in JUMP_RADII_CM:
        enhancement = enhancement + (JUMP_FACTOR - 1.0) * np.exp(
            -(log_radius - np.log10(radius_j)) ** 2 / (2.0 * JUMP_WIDTH_LOG10**2)
        )
    return enhancement


def _positive_for_log(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float).copy()
    arr[(arr <= 0.0) | ~np.isfinite(arr)] = np.nan
    return arr


def _ylim(curves: list[np.ndarray]) -> tuple[float, float]:
    positive = [curve[np.isfinite(curve) & (curve > 0.0)] for curve in curves]
    positive = [curve for curve in positive if curve.size > 0]
    if not positive:
        raise RuntimeError("light-curve panel has no positive finite flux values")
    values = np.concatenate(positive)
    ymax = float(np.nanmax(values))
    shown = values[values >= ymax * 1.0e-8]
    ymin = float(np.nanmin(shown if shown.size else values))
    return 10.0 ** (np.log10(ymin) - 0.25), 10.0 ** (np.log10(ymax) + 0.25)


def _plot_light_curve_panel(
    ax,
    times_s: np.ndarray,
    *,
    band_index: int,
    no_jump_total: np.ndarray,
    jump_total: np.ndarray,
    jump_fwd: np.ndarray,
    jump_rev: np.ndarray,
) -> None:
    curves = [no_jump_total, jump_total, jump_fwd, jump_rev]
    ax.loglog(times_s, _positive_for_log(no_jump_total), color="#666666", lw=1.5, ls="--", label="no-jump total")
    ax.loglog(times_s, _positive_for_log(jump_total), color="#0072B2", lw=2.0, label="triple-jump total")
    ax.loglog(times_s, _positive_for_log(jump_fwd), color="#009E73", lw=1.7, ls="-.", label="triple-jump FS")
    ax.loglog(times_s, _positive_for_log(jump_rev), color="#D55E00", lw=1.7, ls=":", label="triple-jump RS")
    ax.set_ylim(*_ylim(curves))
    ax.set_title(BAND_LABELS[band_index])
    ax.set_xlabel("Observer time [s]")
    ax.grid(True, which="both", alpha=0.25, linestyle=":")


def build_plot(*, mode: str, output: Path, times_count: int | None = None) -> Path:
    grid = dict(MODE_GRIDS[mode])
    if times_count is not None:
        grid["times"] = int(times_count)
    times = np.logspace(0.0, 8.0, grid["times"])
    no_jump = _model(with_jumps=False, grid=grid).flux_density_grid(times, BANDS_HZ)
    triple_jump = _model(with_jumps=True, grid=grid).flux_density_grid(times, BANDS_HZ)

    no_jump_total = np.asarray(no_jump.total, dtype=float)
    jump_total = np.asarray(triple_jump.total, dtype=float)
    jump_fwd = np.asarray(triple_jump.fwd.sync, dtype=float)
    jump_rev = np.asarray(triple_jump.rev.sync, dtype=float)

    fig = plt.figure(figsize=(12.0, 8.0), dpi=220)
    grid_spec = fig.add_gridspec(2, 3, height_ratios=(0.8, 2.2), hspace=0.35, wspace=0.28)
    density_ax = fig.add_subplot(grid_spec[0, :])
    radius = np.logspace(13.5, 16.5, 800)
    density_ax.semilogx(radius, _density_enhancement(radius), color="#000000", lw=1.8)
    for radius_j in JUMP_RADII_CM:
        density_ax.axvline(radius_j, color="#D55E00", lw=1.0, ls=":")
    density_ax.set_ylabel(r"$n/n_0$")
    density_ax.set_xlabel("Radius [cm]")
    density_ax.set_title("Triple ISM density jumps: factor 1000, width 0.01 dex")
    density_ax.grid(True, which="both", alpha=0.25, linestyle=":")

    axes = [fig.add_subplot(grid_spec[1, i]) for i in range(BANDS_HZ.size)]
    for i_band, ax in enumerate(axes):
        _plot_light_curve_panel(
            ax,
            times,
            band_index=i_band,
            no_jump_total=no_jump_total[i_band],
            jump_total=jump_total[i_band],
            jump_fwd=jump_fwd[i_band],
            jump_rev=jump_rev[i_band],
        )
    axes[0].set_ylabel(r"Flux density [erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$]")
    axes[-1].legend(fontsize=7.0, loc="best", frameon=False)
    fig.suptitle("ASGARD RS+FS top-hat light curves with three density jumps", y=0.98)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)

    for i_band, nu_hz in enumerate(BANDS_HZ):
        peak_index = int(np.nanargmax(jump_total[i_band]))
        print(
            f"band={nu_hz:.3e}Hz peak_time={times[peak_index]:.3e}s "
            f"peak_total={jump_total[i_band, peak_index]:.3e}"
        )
    print(f"wrote {output}")
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate the ASGARD triple-density-jump RS+FS top-hat benchmark.")
    parser.add_argument("--mode", choices=sorted(MODE_GRIDS), default="formal")
    parser.add_argument("--times", type=int, default=None)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    build_plot(mode=args.mode, output=args.output, times_count=args.times)


if __name__ == "__main__":
    main()
