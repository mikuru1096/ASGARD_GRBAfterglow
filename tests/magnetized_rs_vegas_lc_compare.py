from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from _repo_path import ensure_repo_root_on_path


ensure_repo_root_on_path()

from asgard_core import Model, Observer, UniformMedium, WindMedium, top_hat_jet
from tests.public_api_builders import hadronic, numerics, observer_grid, radiation, reverse_shock, solver_options

from VegasAfterglow import ISM as VegasISM
from VegasAfterglow import Model as VegasModel
from VegasAfterglow import Observer as VegasObserver
from VegasAfterglow import Radiation as VegasRadiation
from VegasAfterglow import Wind as VegasWind
from VegasAfterglow.VegasAfterglowC import Ejecta as VegasEjecta


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIRS = {
    "ism": ROOT / "output" / "asgard_doc" / "magnetized_rs_sigma_benchmark",
    "wind": ROOT / "output" / "asgard_doc" / "magnetized_rs_sigma_benchmark_wind",
}
SIGMAS = (0.0, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0)
BANDS_HZ = np.array([1.0e9, 1.0e14], dtype=float)
BAND_LABELS = ("1 GHz", "1e14 Hz")
E_ISO = 1.0e54
GAMMA0 = 50.0
THETA_J = 1.0e-1
DURATION_S = 500.0
A_STAR = 5.0e-2
WIND_N_ISM = 1.0e-15
ISM_N = 1.0e-1
Z = 0.4
LUMI_DIST_CM = 1.0e26
FWD_EPS_E = 1.0e-1
FWD_EPS_B = 1.0e-3
FWD_P = 2.5
FWD_XI = 1.0e-1
RVS_EPS_E = 0.3
RVS_EPS_B = 0.1
RVS_P = 2.4
RVS_XI = 1.0


def _top_hat_field(inside: float, outside: float = 0.0):
    def field(_phi: float, theta: float) -> float:
        return inside if theta <= THETA_J else outside

    return field


def _output_path(medium: str) -> Path:
    suffix = "wind" if medium == "wind" else "ism"
    return OUTPUT_DIRS[medium] / f"magnetized_rs_sigma_{suffix}_vegas_lc_compare.png"


def _asgard_medium(medium: str):
    if medium == "ism":
        return UniformMedium(number_density_cm3=ISM_N)
    if medium == "wind":
        return WindMedium(a_star=A_STAR, density_floor_cm3=WIND_N_ISM, density_cap_cm3=None)
    raise ValueError(f"unknown medium {medium!r}")


def _vegas_medium(medium: str):
    if medium == "ism":
        return VegasISM(n_ism=ISM_N)
    if medium == "wind":
        return VegasWind(A_star=A_STAR, n_ism=WIND_N_ISM)
    raise ValueError(f"unknown medium {medium!r}")


def _asgard_model(
    sigma: float,
    medium: str,
    *,
    cooling_mode: str,
    num_r: int,
    num_theta: int,
    num_tobs: int,
) -> Model:
    return Model(
        jet=top_hat_jet(
            energy_iso_erg=E_ISO,
            initial_lorentz_factor=GAMMA0,
            opening_angle_rad=THETA_J,
            shell_duration_s=DURATION_S,
            magnetar=None,
            spreading=False,
        ),
        medium=_asgard_medium(medium),
        observer=Observer(
            z=Z,
            viewing_angle_rad=0.0,
            viewing_azimuth_rad=0.0,
            luminosity_distance_cm=LUMI_DIST_CM,
        ),
        fwd_rad=radiation(
            epsilon_e=FWD_EPS_E,
            epsilon_B=FWD_EPS_B,
            p=FWD_P,
            accelerated_electron_fraction=FWD_XI,
            include_ssc=False,
            include_kn_correction=False,
        ),
        rvs_rad=radiation(
            epsilon_e=RVS_EPS_E,
            epsilon_B=RVS_EPS_B,
            p=RVS_P,
            accelerated_electron_fraction=RVS_XI,
            include_ssc=False,
            include_kn_correction=False,
        ),
        numerics=numerics(
            num_radius=num_r,
            num_theta=num_theta,
            num_phi=1,
            num_observer_time=num_tobs,
            num_electron_gamma=41,
            num_photon_frequency=51,
            num_threads=1,
            initial_radius_cm=1.0e14,
        ),
        observer_grid=observer_grid(time_min_s=1.0, time_max_s=1.0e7),
        solver_options=solver_options(
            electron_solver="fullhide_1d",
            ssc_cooling_mode=cooling_mode,
        ),
        reverse_shock=reverse_shock(
            enabled=True,
            shell_duration_s=DURATION_S,
            upstream_sigma=sigma,
            include_ssc=False,
            include_cross_zone_ic=False,
        ),
        hadronic=hadronic(),
    )


def _vegas_model(sigma: float, medium: str) -> VegasModel:
    return VegasModel(
        jet=VegasEjecta(
            E_iso=_top_hat_field(E_ISO),
            Gamma0=_top_hat_field(GAMMA0, GAMMA0),
            sigma0=_top_hat_field(sigma),
            duration=DURATION_S,
        ),
        medium=_vegas_medium(medium),
        observer=VegasObserver(lumi_dist=LUMI_DIST_CM, z=Z, theta_obs=0.0),
        fwd_rad=VegasRadiation(eps_e=FWD_EPS_E, eps_B=FWD_EPS_B, p=FWD_P, xi_e=FWD_XI, ssc=False),
        rvs_rad=VegasRadiation(eps_e=RVS_EPS_E, eps_B=RVS_EPS_B, p=RVS_P, xi_e=RVS_XI, ssc=False),
        resolutions=(0.06, 0.14, 10),
    )


def _positive(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    return arr[np.isfinite(arr) & (arr > 0.0)]


def _ylim(curves: list[np.ndarray], rel_floor: float = 1.0e-8) -> tuple[float, float]:
    positives = [_positive(curve) for curve in curves]
    panel = np.concatenate([curve for curve in positives if curve.size])
    ymax = float(np.nanmax(panel))
    shown = panel[panel >= ymax * rel_floor]
    ymin = float(np.nanmin(shown if shown.size else panel))
    return 10.0 ** (np.log10(ymin) - 0.25), 10.0 ** (np.log10(ymax) + 0.25)


def _plot_curve(
    ax,
    times: np.ndarray,
    values: np.ndarray,
    *,
    color: str,
    linestyle,
    label: str | None,
    linewidth: float,
) -> None:
    arr = np.asarray(values, dtype=float)
    arr = np.where((arr > 0.0) & np.isfinite(arr), arr, np.nan)
    ax.loglog(times, arr, color=color, linestyle=linestyle, linewidth=linewidth, label=label)


def build_plot(
    times: np.ndarray,
    *,
    medium: str,
    cooling_mode: str,
    num_r: int,
    num_theta: int,
    num_tobs: int,
    output: Path,
) -> None:
    data: dict[tuple[str, float], tuple[np.ndarray, np.ndarray]] = {}
    for sigma in SIGMAS:
        as_flux = _asgard_model(
            sigma,
            medium,
            cooling_mode=cooling_mode,
            num_r=num_r,
            num_theta=num_theta,
            num_tobs=num_tobs,
        ).flux_density_grid(times, BANDS_HZ)
        vg_flux = _vegas_model(sigma, medium).flux_density_grid(times, BANDS_HZ)
        data[("ASGARD", sigma)] = (np.asarray(as_flux.total, dtype=float), np.asarray(as_flux.rev.sync, dtype=float))
        vg_total = np.asarray(vg_flux.total, dtype=float)
        vg_rs = np.asarray(vg_flux.rvs.sync, dtype=float)
        if medium == "ism" and sigma >= 1.0e-5 and not np.any(np.isfinite(vg_rs) & (vg_rs > 0.0)):
            raise RuntimeError(
                "VegasAfterglow returned zero ISM magnetized RS synchrotron; "
                "install the patched 2.0.5 source recorded in doc/patches/"
                "vegasafterglow-2.0.5-rvs-ism-compression.patch."
            )
        data[("Vegas", sigma)] = (vg_total, vg_rs)

    fig, axes = plt.subplots(len(BANDS_HZ), len(SIGMAS), figsize=(24.0, 7.4), dpi=220, sharex=True, squeeze=False)
    colors = {"ASGARD": "#0072B2", "Vegas": "#D55E00"}
    for i_band, band_label in enumerate(BAND_LABELS):
        for i_sigma, sigma in enumerate(SIGMAS):
            ax = axes[i_band, i_sigma]
            as_total, as_rs = data[("ASGARD", sigma)]
            vg_total, vg_rs = data[("Vegas", sigma)]
            _plot_curve(
                ax,
                times,
                as_total[i_band],
                color=colors["ASGARD"],
                linestyle="-",
                linewidth=2.0,
                label="ASGARD total" if i_band == 0 and i_sigma == 0 else None,
            )
            _plot_curve(
                ax,
                times,
                as_rs[i_band],
                color=colors["ASGARD"],
                linestyle=":",
                linewidth=2.2,
                label="ASGARD RS" if i_band == 0 and i_sigma == 0 else None,
            )
            _plot_curve(
                ax,
                times,
                vg_total[i_band],
                color=colors["Vegas"],
                linestyle=(0, (5, 2)),
                linewidth=2.0,
                label="Vegas total" if i_band == 0 and i_sigma == 0 else None,
            )
            _plot_curve(
                ax,
                times,
                vg_rs[i_band],
                color=colors["Vegas"],
                linestyle=(0, (1.2, 1.2)),
                linewidth=2.2,
                label="Vegas RS" if i_band == 0 and i_sigma == 0 else None,
            )
            ax.set_ylim(*_ylim([as_total[i_band], as_rs[i_band], vg_total[i_band], vg_rs[i_band]]))
            ax.grid(True, which="both", alpha=0.25, linestyle=":")
            if i_band == 0:
                ax.set_title(f"sigma={sigma:g}")
            if i_sigma == 0:
                ax.set_ylabel(f"Fnu {band_label}")
    for ax in axes[-1, :]:
        ax.set_xlabel("observer time [s]")
    axes[0, 0].set_xlim(times[0], times[-1])
    fig.legend(loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.99))
    fig.suptitle(f"{medium.upper()} reverse-shock light curves: ASGARD vs VegasAfterglow", y=0.935)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.92))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    plt.close(fig)

    for sigma in SIGMAS:
        for backend in ("ASGARD", "Vegas"):
            total, rs = data[(backend, sigma)]
            msg = []
            for i_band, band_label in enumerate(BAND_LABELS):
                idx = int(np.nanargmax(rs[i_band]))
                ratio = float(rs[i_band, idx] / total[i_band, idx]) if total[i_band, idx] > 0.0 else float("nan")
                msg.append(f"{band_label} RS_peak={rs[i_band, idx]:.3e} ratio_at_rs_peak={ratio:.3g}")
            print(f"{medium} sigma={sigma:g} {backend}: " + "; ".join(msg))
    print(f"wrote {output}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--medium", choices=("ism", "wind", "both"), default="both")
    parser.add_argument("--times", type=int, default=240)
    parser.add_argument("--num-r", type=int, default=80)
    parser.add_argument("--num-tobs", type=int, default=48)
    parser.add_argument("--num-theta", type=int, default=20)
    parser.add_argument(
        "--cooling-mode",
        choices=("none", "numeric_ic_kn", "nakar_y_thomson"),
        default="nakar_y_thomson",
    )
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    if args.output is not None and args.medium == "both":
        raise ValueError("--output is only valid when --medium is ism or wind")
    media = ("ism", "wind") if args.medium == "both" else (args.medium,)
    times = np.logspace(0.0, 7.0, args.times)
    for medium in media:
        output = _output_path(medium) if args.output is None else args.output
        build_plot(
            times,
            medium=medium,
            cooling_mode=args.cooling_mode,
            num_r=args.num_r,
            num_theta=args.num_theta,
            num_tobs=args.num_tobs,
            output=output,
        )


if __name__ == "__main__":
    main()
