from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from _repo_path import ensure_repo_root_on_path


ensure_repo_root_on_path()

from asgard_core import Observer, top_hat_jet
from tests.public_api_builders import numerics, observer_grid, radiation, solver_options, top_hat_model


SOLVERS = ("charint_1d", "charint_2d")
FREQUENCIES_HZ = np.array([1.0e9, 1.0e14, 1.0e18], dtype=float)
TIMES_S = np.logspace(2.0, 6.0, 32)


def _build_model(solver: str):
    return top_hat_model(
        jet=top_hat_jet(
            energy_iso_erg=1.0e52,
            initial_lorentz_factor=180.0,
            opening_angle_rad=0.12,
            shell_duration_s=None,
            magnetar=None,
            spreading=False,
        ),
        observer=Observer(z=0.08, viewing_angle_rad=0.04, viewing_azimuth_rad=0.0, luminosity_distance_cm=1.0e26),
        fwd_rad=radiation(epsilon_B=3.0e-3, accelerated_electron_fraction=0.2, include_ssc=False),
        numerics=numerics(
            num_radius=52,
            num_theta=48,
            num_phi=12,
            num_observer_time=52,
            num_electron_gamma=33,
            num_photon_frequency=45,
            num_chi=6 if solver.endswith("_2d") else None,
            electron_adaptive_substeps=True,
            electron_substep_rtol=0.04,
            electron_substep_min=8,
            electron_substep_max=128,
        ),
        observer_grid=observer_grid(time_min_s=float(TIMES_S[0]), time_max_s=float(TIMES_S[-1])),
        solver_options=solver_options(electron_solver=solver, ssc_cooling_mode="none"),
    )


def _compute_solver(solver: str) -> dict[str, np.ndarray]:
    model = _build_model(solver)
    flux = model.flux_density_grid(TIMES_S, FREQUENCIES_HZ).total
    details = model.details()
    gamma = np.asarray(details.fwd.gamma_e, dtype=float)
    dnde = np.asarray(details.fwd.dN_dgamma_e[:, -1], dtype=float)
    return {
        f"{solver}_flux": flux,
        f"{solver}_gamma": gamma,
        f"{solver}_dnde_final": dnde,
        f"{solver}_nu_m": np.asarray(details.fwd.nu_m, dtype=float),
        f"{solver}_nu_c": np.asarray(details.fwd.nu_c, dtype=float),
        f"{solver}_nu_a": np.asarray(details.fwd.nu_a, dtype=float),
        f"{solver}_electron_number": np.array([np.trapezoid(dnde, gamma)], dtype=float),
        f"{solver}_electron_min": np.array([np.min(dnde)], dtype=float),
    }


def _compute() -> dict[str, np.ndarray]:
    out: dict[str, np.ndarray] = {
        "time_s": TIMES_S,
        "frequency_hz": FREQUENCIES_HZ,
    }
    for solver in SOLVERS:
        out.update(_compute_solver(solver))
    return out


def _load_npz(path: Path) -> dict[str, np.ndarray]:
    with np.load(path) as data:
        return {key: data[key] for key in data.files}


def _plot_relative_change(ax, time: np.ndarray, before: np.ndarray, after: np.ndarray, label: str, color: str) -> None:
    mask = before > 0.0
    ax.semilogx(time[mask], after[mask] / before[mask] - 1.0, color=color, lw=1.3, label=label)


def _write_plot(before: dict[str, np.ndarray], after: dict[str, np.ndarray], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    time = after["time_s"]
    fig = plt.figure(figsize=(10.0, 7.6))
    gs = fig.add_gridspec(2, 2, height_ratios=(1.2, 1.0))
    ax_flux = fig.add_subplot(gs[0, :])
    ax_rel = fig.add_subplot(gs[1, 0])
    ax_spec = fig.add_subplot(gs[1, 1])
    colors = {"charint_1d": "#0173b2", "charint_2d": "#de8f05"}
    band_index = 1

    for solver in SOLVERS:
        color = colors[solver]
        ax_flux.loglog(
            time,
            before[f"{solver}_flux"][band_index],
            color=color,
            ls=":",
            lw=1.4,
            label=f"{solver} before",
        )
        ax_flux.loglog(
            time,
            after[f"{solver}_flux"][band_index],
            color=color,
            lw=1.6,
            label=f"{solver} after",
        )
        _plot_relative_change(
            ax_rel,
            time,
            before[f"{solver}_flux"][band_index],
            after[f"{solver}_flux"][band_index],
            solver,
            color,
        )
        ax_spec.loglog(
            before[f"{solver}_gamma"],
            before[f"{solver}_dnde_final"],
            color=color,
            ls=":",
            lw=1.4,
        )
        ax_spec.loglog(
            after[f"{solver}_gamma"],
            after[f"{solver}_dnde_final"],
            color=color,
            lw=1.6,
            label=solver,
        )

    ax_flux.set_ylabel(r"$F_\nu(10^{14}\,{\rm Hz})$")
    ax_flux.legend(frameon=False, fontsize=8, ncol=4)
    ax_flux.grid(True, which="both", alpha=0.25)
    ax_rel.axhline(0.0, color="black", lw=0.8)
    ax_rel.set_xlabel(r"$t_{\rm obs}$ (s)")
    ax_rel.set_ylabel("after / before - 1")
    ax_rel.legend(frameon=False, fontsize=8)
    ax_rel.grid(True, which="both", alpha=0.25)
    ax_spec.set_xlabel(r"$\gamma_e$")
    ax_spec.set_ylabel(r"$dN_e/d\gamma_e$ at final shell")
    ax_spec.legend(frameon=False, fontsize=8)
    ax_spec.grid(True, which="both", alpha=0.25)
    fig.suptitle("Characteristic electron transport: before/after comparison")
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--save-baseline", type=Path)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    current = _compute()
    if args.save_baseline is not None:
        args.save_baseline.parent.mkdir(parents=True, exist_ok=True)
        np.savez(args.save_baseline, **current)
    if args.baseline is not None:
        if args.output is None:
            raise SystemExit("--output is required with --baseline")
        _write_plot(_load_npz(args.baseline), current, args.output)


if __name__ == "__main__":
    main()
