from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np

from _repo_path import ensure_repo_root_on_path


ensure_repo_root_on_path()

from tests.public_api_builders import numerics, radiation, solver_options, top_hat_model


ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper"
FIG_DIR = PAPER / "figures"
DATA_DIR = PAPER / "source_data"

BLUE = "#24588D"
TEAL = "#2E8B8B"
RED = "#B34A48"
GREEN = "#3E8F4D"
VIOLET = "#7750A6"
GOLD = "#B8860B"
NEUTRAL = "#4A4A4A"
LIGHT = "#E6E6E6"

C_CGS = 2.99792458e10
M_E_G = 9.1093837015e-28
M_P_G = 1.67262192369e-24
M_PI0_G = 2.406176653e-25
MEC2_ERG = M_E_G * C_CGS**2
EV_TO_ERG = 1.602176634e-12
M_E_C2_EV = MEC2_ERG / EV_TO_ERG
M_P_C2_EV = M_P_G * C_CGS**2 / EV_TO_ERG
M_PI0_C2_EV = M_PI0_G * C_CGS**2 / EV_TO_ERG
EPS_PG_THR_EV = M_PI0_C2_EV * (1.0 + M_PI0_C2_EV / (2.0 * M_P_C2_EV))
EPS_BH_THR_EV = 2.0 * M_E_C2_EV
TARGET_PHOTON_MIN_EV = 1.0e-6
TARGET_PHOTON_MAX_EV = 1.0e9
OPTICAL_DEPTH_MIN = 1.0e-3
OPTICAL_DEPTH_MAX = 1.0e2
ATTENUATION_DISPLAY_MIN = np.exp(-OPTICAL_DEPTH_MAX**0.5)
HEATMAP_FLOOR_FRACTION = 1.0e-18
DISPLAY_FLOOR_FRACTION = 1.0e-12


plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans", "sans-serif"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.size": 7,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.8,
        "legend.frameon": False,
    }
)


def ensure_dirs() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)


def save_pub(fig: plt.Figure, stem: str) -> None:
    for suffix, kwargs in (("svg", {}), ("pdf", {}), ("tiff", {"dpi": 600})):
        path = FIG_DIR / f"{stem}.{suffix}"
        fig.savefig(path, bbox_inches="tight", **kwargs)
        if suffix == "svg":
            strip_trailing_whitespace(path)
    plt.close(fig)


def strip_trailing_whitespace(path: Path) -> None:
    text = path.read_text(encoding="utf-8")
    cleaned = "\n".join(line.rstrip() for line in text.splitlines()) + "\n"
    path.write_text(cleaned, encoding="utf-8", newline="\n")


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"no rows for {path}")
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


def add_panel(ax: plt.Axes, label: str) -> None:
    ax.text(-0.10, 1.04, label, transform=ax.transAxes, fontweight="bold", fontsize=9)


def draw_flow_box(ax: plt.Axes, xy: tuple[float, float], text: str, color: str) -> None:
    x, y = xy
    patch = FancyBboxPatch(
        (x, y),
        0.31,
        0.11,
        boxstyle="round,pad=0.015,rounding_size=0.02",
        facecolor=color,
        edgecolor="none",
        alpha=0.14,
    )
    ax.add_patch(patch)
    ax.text(x + 0.155, y + 0.055, text, ha="center", va="center", fontsize=6.4, color=NEUTRAL)


def draw_flow_arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float]) -> None:
    arrow = FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=8, lw=0.8, color=NEUTRAL)
    ax.add_patch(arrow)


def diagnostic_model():
    return top_hat_model(
        fwd_rad=radiation(include_ssc=True, include_kn_correction=False),
        numerics=numerics(
            num_radius=24,
            eats_num_theta=12,
            num_observer_time=24,
            num_electron_gamma=32,
            num_photon_frequency=32,
            num_threads=1,
        ),
    )


def fig1_radius_state() -> None:
    model = diagnostic_model()
    details = model.details(1.0e4, 2.0e4).fwd
    radius = np.asarray(details.radius, dtype=float)
    gamma = np.asarray(details.Gamma, dtype=float)
    b_comv = np.asarray(details.B_comv, dtype=float)
    beta = np.sqrt(1.0 - gamma ** -2)
    theta_j = 0.10
    selected = np.linspace(0, radius.size - 1, min(6, radius.size), dtype=int)
    mu_axis = 1.0
    mu_edge = np.cos(theta_j)
    t_axis = radius * (1.0 / beta - mu_axis) / C_CGS
    dt_on = (1.0 / beta - 1.0) / C_CGS
    dt_edge = (1.0 / beta - mu_edge) / C_CGS

    rows: list[dict[str, object]] = []
    flow_nodes = [
        ("public inputs", "public inputs"),
        ("runtime", "runtime configuration"),
        ("dynamics", "dynamics"),
        ("electron_photon", "electron and photon state"),
        ("optional_branches", "optional branches"),
        ("observer", "observer projection"),
    ]
    for key, label in flow_nodes:
        rows.append(
            {
                "kind": "contract_node",
                "radius_cm": "",
                "Gamma": "",
                "B_comoving_G": "",
                "tobs_axis_s": "",
                "dtobs_dR_on_axis_s_per_cm": "",
                "dtobs_dR_edge_s_per_cm": "",
                "node": key,
                "label": label,
            }
        )
    for r, g, b, ta, do, de in zip(radius, gamma, b_comv, t_axis, dt_on, dt_edge):
        rows.append(
            {
                "kind": "asgard_state",
                "radius_cm": f"{r:.8e}",
                "Gamma": f"{g:.8e}",
                "B_comoving_G": f"{b:.8e}",
                "tobs_axis_s": f"{ta:.8e}",
                "dtobs_dR_on_axis_s_per_cm": f"{do:.8e}",
                "dtobs_dR_edge_s_per_cm": f"{de:.8e}",
                "node": "",
                "label": "",
            }
        )
    write_rows(DATA_DIR / "fig1_radius_state.csv", rows)

    fig = plt.figure(figsize=(7.1, 3.25))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.08, 1.0], height_ratios=[1.0, 1.0])
    fig.subplots_adjust(left=0.065, right=0.985, bottom=0.15, top=0.91, wspace=0.34, hspace=0.78)
    ax_flow = fig.add_subplot(gs[:, 0])
    ax_state = fig.add_subplot(gs[0, 1])
    ax_time = fig.add_subplot(gs[1, 1])

    ax_flow.set_axis_off()
    add_panel(ax_flow, "a")
    ax_flow.set_title("Radius-ordered calculation contract", pad=8)
    positions = [(0.04, 0.72), (0.40, 0.72), (0.04, 0.48), (0.40, 0.48), (0.04, 0.24), (0.40, 0.24)]
    labels = [
        "public\ninputs",
        "runtime\nconfiguration",
        "dynamics\n$\\Gamma(R),B'(R)$",
        "electron + photon\n$N_e,n_\\gamma$",
        "bounded optional\nRS/hadronic",
        "observer\n$F_\\nu(t_{\\rm obs})$",
    ]
    colors = [BLUE, BLUE, TEAL, TEAL, GOLD, RED]
    for pos, label, color in zip(positions, labels, colors):
        draw_flow_box(ax_flow, pos, label, color)
    arrow_pairs = [
        ((0.35, 0.775), (0.40, 0.775)),
        ((0.555, 0.72), (0.195, 0.59)),
        ((0.35, 0.535), (0.40, 0.535)),
        ((0.555, 0.48), (0.195, 0.35)),
        ((0.35, 0.295), (0.40, 0.295)),
    ]
    for start, end in arrow_pairs:
        draw_flow_arrow(ax_flow, start, end)
    ax_flow.text(
        0.50,
        0.08,
        r"local source terms are stored before EATS projection",
        ha="center",
        va="center",
        fontsize=6.3,
        color=NEUTRAL,
    )
    ax_flow.set_xlim(0.0, 0.78)
    ax_flow.set_ylim(0.0, 0.90)

    ax_state.loglog(radius, gamma, color=BLUE, lw=1.7)
    ax0b = ax_state.twinx()
    ax0b.loglog(radius, b_comv, color=TEAL, lw=1.5)
    add_panel(ax_state, "b")
    ax_state.set_xlabel(r"$R$ (cm)")
    ax_state.set_ylabel(r"$\Gamma$", color=BLUE)
    ax0b.set_ylabel(r"$B'$ (G)", color=TEAL, labelpad=8)
    ax_state.set_title("Local blast-wave state")
    ax_state.grid(color=LIGHT, lw=0.5, which="both")

    ax_time.loglog(radius, t_axis, color=NEUTRAL, lw=1.5, label=r"$\mu=1$")
    ax_time.scatter(radius[selected], t_axis[selected], color=RED, s=14, zorder=3)
    ax_time_t = ax_time.twinx()
    ax_time_t.semilogx(radius, C_CGS * dt_edge, color=VIOLET, lw=1.2, ls="--", label=r"$\mu=\cos\theta_j$")
    add_panel(ax_time, "c")
    ax_time.set_xlabel(r"$R$ (cm)")
    ax_time.set_ylabel(r"$t_{\rm obs}$ (s)")
    ax_time_t.set_ylabel(r"$c\,\partial t_{\rm obs}/\partial R$", color=VIOLET)
    ax_time.set_title("EATS mapping")
    ax_time.grid(color=LIGHT, lw=0.5, which="both")
    save_pub(fig, "fig1_radius_state")


def fig2_forward_spectrum() -> None:
    model = diagnostic_model()
    details = model.details(1.0e4, 2.0e4).fwd
    radius = np.asarray(details.radius, dtype=float)
    gamma = np.asarray(details.Gamma, dtype=float)
    b_comv = np.asarray(details.B_comv, dtype=float)
    times = np.logspace(2.0, 7.0, 42)
    freqs = np.array([1.0e9, 1.0e14, 1.0e18])
    spec_freq = np.logspace(7.0, 21.0, 80)
    spec_times = np.array([1.0e3, 1.0e5, 1.0e7])
    flux = model.flux_density_grid(times, freqs)
    spec_flux = model.flux_density_grid(spec_times, spec_freq)

    rows: list[dict[str, object]] = []
    for r, g, b in zip(radius, gamma, b_comv):
        rows.append(
            {
                "kind": "dynamics_state",
                "radius_cm": f"{r:.8e}",
                "Gamma": f"{g:.8e}",
                "B_comoving_G": f"{b:.8e}",
                "time_s": "",
                "frequency_hz": "",
                "total_fnu_cgs": "",
                "fs_synch_fnu_cgs": "",
                "fs_ssc_fnu_cgs": "",
            }
        )
    for i, nu in enumerate(freqs):
        for t, total, sync, ssc in zip(times, flux.total[i], flux.fwd.sync[i], flux.fwd.ssc[i]):
            rows.append(
                {
                    "kind": "lightcurve",
                    "radius_cm": "",
                    "Gamma": "",
                    "B_comoving_G": "",
                    "time_s": f"{t:.8e}",
                    "frequency_hz": f"{nu:.8e}",
                    "total_fnu_cgs": f"{total:.8e}",
                    "fs_synch_fnu_cgs": f"{sync:.8e}",
                    "fs_ssc_fnu_cgs": f"{ssc:.8e}",
                }
            )
    for j, t in enumerate(spec_times):
        for nu, total, sync, ssc in zip(
            spec_freq,
            spec_flux.total[:, j],
            spec_flux.fwd.sync[:, j],
            spec_flux.fwd.ssc[:, j],
        ):
            rows.append(
                {
                    "kind": f"spectrum_t{t:.0e}s",
                    "radius_cm": "",
                    "Gamma": "",
                    "B_comoving_G": "",
                    "time_s": f"{t:.8e}",
                    "frequency_hz": f"{nu:.8e}",
                    "total_fnu_cgs": f"{total:.8e}",
                    "fs_synch_fnu_cgs": f"{sync:.8e}",
                    "fs_ssc_fnu_cgs": f"{ssc:.8e}",
                }
            )
    write_rows(DATA_DIR / "fig2_forward_api.csv", rows)

    fig, axes = plt.subplots(2, 2, figsize=(7.1, 4.35))
    fig.subplots_adjust(left=0.08, right=0.985, bottom=0.12, top=0.92, wspace=0.45, hspace=0.52)
    colors = [BLUE, TEAL, RED]
    labels = ["1 GHz", "optical", "X-ray"]
    ax = axes[0, 0]
    ax.loglog(radius, gamma, color=BLUE, lw=1.6)
    axb = ax.twinx()
    axb.loglog(radius, b_comv, color=TEAL, lw=1.4)
    add_panel(ax, "a")
    ax.set_xlabel(r"$R$ (cm)")
    ax.set_ylabel(r"$\Gamma$", color=BLUE)
    axb.set_ylabel(r"$B'$ (G)", color=TEAL)
    ax.set_title("State entering observer products")
    ax.grid(color=LIGHT, lw=0.5, which="both")

    ax = axes[0, 1]
    for i, (color, label) in enumerate(zip(colors, labels)):
        ax.loglog(times, flux.total[i], color=color, lw=1.6, label=label)
    add_panel(ax, "b")
    ax.set_xlabel(r"$t_{\rm obs}$ (s)")
    ax.set_ylabel(r"$F_\nu$ (cgs)")
    ax.set_title("Multi-frequency light curves")
    ax.legend(fontsize=6)
    ax.grid(color=LIGHT, lw=0.5, which="both")

    ax = axes[1, 0]
    for j, (t, color) in enumerate(zip(spec_times, [BLUE, TEAL, RED])):
        ax.loglog(spec_freq, spec_flux.total[:, j], color=color, lw=1.5, label=rf"$10^{{{int(np.log10(t))}}}$ s")
    add_panel(ax, "c")
    ax.set_xlabel(r"$\nu$ (Hz)")
    ax.set_ylabel(r"$F_\nu$ (cgs)")
    ax.set_title("Spectral evolution")
    ax.legend(fontsize=6)
    ax.grid(color=LIGHT, lw=0.5, which="both")

    ssc_fraction = np.divide(
        spec_flux.fwd.ssc[:, 1],
        spec_flux.total[:, 1],
        out=np.zeros_like(spec_flux.fwd.ssc[:, 1]),
        where=spec_flux.total[:, 1] > 0.0,
    )
    ax = axes[1, 1]
    ax.semilogx(spec_freq, ssc_fraction, color=VIOLET, lw=1.7)
    add_panel(ax, "d")
    ax.set_xlabel(r"$\nu$ (Hz)")
    ax.set_ylabel(r"$F_\nu^{\rm SSC}/F_\nu$")
    ax.set_ylim(0.0, min(1.0, float(np.nanmax(ssc_fraction)) * 1.2 + 0.02))
    ax.set_title(r"Component role at $10^5$ s")
    ax.grid(color=LIGHT, lw=0.5, which="both")
    save_pub(fig, "fig2_forward_api")


def fig3_electron_transport() -> None:
    model = top_hat_model(
        fwd_rad=radiation(include_ssc=False),
        numerics=numerics(
            num_radius=12,
            eats_num_theta=8,
            num_observer_time=12,
            num_electron_gamma=16,
            num_photon_frequency=16,
            num_threads=1,
        ),
    )
    details = model.details(1.0e4, 2.0e4).fwd
    chi_model = top_hat_model(
        fwd_rad=radiation(include_ssc=False),
        numerics=numerics(
            num_radius=8,
            eats_num_theta=6,
            num_observer_time=8,
            num_electron_gamma=12,
            num_photon_frequency=12,
            downstream_num_chi=6,
            num_threads=1,
        ),
        solver_options=solver_options(
            electron_solver="fullhide_2d",
            geometry_projection="chi_eats_2d",
        ),
    )
    chi_details = chi_model.details(1.0e4, 2.0e4).fwd
    radius = np.asarray(details.radius, dtype=float)
    gamma_e = np.asarray(details.gamma_e, dtype=float)
    electron = np.asarray(details.dN_dgamma_e, dtype=float)
    chi_radius = np.asarray(chi_details.chi_radius_cm, dtype=float)
    chi_grid = np.asarray(chi_details.chi_grid, dtype=float)
    chi_weight = np.asarray(chi_details.chi_dvolume_weight, dtype=float)
    tau_chi = np.asarray(chi_details.tau_syn_chi, dtype=float)
    seed_chi = np.asarray(chi_details.seed_syn_chi, dtype=float)
    tau_max = np.nanmax(tau_chi, axis=0)
    seed_max = np.nanmax(seed_chi, axis=0)
    xgrid = np.log10(gamma_e)
    dx = np.gradient(xgrid)
    number_per_x = electron * gamma_e[:, None] * np.log(10.0)
    cell_number = number_per_x * dx[:, None]
    energy_per_shell = np.sum((gamma_e[:, None] - 1.0) * cell_number, axis=0) * MEC2_ERG
    total_number = np.sum(cell_number, axis=0)

    rows: list[dict[str, object]] = []
    for j, r in enumerate(radius):
        rows.append(
            {
                "kind": "shell_budget",
                "radius_cm": f"{r:.8e}",
                "gamma_e": "",
                "dN_dgamma_e": "",
                "electron_number": f"{total_number[j]:.8e}",
                "electron_energy_erg": f"{energy_per_shell[j]:.8e}",
                "chi": "",
                "chi_radius_cm": "",
                "chi_dvolume_weight": "",
                "max_tau_syn_chi": "",
                "max_seed_syn_chi": "",
            }
        )
    for i, g in enumerate(gamma_e):
        for j, r in enumerate(radius):
            rows.append(
                {
                    "kind": "electron_spectrum",
                    "radius_cm": f"{r:.8e}",
                    "gamma_e": f"{g:.8e}",
                    "dN_dgamma_e": f"{electron[i, j]:.8e}",
                    "electron_number": "",
                    "electron_energy_erg": "",
                    "chi": "",
                    "chi_radius_cm": "",
                    "chi_dvolume_weight": "",
                    "max_tau_syn_chi": "",
                    "max_seed_syn_chi": "",
                }
            )
    for k, chi in enumerate(chi_grid):
        for j in range(chi_radius.shape[1]):
            rows.append(
                {
                    "kind": "chi_volume_weight",
                    "radius_cm": "",
                    "gamma_e": "",
                    "dN_dgamma_e": "",
                    "electron_number": "",
                    "electron_energy_erg": "",
                    "chi": f"{chi:.8e}",
                    "chi_radius_cm": f"{chi_radius[k, j]:.8e}",
                    "chi_dvolume_weight": f"{chi_weight[k, j]:.8e}",
                    "max_tau_syn_chi": f"{tau_max[k, j]:.8e}",
                    "max_seed_syn_chi": f"{seed_max[k, j]:.8e}",
                }
            )
    write_rows(DATA_DIR / "fig3_electron_transport.csv", rows)

    fig, axes = plt.subplots(2, 2, figsize=(7.1, 4.35))
    fig.subplots_adjust(left=0.08, right=0.985, bottom=0.12, top=0.91, wspace=0.45, hspace=0.58)
    ax = axes[0, 0]
    pcm = ax.pcolormesh(
        np.log10(radius),
        np.log10(gamma_e),
        np.log10(np.maximum(electron, np.max(electron) * HEATMAP_FLOOR_FRACTION)),
        shading="auto",
        cmap="viridis",
    )
    add_panel(ax, "a")
    ax.set_xlabel(r"$\log_{10} R$ (cm)")
    ax.set_ylabel(r"$\log_{10}\gamma_e$")
    ax.set_title(r"$N_e(\gamma_e,R)$")
    fig.colorbar(pcm, ax=ax, fraction=0.045, pad=0.02, label=r"$\log N_e$")

    ax = axes[0, 1]
    for idx, color in zip(np.linspace(0, radius.size - 1, 4, dtype=int), [BLUE, TEAL, RED, VIOLET]):
        spec = gamma_e**2 * np.maximum(electron[:, idx], 0.0)
        norm = np.max(spec)
        spec = spec / norm if norm > 0.0 else spec
        ax.loglog(gamma_e, np.maximum(spec, DISPLAY_FLOOR_FRACTION), color=color, lw=1.4, label=rf"$R={radius[idx]:.1e}$")
    add_panel(ax, "b")
    ax.set_xlabel(r"$\gamma_e$")
    ax.set_ylabel(r"normalized $\gamma_e^2 dN_e/d\gamma_e$")
    ax.set_ylim(1.0e-6, 2.0)
    ax.set_title("Electron spectral aging")
    ax.legend(fontsize=5.4)
    ax.grid(color=LIGHT, lw=0.5, which="both")

    ax = axes[1, 0]
    chi_plot = np.log10(np.maximum(chi_weight, np.nanmax(chi_weight) * DISPLAY_FLOOR_FRACTION))
    pcm = ax.pcolormesh(
        np.log10(np.maximum(chi_radius, np.nanmin(chi_radius[chi_radius > 0]))),
        np.log10(chi_grid)[:, None] * np.ones_like(chi_radius),
        chi_plot,
        shading="auto",
        cmap="magma",
    )
    add_panel(ax, "c")
    ax.set_xlabel(r"$\log_{10} R_\chi$ (cm)")
    ax.set_ylabel(r"$\log_{10}\chi$")
    ax.set_title(r"Finite-thickness volume weight")
    fig.colorbar(pcm, ax=ax, fraction=0.045, pad=0.02, label=r"$\log w_\chi$")

    ax = axes[1, 1]
    tau_plot = np.log10(np.maximum(tau_max, np.nanmax(tau_max) * DISPLAY_FLOOR_FRACTION))
    pcm = ax.pcolormesh(
        np.log10(np.maximum(chi_radius, np.nanmin(chi_radius[chi_radius > 0]))),
        np.log10(chi_grid)[:, None] * np.ones_like(chi_radius),
        tau_plot,
        shading="auto",
        cmap="cividis",
    )
    add_panel(ax, "d")
    ax.set_xlabel(r"$\log_{10} R_\chi$ (cm)")
    ax.set_ylabel(r"$\log_{10}\chi$")
    ax.set_title(r"Stored synchrotron opacity")
    fig.colorbar(pcm, ax=ax, fraction=0.045, pad=0.02, label=r"$\log\max_\nu\tau_\nu$")
    save_pub(fig, "fig3_transport_projection")


def figA1_hadronic_thresholds() -> None:
    eps_ev = np.logspace(np.log10(TARGET_PHOTON_MIN_EV), np.log10(TARGET_PHOTON_MAX_EV), 240)
    eps_erg = eps_ev * EV_TO_ERG
    gamma_p_pg = EPS_PG_THR_EV / (2.0 * eps_ev)
    gamma_p_bh = EPS_BH_THR_EV / (2.0 * eps_ev)
    e_gamma_gg_ev = M_E_C2_EV**2 / eps_ev
    tau = np.logspace(np.log10(OPTICAL_DEPTH_MIN), np.log10(OPTICAL_DEPTH_MAX), 220)
    survival = np.exp(-tau)
    cell_transfer = (1.0 - np.exp(-tau)) / tau
    broken = np.where(eps_ev < 1.0, eps_ev ** 0.4, eps_ev ** -0.8)
    photon_density = broken / eps_erg
    photon_density = photon_density / np.nanmax(photon_density)

    rows: list[dict[str, object]] = []
    for e, gp, gbh, egg, nd in zip(eps_ev, gamma_p_pg, gamma_p_bh, e_gamma_gg_ev, photon_density):
        rows.append(
            {
                "kind": "threshold",
                "target_photon_ev": f"{e:.8e}",
                "pgamma_gamma_p_threshold": f"{gp:.8e}",
                "bh_gamma_p_threshold": f"{gbh:.8e}",
                "gamma_gamma_partner_ev": f"{egg:.8e}",
                "normalized_n_epsilon": f"{nd:.8e}",
            }
        )
    for t, s, tr in zip(tau, survival, cell_transfer):
        rows.append(
            {
                "kind": "survival",
                "target_photon_ev": f"{t:.8e}",
                "pgamma_gamma_p_threshold": f"{s:.8e}",
                "bh_gamma_p_threshold": f"{tr:.8e}",
                "gamma_gamma_partner_ev": "",
                "normalized_n_epsilon": "",
            }
        )
    write_rows(DATA_DIR / "figA1_hadronic_thresholds.csv", rows)

    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.45))
    fig.subplots_adjust(left=0.08, right=0.985, bottom=0.24, top=0.84, wspace=0.58)

    axes[0].loglog(eps_ev, photon_density, color=BLUE, lw=1.7)
    add_panel(axes[0], "a")
    axes[0].set_xlabel(r"target photon energy $\epsilon'$ (eV)")
    axes[0].set_ylabel(r"normalized $n'_{\epsilon'}$")
    axes[0].set_title(r"$L_{\nu}\rightarrow u_\nu\rightarrow n_\epsilon$")
    axes[0].grid(color=LIGHT, lw=0.5, which="both")

    axes[1].loglog(eps_ev, gamma_p_pg, color=RED, lw=1.6, label=r"$p\gamma$")
    axes[1].loglog(eps_ev, gamma_p_bh, color=TEAL, lw=1.6, label="BH")
    axes[1].loglog(eps_ev, e_gamma_gg_ev / M_P_C2_EV, color=VIOLET, lw=1.4, label=r"$E_\gamma/m_pc^2$")
    add_panel(axes[1], "b")
    axes[1].set_xlabel(r"$\epsilon'$ (eV)")
    axes[1].set_ylabel("threshold scale")
    axes[1].set_title("Interaction thresholds")
    axes[1].legend(fontsize=5.8)
    axes[1].grid(color=LIGHT, lw=0.5, which="both")

    axes[2].loglog(tau, survival, color=RED, lw=1.7, label=r"$e^{-\tau}$")
    axes[2].loglog(tau, cell_transfer, color=BLUE, lw=1.7, label=r"$(1-e^{-\tau})/\tau$")
    add_panel(axes[2], "c")
    axes[2].set_xlabel(r"optical depth $\tau$")
    axes[2].set_ylabel("survival / transfer")
    axes[2].set_ylim(ATTENUATION_DISPLAY_MIN, 1.2)
    axes[2].set_title("Photon sink semantics")
    axes[2].legend(fontsize=6)
    axes[2].grid(color=LIGHT, lw=0.5, which="both")
    save_pub(fig, "figA1_hadronic_feedback")


def fig4_reverse_shock() -> None:
    sigma = read_csv(DATA_DIR / "fig4_reverse_shock_sigma_summary.csv")
    lc = read_csv(DATA_DIR / "fig4_reverse_shock_lightcurve_summary.csv")
    events = read_csv(DATA_DIR / "fig4_secondary_rs_events.csv")
    energy = read_csv(DATA_DIR / "fig4_secondary_rs_energy.csv")

    sig = np.array([float(r["sigma"]) for r in sigma])
    sigma_floor = np.nanmin(sig[sig > 0.0])
    sig_plot = np.where(sig > 0.0, sig, sigma_floor / 2.0)
    max_b = np.array([float(r["max_B3_G"]) for r in sigma])
    gamma34 = np.array([float(r["max_gamma34"]) for r in sigma])
    optical = [r for r in lc if float(r["nu_hz"]) == 1.0e14]
    sig_lc = np.array([float(r["sigma"]) for r in optical])
    sigma_floor_lc = np.nanmin(sig_lc[sig_lc > 0.0])
    sig_lc_plot = np.where(sig_lc > 0.0, sig_lc, sigma_floor_lc / 2.0)
    rs_frac = np.array([float(r["rs_to_total_at_total_peak"]) for r in optical])
    start = np.array([float(r["start_tobs_axis_s"]) for r in events])
    end = np.array([float(r["shock_end_tobs_axis_s"]) for r in events])
    inj = float(energy[0]["secondary_rs_electron_injected_energy_erg"])
    diss = float(energy[0]["secondary_rs_dissipated_energy_erg"])

    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.35))
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.23, top=0.84, wspace=0.70)
    sigma_ticks = [sigma_floor / 2.0, 1.0e-2, 1.0, 1.0e2]
    sigma_ticklabels = ["0", r"$10^{-2}$", "1", r"$10^2$"]
    axes[0].plot(sig_plot, max_b, color=BLUE, marker="o", label=r"$B_3$")
    axes[0].set_xscale("log")
    axes[0].set_xticks(sigma_ticks)
    axes[0].set_xticklabels(sigma_ticklabels)
    ax2 = axes[0].twinx()
    ax2.plot(sig_plot, gamma34, color=TEAL, marker="s", label=r"$\gamma_{34}$")
    ax2.set_xscale("log")
    ax2.set_xticks(sigma_ticks)
    add_panel(axes[0], "a")
    axes[0].set_xlabel(r"$\sigma$")
    axes[0].set_ylabel(r"max $B_3$ (G)", color=BLUE)
    ax2.set_ylabel(r"max $\gamma_{34}$", color=TEAL, labelpad=1)
    axes[0].set_title("Magnetized RS state")
    axes[0].grid(color=LIGHT, which="both", lw=0.5)

    axes[1].plot(sig_lc_plot, rs_frac, color=RED, marker="o")
    axes[1].set_xscale("log")
    axes[1].set_xticks(sigma_ticks)
    axes[1].set_xticklabels(sigma_ticklabels)
    add_panel(axes[1], "b")
    axes[1].set_xlabel(r"$\sigma$")
    axes[1].set_ylabel("RS fraction at peak", labelpad=2)
    axes[1].set_ylim(0, 1.08)
    axes[1].set_title("RS contribution")
    axes[1].grid(color=LIGHT, which="both", lw=0.5)

    y = np.arange(start.size)
    axes[2].barh(y, end - start, left=start, color=VIOLET, edgecolor="black", linewidth=0.6)
    axes[2].set_xscale("log")
    axes[2].set_yticks(y)
    axes[2].set_yticklabels([f"jump {r['jump_index']}" for r in events])
    axes[2].set_xlabel(r"$t_{\rm obs}$ (s)")
    axes[2].set_title("Secondary RS windows")
    axes[2].text(0.52, 0.08, rf"$E_e/E_{{\rm diss}}={inj/diss:.2f}$", transform=axes[2].transAxes, ha="center")
    add_panel(axes[2], "c")
    axes[2].grid(axis="x", color=LIGHT, lw=0.5, which="both")
    save_pub(fig, "fig4_reverse_shock")


def fig5_validation_benchmark() -> None:
    rows = read_csv(DATA_DIR / "fig5_validation_benchmark.csv")
    dynamics = [row for row in rows if row["kind"] == "dynamics_event_split"]
    angular = [row for row in rows if row["kind"] == "angular_projection"]
    rs = [row for row in rows if row["kind"] == "reverse_shock_adaptive"]
    transport = [row for row in rows if row["kind"] == "transport_remap"]

    dyn_tab = [row for row in dynamics if row["case"] == "tabulated_csm"]
    dyn_jump = [row for row in dynamics if row["case"] == "multi_jump"]
    nr_tab = np.array([float(row["coordinate"]) for row in dyn_tab])
    nr_jump = np.array([float(row["coordinate"]) for row in dyn_jump])
    gamma_tab = np.array([float(row["metric_a"]) for row in dyn_tab])
    mass_tab = np.array([float(row["metric_b"]) for row in dyn_tab])
    gamma_jump = np.array([float(row["metric_a"]) for row in dyn_jump])
    mass_jump = np.array([float(row["metric_b"]) for row in dyn_jump])

    theta_values = np.array([float(row["coordinate"]) for row in angular])
    angular_p95 = np.array([float(row["metric_b"]) for row in angular])
    angular_median = np.array([float(row["metric_c"]) for row in angular])

    rs_freq = np.array([float(row["case"]) for row in rs])
    peak_flux = np.array([float(row["metric_b"]) for row in rs])
    integral = np.array([float(row["metric_c"]) for row in rs])
    freq_labels = ["1 GHz" if nu == 1.0e9 else "opt." if nu == 1.0e14 else "X-ray" for nu in rs_freq]

    transport_dims = [row["case"] for row in transport]
    transport_max = [float(row["metric_a"]) for row in transport]
    transport_p95 = [float(row["metric_b"]) for row in transport]

    fig, axes = plt.subplots(2, 2, figsize=(7.1, 4.75))
    fig.subplots_adjust(left=0.08, right=0.985, bottom=0.12, top=0.92, wspace=0.42, hspace=0.48)
    ax = axes[0, 0]

    ax.loglog(nr_tab, gamma_tab, color=BLUE, marker="o", lw=1.4, label=r"$\Gamma$, tabulated")
    ax.loglog(nr_tab, mass_tab, color=TEAL, marker="s", lw=1.4, label=r"$M_{\rm sw}$, tabulated")
    ax.loglog(nr_jump, gamma_jump, color=RED, marker="o", lw=1.1, ls=":", label=r"$\Gamma$, jumps")
    ax.loglog(nr_jump, mass_jump, color=VIOLET, marker="s", lw=1.1, ls=":", label=r"$M_{\rm sw}$, jumps")
    add_panel(ax, "a")
    ax.set_xlabel(r"$N_R$")
    ax.set_ylabel("relative error")
    ax.set_title("Dynamics event split")
    ax.legend(fontsize=5.2)
    ax.grid(color=LIGHT, lw=0.5, which="both")

    ax = axes[0, 1]
    ax.plot(theta_values, 100.0 * angular_p95, color=BLUE, marker="o", lw=1.4, label="p95")
    ax.plot(theta_values, 100.0 * angular_median, color=TEAL, marker="s", lw=1.4, label="median")
    add_panel(ax, "b")
    ax.set_xlabel(r"$\theta_{\rm obs}/\theta_{\rm c}$")
    ax.set_ylabel("|relative error| (%)")
    ax.set_title("Projection vs 300x48")
    ax.legend(fontsize=5.8)
    ax.grid(color=LIGHT, lw=0.5)

    ax = axes[1, 0]
    y = np.arange(len(freq_labels))
    ax.barh(y - 0.18, 100.0 * peak_flux, height=0.32, color=RED, label="peak flux")
    ax.barh(y + 0.18, 100.0 * integral, height=0.32, color=VIOLET, label="integral")
    ax.axvline(0.0, color=NEUTRAL, lw=0.7)
    add_panel(ax, "c")
    ax.set_yticks(y)
    ax.set_yticklabels(freq_labels)
    ax.set_xlabel("fractional difference (%)")
    ax.set_title("RS adaptive vs direct")
    ax.legend(fontsize=5.8)
    ax.grid(axis="x", color=LIGHT, lw=0.5)

    ax = axes[1, 1]
    x = np.arange(len(transport_dims))
    ax.bar(x - 0.18, 100.0 * np.array(transport_max), width=0.32, color=GOLD, label="max flux")
    ax.bar(x + 0.18, 100.0 * np.array(transport_p95), width=0.32, color=GREEN, label="p95 flux")
    add_panel(ax, "d")
    ax.set_xticks(x)
    ax.set_xticklabels(["1D", "2D"])
    ax.set_ylabel("|before-after| (%)")
    ax.set_title("Transport remap")
    ax.legend(fontsize=5.8)
    ax.grid(axis="y", color=LIGHT, lw=0.5)
    save_pub(fig, "fig5_validation_benchmark")


def main() -> None:
    ensure_dirs()
    fig1_radius_state()
    fig2_forward_spectrum()
    fig3_electron_transport()
    fig4_reverse_shock()
    fig5_validation_benchmark()
    figA1_hadronic_thresholds()


if __name__ == "__main__":
    main()
