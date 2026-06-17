from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from _repo_path import ensure_repo_root_on_path


ensure_repo_root_on_path()

from tests.public_api_builders import numerics, radiation, top_hat_model


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
MEC2_ERG = 8.18710565e-7
EV_TO_ERG = 1.602176634e-12


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


def diagnostic_model():
    return top_hat_model(
        fwd_rad=radiation(include_ssc=True, include_kn_correction=False),
        numerics=numerics(
            num_radius=24,
            num_theta=12,
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
    theta = np.linspace(0.0, 0.35, 160)
    theta_j = 0.10
    solid_weight = np.sin(theta)
    solid_weight = solid_weight / np.trapezoid(solid_weight, theta)
    selected = np.linspace(0, radius.size - 1, min(6, radius.size), dtype=int)
    t_axis = radius / (2.0 * C_CGS * gamma**2)
    dt_on = (1.0 / beta - 1.0) / C_CGS
    dt_edge = (1.0 / beta - np.cos(theta_j)) / C_CGS

    rows: list[dict[str, object]] = []
    for r, g, b, ta, do, de in zip(radius, gamma, b_comv, t_axis, dt_on, dt_edge):
        rows.append(
            {
                "kind": "asgard_state",
                "radius_cm": f"{r:.8e}",
                "Gamma": f"{g:.8e}",
                "B_comoving_G": f"{b:.8e}",
                "thin_axis_tobs_s": f"{ta:.8e}",
                "dtobs_dR_on_axis_s_per_cm": f"{do:.8e}",
                "dtobs_dR_edge_s_per_cm": f"{de:.8e}",
            }
        )
    for th, weight in zip(theta, solid_weight):
        rows.append(
            {
                "kind": "angular_weight",
                "radius_cm": "",
                "Gamma": "",
                "B_comoving_G": "",
                "thin_axis_tobs_s": f"{th:.8e}",
                "dtobs_dR_on_axis_s_per_cm": f"{weight:.8e}",
                "dtobs_dR_edge_s_per_cm": "",
            }
        )
    write_rows(DATA_DIR / "fig1_radius_state.csv", rows)

    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.35))
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.24, top=0.84, wspace=0.78)

    axes[0].loglog(radius, gamma, color=BLUE, lw=1.7)
    ax0b = axes[0].twinx()
    ax0b.loglog(radius, b_comv, color=TEAL, lw=1.5)
    add_panel(axes[0], "a")
    axes[0].set_xlabel(r"$R$ (cm)")
    axes[0].set_ylabel(r"$\Gamma$", color=BLUE)
    ax0b.set_ylabel(r"$B'$ (G)", color=TEAL, labelpad=8)
    axes[0].set_title("Local blast-wave state")
    axes[0].grid(color=LIGHT, lw=0.5, which="both")

    axes[1].loglog(radius, t_axis, color=NEUTRAL, lw=1.5, label=r"$R/2c\Gamma^2$")
    axes[1].scatter(radius[selected], t_axis[selected], color=RED, s=14, zorder=3)
    add_panel(axes[1], "b")
    axes[1].set_xlabel(r"$R$ (cm)")
    axes[1].set_ylabel(r"$t_{\rm obs}$ scale (s)")
    axes[1].set_title("Observer time is derived")
    axes[1].grid(color=LIGHT, lw=0.5, which="both")

    axes[2].semilogy(radius, dt_on * C_CGS, color=BLUE, lw=1.6, label="on axis")
    axes[2].semilogy(radius, dt_edge * C_CGS, color=VIOLET, lw=1.6, label=r"$\theta=0.1$")
    axes[2].text(-0.16, 1.06, "c", transform=axes[2].transAxes, fontweight="bold", fontsize=9)
    axes[2].set_xlabel(r"$R$ (cm)")
    axes[2].set_ylabel(r"$c\,dt_{\rm obs}/dR$")
    axes[2].set_title("EATS stiffness")
    axes[2].legend(fontsize=6)
    axes[2].grid(color=LIGHT, lw=0.5, which="both")
    save_pub(fig, "fig1_radius_state")


def fig2_forward_spectrum() -> None:
    model = diagnostic_model()
    times = np.logspace(2.0, 7.0, 42)
    freqs = np.array([1.0e9, 1.0e14, 1.0e18])
    spec_freq = np.logspace(7.0, 21.0, 80)
    spec_times = np.array([1.0e3, 1.0e5, 1.0e7])
    flux = model.flux_density_grid(times, freqs)
    spec_flux = model.flux_density_grid(spec_times, spec_freq)

    rows: list[dict[str, object]] = []
    for i, nu in enumerate(freqs):
        for t, total, sync, ssc in zip(times, flux.total[i], flux.fwd.sync[i], flux.fwd.ssc[i]):
            rows.append(
                {
                    "kind": "lightcurve",
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
                    "time_s": f"{t:.8e}",
                    "frequency_hz": f"{nu:.8e}",
                    "total_fnu_cgs": f"{total:.8e}",
                    "fs_synch_fnu_cgs": f"{sync:.8e}",
                    "fs_ssc_fnu_cgs": f"{ssc:.8e}",
                }
            )
    write_rows(DATA_DIR / "fig2_forward_api.csv", rows)

    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.45))
    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.24, top=0.84, wspace=0.58)
    colors = [BLUE, TEAL, RED]
    labels = ["1 GHz", "optical", "X-ray"]
    for i, (color, label) in enumerate(zip(colors, labels)):
        axes[0].loglog(times, flux.total[i], color=color, lw=1.6, label=label)
    add_panel(axes[0], "a")
    axes[0].set_xlabel(r"$t_{\rm obs}$ (s)")
    axes[0].set_ylabel(r"$F_\nu$ (cgs)")
    axes[0].set_title("Multi-frequency output")
    axes[0].legend(fontsize=6)
    axes[0].grid(color=LIGHT, lw=0.5, which="both")

    for j, (t, color) in enumerate(zip(spec_times, [BLUE, TEAL, RED])):
        axes[1].loglog(spec_freq, spec_flux.total[:, j], color=color, lw=1.5, label=rf"$10^{{{int(np.log10(t))}}}$ s")
    add_panel(axes[1], "b")
    axes[1].set_xlabel(r"$\nu$ (Hz)")
    axes[1].set_ylabel(r"$F_\nu$ (cgs)")
    axes[1].set_title("Spectral evolution")
    axes[1].legend(fontsize=6)
    axes[1].grid(color=LIGHT, lw=0.5, which="both")

    ssc_fraction = np.divide(
        spec_flux.fwd.ssc[:, 1],
        spec_flux.total[:, 1],
        out=np.zeros_like(spec_flux.fwd.ssc[:, 1]),
        where=spec_flux.total[:, 1] > 0.0,
    )
    axes[2].semilogx(spec_freq, ssc_fraction, color=VIOLET, lw=1.7)
    add_panel(axes[2], "c")
    axes[2].set_xlabel(r"$\nu$ (Hz)")
    axes[2].set_ylabel("SSC / total")
    axes[2].set_ylim(0.0, min(1.0, float(np.nanmax(ssc_fraction)) * 1.2 + 0.02))
    axes[2].set_title(r"Component role at $10^5$ s")
    axes[2].grid(color=LIGHT, lw=0.5, which="both")
    save_pub(fig, "fig2_forward_api")


def fig3_electron_transport() -> None:
    model = top_hat_model(
        fwd_rad=radiation(include_ssc=False),
        numerics=numerics(
            num_radius=12,
            num_theta=8,
            num_observer_time=12,
            num_electron_gamma=16,
            num_photon_frequency=16,
            num_threads=1,
        ),
    )
    details = model.details(1.0e4, 2.0e4).fwd
    radius = np.asarray(details.radius, dtype=float)
    gamma_e = np.asarray(details.gamma_e, dtype=float)
    electron = np.asarray(details.dN_dgamma_e, dtype=float)
    xgrid = np.log10(gamma_e)
    dx = np.gradient(xgrid)
    number_per_x = electron * gamma_e[:, None] * np.log(10.0)
    cell_number = number_per_x * dx[:, None]
    energy_per_shell = np.sum((gamma_e[:, None] - 1.0) * cell_number, axis=0) * MEC2_ERG
    total_number = np.sum(cell_number, axis=0)

    rnorm = np.linspace(0.0, 1.0, 200)
    cool_a = 1.4e-4
    char_rows: list[dict[str, object]] = []
    for g0 in (1.0e3, 1.0e4, 1.0e5, 1.0e6):
        g = g0 / (1.0 + cool_a * g0 * rnorm)
        for x, y in zip(rnorm, g):
            char_rows.append({"gamma0": f"{g0:.8e}", "normalized_radius_step": f"{x:.8e}", "gamma": f"{y:.8e}"})
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
                }
            )
    write_rows(DATA_DIR / "fig3_electron_transport.csv", rows)
    write_rows(DATA_DIR / "fig3_electron_characteristics.csv", char_rows)

    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.45))
    fig.subplots_adjust(left=0.08, right=0.985, bottom=0.24, top=0.84, wspace=0.58)
    pcm = axes[0].pcolormesh(
        np.log10(radius),
        np.log10(gamma_e),
        np.log10(np.maximum(electron, np.max(electron) * 1.0e-18)),
        shading="auto",
        cmap="viridis",
    )
    add_panel(axes[0], "a")
    axes[0].set_xlabel(r"$\log_{10} R$ (cm)")
    axes[0].set_ylabel(r"$\log_{10}\gamma_e$")
    axes[0].set_title(r"$N_e(\gamma,R)$")
    fig.colorbar(pcm, ax=axes[0], fraction=0.045, pad=0.02, label=r"$\log N_e$")

    for idx, color in zip(np.linspace(0, radius.size - 1, 4, dtype=int), [BLUE, TEAL, RED, VIOLET]):
        spec = gamma_e**2 * np.maximum(electron[:, idx], 0.0)
        norm = np.max(spec)
        spec = spec / norm if norm > 0.0 else spec
        axes[1].loglog(gamma_e, np.maximum(spec, 1.0e-10), color=color, lw=1.4, label=rf"$R={radius[idx]:.1e}$")
    add_panel(axes[1], "b")
    axes[1].set_xlabel(r"$\gamma_e$")
    axes[1].set_ylabel(r"normalized $\gamma_e^2 dN/d\gamma_e$")
    axes[1].set_ylim(1.0e-6, 2.0)
    axes[1].set_title("Spectral aging")
    axes[1].legend(fontsize=5.4)
    axes[1].grid(color=LIGHT, lw=0.5, which="both")

    for g0, color in zip((1.0e3, 1.0e4, 1.0e5, 1.0e6), [BLUE, TEAL, RED, VIOLET]):
        g = g0 / (1.0 + cool_a * g0 * rnorm)
        axes[2].semilogy(rnorm, g, color=color, lw=1.4)
        axes[2].text(rnorm[-1], g[-1], rf"$10^{{{int(np.log10(g0))}}}$", color=color, fontsize=6)
    add_panel(axes[2], "c")
    axes[2].set_xlabel(r"normalized $\Delta R$")
    axes[2].set_ylabel(r"cooling characteristic $\gamma_e$")
    axes[2].set_title(r"$d\gamma/dR\propto-\gamma^2$")
    axes[2].grid(color=LIGHT, lw=0.5, which="both")
    save_pub(fig, "fig3_transport_projection")


def fig4_hadronic_thresholds() -> None:
    eps_ev = np.logspace(-6.0, 9.0, 240)
    eps_erg = eps_ev * EV_TO_ERG
    # Approximate threshold products in the comoving frame.
    pgamma_threshold_ev2 = 0.15e18
    bh_threshold_ev2 = 2.0 * 0.511e6 * 0.938e9
    gamma_p_pg = pgamma_threshold_ev2 / eps_ev / 0.938e9
    gamma_p_bh = bh_threshold_ev2 / eps_ev / 0.938e9
    e_gamma_gg_ev = (0.511e6) ** 2 / eps_ev
    tau = np.logspace(-3.0, 2.0, 220)
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
    write_rows(DATA_DIR / "fig4_hadronic_thresholds.csv", rows)

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
    axes[1].loglog(eps_ev, e_gamma_gg_ev / 1.0e9, color=VIOLET, lw=1.4, label=r"$\gamma\gamma$ partner (GeV)")
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
    axes[2].set_title("Photon sink semantics")
    axes[2].legend(fontsize=6)
    axes[2].grid(color=LIGHT, lw=0.5, which="both")
    save_pub(fig, "fig4_hadronic_feedback")


def fig5_reverse_shock() -> None:
    sigma = read_csv(ROOT / "output/asgard_doc/magnetized_rs_sigma_benchmark/sigma_scan_summary.csv")
    lc = read_csv(ROOT / "output/asgard_doc/magnetized_rs_sigma_benchmark/sigma_scan_lightcurve_summary.csv")
    events = read_csv(ROOT / "output/asgard_doc/reverse_density_jump_tests/triple_density_jump_rs_fs_tophat_secondary_rs_events.csv")
    energy = read_csv(ROOT / "output/asgard_doc/reverse_density_jump_tests/triple_density_jump_rs_fs_tophat_secondary_rs_energy.csv")
    write_rows(DATA_DIR / "fig5_reverse_shock_sigma_summary.csv", sigma)
    write_rows(DATA_DIR / "fig5_reverse_shock_lightcurve_summary.csv", lc)
    write_rows(DATA_DIR / "fig5_secondary_rs_events.csv", events)
    write_rows(DATA_DIR / "fig5_secondary_rs_energy.csv", energy)

    sig = np.array([float(r["sigma"]) for r in sigma])
    max_b = np.array([float(r["max_B3_G"]) for r in sigma])
    gamma34 = np.array([float(r["max_gamma34"]) for r in sigma])
    optical = [r for r in lc if float(r["nu_hz"]) == 1.0e14]
    sig_lc = np.array([float(r["sigma"]) for r in optical])
    rs_frac = np.array([float(r["rs_to_total_at_total_peak"]) for r in optical])
    start = np.array([float(r["start_tobs_axis_s"]) for r in events])
    end = np.array([float(r["shock_end_tobs_axis_s"]) for r in events])
    inj = float(energy[0]["secondary_rs_electron_injected_energy_erg"])
    diss = float(energy[0]["secondary_rs_dissipated_energy_erg"])

    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.35))
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.23, top=0.84, wspace=0.70)
    axes[0].loglog(sig + 1.0e-6, max_b, color=BLUE, marker="o", label=r"$B_3$")
    ax2 = axes[0].twinx()
    ax2.semilogx(sig + 1.0e-6, gamma34, color=TEAL, marker="s", label=r"$\gamma_{34}$")
    add_panel(axes[0], "a")
    axes[0].set_xlabel(r"$\sigma + 10^{-6}$")
    axes[0].set_ylabel(r"max $B_3$ (G)", color=BLUE)
    ax2.set_ylabel(r"max $\gamma_{34}$", color=TEAL, labelpad=1)
    axes[0].set_title("Magnetized RS state")
    axes[0].grid(color=LIGHT, which="both", lw=0.5)

    axes[1].semilogx(sig_lc + 1.0e-6, rs_frac, color=RED, marker="o")
    add_panel(axes[1], "b")
    axes[1].set_xlabel(r"$\sigma + 10^{-6}$")
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
    save_pub(fig, "fig5_reverse_shock")


def fig6_validation_benchmark() -> None:
    dynamics = read_csv(ROOT / "output/asgard_doc/dynamics_event_split_tests/multi_jump_tabulated_before_after_delta_metrics.csv")
    angular = read_csv(ROOT / "output/asgard_doc/angular_sampling_compare/fullref300x48_vs_dominant-region-ioka-time-v1_20x15_metrics.csv")
    rs = read_csv(ROOT / "output/asgard_doc/reverse_density_jump_tests/triple_density_jump_rs_fs_tophat_adaptive_convergence.csv")
    before = np.load(ROOT / "output/asgard_doc/charint_structure_preserving/before_charint_ppm_clip.npz")
    after = np.load(ROOT / "output/asgard_doc/charint_structure_preserving/after_pfc_ppm_charint.npz")

    rows: list[dict[str, object]] = []
    for row in dynamics:
        rows.append(
            {
                "kind": "dynamics_event_split",
                "case": row["case"],
                "coordinate": row["num_r"],
                "metric_a": row["max_gamma_direct_rel"],
                "metric_b": row["max_swept_mass_direct_rel"],
                "metric_c": row["rms_gamma_direct_rel"],
                "metric_d": row["rms_swept_mass_direct_rel"],
            }
        )
    for row in angular:
        rows.append(
            {
                "kind": "angular_projection",
                "case": "dominant_region_20x15_vs_fullref_300x48",
                "coordinate": row["theta_obs_over_theta_c"],
                "metric_a": row["max_abs"],
                "metric_b": row["p95_abs"],
                "metric_c": row["median_abs"],
                "metric_d": "",
            }
        )
    for row in rs:
        rows.append(
            {
                "kind": "reverse_shock_adaptive",
                "case": row["band_hz"],
                "coordinate": row["adaptive_time_count"],
                "metric_a": row["peak_time_ratio"],
                "metric_b": row["peak_flux_fractional_difference"],
                "metric_c": row["integral_fractional_difference"],
                "metric_d": row["user_time_count"],
            }
        )
    for dim in ("1d", "2d"):
        rel = np.abs((after[f"charint_{dim}_flux"] - before[f"charint_{dim}_flux"]) / before[f"charint_{dim}_flux"])
        number_rel = (
            after[f"charint_{dim}_electron_number"][0] - before[f"charint_{dim}_electron_number"][0]
        ) / before[f"charint_{dim}_electron_number"][0]
        rows.append(
            {
                "kind": "transport_remap",
                "case": dim,
                "coordinate": "all_frequencies",
                "metric_a": f"{np.nanmax(rel):.8e}",
                "metric_b": f"{np.nanpercentile(rel, 95.0):.8e}",
                "metric_c": f"{np.nanmedian(rel):.8e}",
                "metric_d": f"{number_rel:.8e}",
            }
        )
    write_rows(DATA_DIR / "fig6_validation_benchmark.csv", rows)

    dyn_tab = [row for row in dynamics if row["case"] == "tabulated_csm"]
    dyn_jump = [row for row in dynamics if row["case"] == "multi_jump"]
    nr_tab = np.array([float(row["num_r"]) for row in dyn_tab])
    nr_jump = np.array([float(row["num_r"]) for row in dyn_jump])
    gamma_tab = np.array([float(row["max_gamma_direct_rel"]) for row in dyn_tab])
    mass_tab = np.array([float(row["max_swept_mass_direct_rel"]) for row in dyn_tab])
    gamma_jump = np.array([float(row["max_gamma_direct_rel"]) for row in dyn_jump])
    mass_jump = np.array([float(row["max_swept_mass_direct_rel"]) for row in dyn_jump])

    theta_values = np.array([float(row["theta_obs_over_theta_c"]) for row in angular])
    angular_p95 = np.array([float(row["p95_abs"]) for row in angular])
    angular_median = np.array([float(row["median_abs"]) for row in angular])

    rs_freq = np.array([float(row["band_hz"]) for row in rs])
    peak_flux = np.array([float(row["peak_flux_fractional_difference"]) for row in rs])
    integral = np.array([float(row["integral_fractional_difference"]) for row in rs])
    freq_labels = ["1 GHz" if nu == 1.0e9 else "opt." if nu == 1.0e14 else "X-ray" for nu in rs_freq]

    transport_dims = ["1d", "2d"]
    transport_max = []
    transport_p95 = []
    for dim in transport_dims:
        rel = np.abs((after[f"charint_{dim}_flux"] - before[f"charint_{dim}_flux"]) / before[f"charint_{dim}_flux"])
        transport_max.append(np.nanmax(rel))
        transport_p95.append(np.nanpercentile(rel, 95.0))

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
    save_pub(fig, "fig6_validation_benchmark")


def main() -> None:
    ensure_dirs()
    fig1_radius_state()
    fig2_forward_spectrum()
    fig3_electron_transport()
    fig4_hadronic_thresholds()
    fig5_reverse_shock()
    fig6_validation_benchmark()


if __name__ == "__main__":
    main()
