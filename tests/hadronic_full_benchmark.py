"""ASGARD 全量强子过程 benchmark：运行所有强子通道并绘制诊断图。"""
from __future__ import annotations

from pathlib import Path
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
THIS_DIR = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet

OUTPUT_DIR = ROOT / "output" / "hadronic_full_benchmark"
PROCESS_LABELS = ("photopion", "pion_decay", "muon_decay")
CHANNEL_COLORS = {
    "proton_synch": "#1f77b4",
    "pg_gamma": "#ff7f0e",
    "bethe_heitler": "#2ca02c",
    "hadronic_ic": "#d62728",
    "pion_synch": "#9467bd",
    "muon_synch": "#8c564b",
    "pion_ic": "#e377c2",
    "muon_ic": "#7f7f7f",
}


def _save_fig(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".pdf"), dpi=200, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def _build_full_model() -> Model:
    return Model(
        TophatJet(0.1, 1.0e52, 300.0),
        ISM(1.0),
        Observer(1.0e26, 0.1, 0.0),
        Radiation(
            0.1, 1.0e-3, 2.3,
            epsilon_p=0.2, ssc=True,
            proton_synch=True, pg=True, neutrino=True,
            bethe_heitler=True, hadronic_inverse_compton=True, pp=True,
            pgamma_scheme="hummer_2010_response",
        ),
        setups=Setups(
            electron_solver="fullhide_1d",
            hadronic_enabled=True, hadronic_solver="am3_1d",
            pgamma_scheme="hummer_2010_response",
            num_gam_e=24, num_gam_p=40, num_nu=40, num_nu_nu=24,
            num_r=24, num_theta=16, num_tobs=24,
        ),
    )


def plot_proton_spectrum(details, output_dir: Path) -> None:
    fwd = details.fwd
    if fwd.gamma_p is None or fwd.dN_dgamma_p is None:
        print("  skip proton_spectrum: no proton data")
        return
    gam_p = fwd.gamma_p
    dN = fwd.dN_dgamma_p
    num_r = dN.shape[1]
    shell_indices = [0, max(0, num_r // 3), num_r // 2, max(0, 2 * num_r // 3), num_r - 1]

    fig, ax = plt.subplots(figsize=(8, 5))
    for idx in shell_indices:
        ax.loglog(gam_p, np.maximum(dN[:, idx], 1e-99), label=f"shell {idx}")
    ax.set_xlabel(r"proton Lorentz factor $\gamma_p$")
    ax.set_ylabel("$dN/d\\gamma_p$")
    ax.set_title("Proton Spectrum Evolution")
    ax.grid(True, which="both", alpha=0.25, linestyle=":")
    ax.legend(fontsize=8)
    _save_fig(fig, output_dir / "proton_spectrum")


def plot_hadronic_sed(details, output_dir: Path) -> None:
    fwd = details.fwd
    freq = fwd.seed_frequency_hz
    if freq is None:
        print("  skip hadronic_sed: no seed_frequency_hz")
        return
    num_r = fwd.l_had_syn_spec.shape[1] if fwd.l_had_syn_spec is not None else 0
    shell_indices = [0, max(0, num_r // 2), max(0, num_r - 1)]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)

    # Proton synchrotron
    ax = axes[0, 0]
    for idx in shell_indices:
        if fwd.l_had_syn_spec is not None:
            ax.loglog(freq, np.maximum(fwd.l_had_syn_spec[:, idx], 1e-99), label=f"shell {idx}")
    ax.set_xlabel("$\\nu$ [Hz]"); ax.set_ylabel("$\\nu L_\\nu$ [erg/s]")
    ax.set_title("Proton Synchrotron"); ax.grid(True, which="both", alpha=0.25, linestyle=":")
    ax.legend(fontsize=8)

    # P-gamma
    ax = axes[0, 1]
    for idx in shell_indices:
        if fwd.l_had_pg_gamma_spec is not None:
            ax.loglog(freq, np.maximum(fwd.l_had_pg_gamma_spec[:, idx], 1e-99), label=f"shell {idx}")
    ax.set_xlabel("$\\nu$ [Hz]"); ax.set_ylabel("$\\nu L_\\nu$ [erg/s]")
    ax.set_title("P-Gamma Photons"); ax.grid(True, which="both", alpha=0.25, linestyle=":")
    ax.legend(fontsize=8)

    # Bethe-Heitler
    ax = axes[1, 0]
    for idx in shell_indices:
        if fwd.l_had_bethe_heitler_spec is not None:
            ax.loglog(freq, np.maximum(fwd.l_had_bethe_heitler_spec[:, idx], 1e-99), label=f"shell {idx}")
    ax.set_xlabel("$\\nu$ [Hz]"); ax.set_ylabel("$\\nu L_\\nu$ [erg/s]")
    ax.set_title("Bethe-Heitler Pair Synchrotron"); ax.grid(True, which="both", alpha=0.25, linestyle=":")
    ax.legend(fontsize=8)

    # Hadronic IC
    ax = axes[1, 1]
    for idx in shell_indices:
        if fwd.l_had_hadronic_ic_spec is not None:
            ax.loglog(freq, np.maximum(fwd.l_had_hadronic_ic_spec[:, idx], 1e-99), label=f"shell {idx}")
    ax.set_xlabel("$\\nu$ [Hz]"); ax.set_ylabel("$\\nu L_\\nu$ [erg/s]")
    ax.set_title("Hadronic Inverse Compton"); ax.grid(True, which="both", alpha=0.25, linestyle=":")
    ax.legend(fontsize=8)

    _save_fig(fig, output_dir / "hadronic_spectral_components")


def plot_secondary_radiation(details, output_dir: Path) -> None:
    fwd = details.fwd
    freq = fwd.seed_frequency_hz
    num_r = fwd.hadronic_pion_synch.shape[1] if fwd.hadronic_pion_synch is not None else 0
    idx = max(0, num_r // 2)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)
    for ax, data, name in [
        (axes[0, 0], fwd.hadronic_pion_synch, "Pion Synchrotron"),
        (axes[0, 1], fwd.hadronic_muon_synch, "Muon Synchrotron"),
        (axes[1, 0], fwd.hadronic_pion_inverse_compton, "Pion IC"),
        (axes[1, 1], fwd.hadronic_muon_inverse_compton, "Muon IC"),
    ]:
        if data is not None and freq is not None:
            ax.loglog(freq, np.maximum(data[:, idx], 1e-99))
        ax.set_xlabel("$\\nu$ [Hz]"); ax.set_ylabel("$\\nu L_\\nu$ [erg/s]")
        ax.set_title(name); ax.grid(True, which="both", alpha=0.25, linestyle=":")

    _save_fig(fig, output_dir / "secondary_radiation")


def plot_species_distributions(details, output_dir: Path) -> None:
    fwd = details.fwd
    gam_sec = fwd.gamma_secondary
    num_r = fwd.dN_dgamma_n.shape[1] if fwd.dN_dgamma_n is not None else 0
    idx = max(0, num_r // 2)

    fig, axes = plt.subplots(2, 4, figsize=(20, 10), constrained_layout=True)
    species = [
        (axes[0, 0], fwd.dN_dgamma_n, "Neutron"),
        (axes[0, 1], fwd.dN_dgamma_pi_plus, "Pion $\\pi^+$"),
        (axes[0, 2], fwd.dN_dgamma_pi_minus, "Pion $\\pi^-$"),
        (axes[0, 3], fwd.dN_dgamma_p, "Proton (ref)"),
        (axes[1, 0], fwd.dN_dgamma_mu_plus_left, "Muon $\\mu^+_L$"),
        (axes[1, 1], fwd.dN_dgamma_mu_plus_right, "Muon $\\mu^+_R$"),
        (axes[1, 2], fwd.dN_dgamma_mu_minus_left, "Muon $\\mu^-_L$"),
        (axes[1, 3], fwd.dN_dgamma_mu_minus_right, "Muon $\\mu^-_R$"),
    ]
    for ax, data, name in species:
        if data is not None and gam_sec is not None:
            ax.loglog(gam_sec, np.maximum(data[:, idx], 1e-99))
        ax.set_xlabel("$\\gamma$"); ax.set_ylabel("$dN/d\\gamma$")
        ax.set_title(name); ax.grid(True, which="both", alpha=0.25, linestyle=":")

    _save_fig(fig, output_dir / "species_distributions")


def plot_process_power(details, output_dir: Path) -> None:
    fwd = details.fwd
    if fwd.am3_process_power is None:
        print("  skip process_power: no AM3 data")
        return
    num_r = fwd.am3_process_power.shape[2]
    idx = max(0, num_r // 2)

    fig, ax = plt.subplots(figsize=(10, 6))
    for i_proc, label in enumerate(PROCESS_LABELS):
        gam_p = np.arange(fwd.am3_process_power.shape[1])
        power = np.maximum(fwd.am3_process_power[i_proc, :, idx], 1e-99)
        ax.semilogy(gam_p, power, label=label)
    ax.set_xlabel("proton grid index"); ax.set_ylabel("power [erg/s]")
    ax.set_title("AM3 Process Power Breakdown"); ax.grid(True, which="both", alpha=0.25, linestyle=":")
    ax.legend()
    _save_fig(fig, output_dir / "process_power")


def plot_neutrino_luminosity(details, output_dir: Path) -> None:
    fwd = details.fwd
    if fwd.neutrino_luminosity is None or fwd.neutrino_frequency_hz is None:
        print("  skip neutrino: no data")
        return
    nu_freq = fwd.neutrino_frequency_hz
    nu_lum = fwd.neutrino_luminosity
    num_r = nu_lum.shape[1]
    shell_indices = [0, max(0, num_r // 2), max(0, num_r - 1)]

    fig, ax = plt.subplots(figsize=(8, 5))
    for idx in shell_indices:
        ax.loglog(nu_freq, np.maximum(nu_lum[:, idx], 1e-99), label=f"shell {idx}")
    ax.set_xlabel("$\\nu$ [Hz]"); ax.set_ylabel("$\\nu L_\\nu$ [erg/s]")
    ax.set_title("Neutrino Luminosity"); ax.grid(True, which="both", alpha=0.25, linestyle=":")
    ax.legend(fontsize=8)
    _save_fig(fig, output_dir / "neutrino_luminosity")


def plot_timing_breakdown(details, output_dir: Path) -> None:
    fwd = details.fwd
    timings = fwd.timings if hasattr(fwd, 'timings') else {}
    if not timings or timings.get("total", 0.0) <= 0.0:
        print("  skip timing: no data")
        return

    entries = sorted(
        [(k, v) for k, v in timings.items() if k != "total"],
        key=lambda x: -x[1],
    )
    if not entries:
        return
    labels = [e[0] for e in entries]
    values = [e[1] for e in entries]

    fig, ax = plt.subplots(figsize=(10, 5))
    bars = ax.bar(labels, values, color=[CHANNEL_COLORS.get(l, "#888888") for l in labels])
    ax.set_ylabel("wall-clock time [s]")
    ax.set_title(f"Hadronic Process Timing Breakdown (total={timings['total']:.2f}s)")
    ax.tick_params(axis="x", rotation=45)
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() * 1.01,
                f"{val:.2f}s", ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    _save_fig(fig, output_dir / "timing_breakdown")


def plot_tau_pg(details, output_dir: Path) -> None:
    fwd = details.fwd
    if fwd.tau_pg is None or fwd.pg_photon_survival is None:
        print("  skip tau_pg: no data")
        return
    freq = fwd.seed_frequency_hz
    num_r = fwd.tau_pg.shape[1]
    idx = max(0, num_r // 2)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    if freq is not None:
        ax1.loglog(freq, np.maximum(fwd.tau_pg[:, idx], 1e-99))
    ax1.set_xlabel("$\\nu$ [Hz]"); ax1.set_ylabel("$\\tau_{p\\gamma}$")
    ax1.set_title("P-Gamma Optical Depth"); ax1.grid(True, which="both", alpha=0.25, linestyle=":")

    if freq is not None:
        ax2.plot(freq, fwd.pg_photon_survival[:, idx])
    ax2.set_xscale("log"); ax2.set_xlabel("$\\nu$ [Hz]"); ax2.set_ylabel("survival factor")
    ax2.set_title("Photon Survival Factor"); ax2.grid(True, which="both", alpha=0.25, linestyle=":")

    _save_fig(fig, output_dir / "tau_pg")


def main() -> None:
    output_dir = OUTPUT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)
    print("Building full hadronic model...")
    model = _build_full_model()
    print(f"  num_gam_p={model.setups.num_gam_p}, num_r={model.setups.num_r}")

    t_start = time.perf_counter()
    print("Running model.details(1e3, 1e6)...")
    details = model.details(1.0e3, 1.0e6)
    t_elapsed = time.perf_counter() - t_start
    print(f"  total wall time: {t_elapsed:.2f}s")

    fwd = details.fwd
    print(f"  gam_p range: [{fwd.gamma_p[0]:.1f}, {fwd.gamma_p[-1]:.1f}]")
    if fwd.gamma_secondary is not None:
        print(f"  gamma_secondary range: [{fwd.gamma_secondary[0]:.1f}, {fwd.gamma_secondary[-1]:.1f}]")
    if fwd.neutrino_luminosity is not None:
        print(f"  neutrino_luminosity shape: {fwd.neutrino_luminosity.shape}")
    if hasattr(fwd, 'timings') and fwd.timings:
        print(f"  hadronic timings: {', '.join(f'{k}={v:.2f}s' for k, v in sorted(fwd.timings.items()))}")

    print("\nGenerating diagnostic plots...")
    plot_proton_spectrum(details, output_dir)
    plot_hadronic_sed(details, output_dir)
    plot_secondary_radiation(details, output_dir)
    plot_species_distributions(details, output_dir)
    plot_process_power(details, output_dir)
    plot_neutrino_luminosity(details, output_dir)
    plot_timing_breakdown(details, output_dir)
    plot_tau_pg(details, output_dir)

    print(f"\nBenchmark complete. Output saved to {output_dir}")
    for f in sorted(output_dir.iterdir()):
        print(f"  {f.name}")


if __name__ == "__main__":
    main()
