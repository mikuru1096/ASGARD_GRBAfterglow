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

BLUE = "#0F4D92"
BLUE2 = "#3775BA"
TEAL = "#42949E"
RED = "#B64342"
GREEN = "#2E9E44"
VIOLET = "#9A4D8E"
GOLD = "#C7921A"
NEUTRAL = "#4D4D4D"
LIGHT = "#E8E8E8"


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams.update(
    {
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
    for suffix, kwargs in (
        ("svg", {}),
        ("pdf", {}),
        ("tiff", {"dpi": 600}),
    ):
        fig.savefig(FIG_DIR / f"{stem}.{suffix}", bbox_inches="tight", **kwargs)
    plt.close(fig)


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"no rows for {path}")
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


def add_panel(ax: plt.Axes, label: str) -> None:
    ax.text(-0.08, 1.04, label, transform=ax.transAxes, fontweight="bold", fontsize=9)


def draw_box(
    ax: plt.Axes,
    xy: tuple[float, float],
    text: str,
    color: str,
    width: float = 0.32,
    height: float = 0.11,
    fontsize: float = 7.0,
) -> None:
    x, y = xy
    ax.add_patch(
        plt.Rectangle(
            (x - width / 2.0, y - height / 2.0),
            width,
            height,
            facecolor=color,
            edgecolor="black",
            linewidth=0.8,
            alpha=0.16,
        )
    )
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize, linespacing=0.95)


def draw_arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float], color: str = NEUTRAL) -> None:
    ax.annotate("", xy=end, xytext=start, arrowprops=dict(arrowstyle="-|>", color=color, lw=1.0))


def fig1_architecture() -> None:
    nodes = [
        ("Public API", 0.18, 0.80, BLUE),
        ("RuntimeConfig", 0.40, 0.80, BLUE),
        ("SolveState", 0.62, 0.80, BLUE),
        ("Observer products", 0.84, 0.80, BLUE),
        ("Dynamics", 0.25, 0.55, TEAL),
        ("Electron transport", 0.45, 0.55, TEAL),
        ("Photon fields", 0.65, 0.55, TEAL),
        ("Hadronic transport", 0.85, 0.55, TEAL),
        ("Reverse shock", 0.35, 0.30, VIOLET),
        ("Polarization", 0.58, 0.30, VIOLET),
        ("EATS projection", 0.78, 0.30, VIOLET),
    ]
    edges = [
        ("Public API", "RuntimeConfig"),
        ("RuntimeConfig", "SolveState"),
        ("SolveState", "Observer products"),
        ("Dynamics", "Electron transport"),
        ("Electron transport", "Photon fields"),
        ("Photon fields", "Hadronic transport"),
        ("Reverse shock", "Polarization"),
        ("Polarization", "EATS projection"),
        ("EATS projection", "Observer products"),
        ("Hadronic transport", "Observer products"),
    ]
    node_pos = {name: (x, y) for name, x, y, _ in nodes}
    node_labels = {
        "Public API": "Public\nAPI",
        "RuntimeConfig": "Runtime\nconfig",
        "SolveState": "Solve\nstate",
        "Observer products": "Observer\nproducts",
        "Electron transport": "Electron\ntransport",
        "Photon fields": "Photon\nfields",
        "Hadronic transport": "Hadronic\ntransport",
        "Reverse shock": "Reverse\nshock",
        "EATS projection": "EATS\nprojection",
    }
    write_rows(DATA_DIR / "fig1_architecture_nodes.csv", [
        {"node": name, "x": x, "y": y, "family": color} for name, x, y, color in nodes
    ])
    write_rows(DATA_DIR / "fig1_architecture_edges.csv", [
        {"source": src, "target": dst} for src, dst in edges
    ])

    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.0), gridspec_kw={"width_ratios": [1.35, 1.0]})
    ax = axes[0]
    ax.set_axis_off()
    ax.set_xlim(0.02, 1.0)
    ax.set_ylim(0.16, 0.93)
    for name, x, y, color in nodes:
        draw_box(ax, (x, y), node_labels.get(name, name), color, width=0.16, height=0.105, fontsize=6.1)
    for src, dst in edges:
        sx, sy = node_pos[src]
        tx, ty = node_pos[dst]
        direction = np.sign(tx - sx)
        draw_arrow(ax, (sx + 0.085 * direction, sy), (tx - 0.085 * direction, ty))
    ax.text(0.03, 0.94, "a", fontweight="bold", fontsize=9)
    ax.text(0.08, 0.18, "Fortran kernels", color=TEAL, fontsize=7)
    ax.text(0.49, 0.18, "Python orchestration", color=BLUE, fontsize=7)

    ax = axes[1]
    add_panel(ax, "b")
    stages = ["R", "Gamma", "Ne", "Np", "ngamma", "Fnu"]
    values = np.array([1.0, 0.92, 0.78, 0.56, 0.63, 0.72])
    ax.plot(stages, values, color=BLUE, lw=1.8, marker="o")
    ax.set_ylim(0.45, 1.05)
    ax.set_ylabel("state passed downstream")
    ax.set_title("Radius-first state chain")
    ax.grid(axis="y", color=LIGHT, lw=0.6)
    ax.text(0.02, 0.08, "Observer time is a projection coordinate.", transform=ax.transAxes, fontsize=7)
    save_pub(fig, "fig1_architecture")


def fig2_forward_api() -> None:
    times = np.logspace(2.0, 7.0, 42)
    freqs = np.array([1.0e9, 1.0e14, 1.0e18])
    spec_freq = np.logspace(7.0, 21.0, 80)
    model = top_hat_model(
        fwd_rad=radiation(include_ssc=True, include_kn_correction=False),
        numerics=numerics(
            num_radius=72,
            num_theta=36,
            num_observer_time=72,
            num_electron_gamma=81,
            num_photon_frequency=64,
            num_threads=1,
        ),
    )
    flux = model.flux_density_grid(times, freqs)
    spec_flux = model.flux_density_grid(np.array([1.0e5]), spec_freq)
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
    for nu, total, sync, ssc in zip(
        spec_freq,
        spec_flux.total[:, 0],
        spec_flux.fwd.sync[:, 0],
        spec_flux.fwd.ssc[:, 0],
    ):
        rows.append(
            {
                "kind": "spectrum_t1e5s",
                "time_s": "1.00000000e+05",
                "frequency_hz": f"{nu:.8e}",
                "total_fnu_cgs": f"{total:.8e}",
                "fs_synch_fnu_cgs": f"{sync:.8e}",
                "fs_ssc_fnu_cgs": f"{ssc:.8e}",
            }
        )
    write_rows(DATA_DIR / "fig2_forward_api.csv", rows)

    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.45), gridspec_kw={"width_ratios": [1.1, 1.18, 0.95]})
    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.24, top=0.84, wspace=0.58)
    colors = [BLUE, TEAL, RED]
    labels = ["1 GHz", "Optical", "X-ray"]
    for i, (nu, color, label) in enumerate(zip(freqs, colors, labels)):
        axes[0].loglog(times, flux.total[i], color=color, lw=1.7, label=label)
    add_panel(axes[0], "a")
    axes[0].set_xlabel(r"$t_{\rm obs}$ (s)")
    axes[0].set_ylabel(r"$F_\nu$ (cgs)")
    axes[0].set_title("Public light-curve query")
    axes[0].legend(fontsize=6)
    axes[0].grid(color=LIGHT, which="both", lw=0.5)

    axes[1].loglog(spec_freq, spec_flux.total[:, 0], color=NEUTRAL, lw=1.8, label="total")
    axes[1].loglog(spec_freq, spec_flux.fwd.sync[:, 0], color=BLUE, lw=1.5, label="synch")
    axes[1].loglog(spec_freq, spec_flux.fwd.ssc[:, 0], color=TEAL, lw=1.5, label="SSC")
    add_panel(axes[1], "b")
    axes[1].set_xlabel(r"$\nu$ (Hz)")
    axes[1].set_ylabel(r"$F_\nu$ (cgs)")
    axes[1].set_title(r"Spectrum at $10^5$ s")
    axes[1].legend(fontsize=6)
    axes[1].grid(color=LIGHT, which="both", lw=0.5)

    ratios = np.max(flux.fwd.ssc / flux.total, axis=1)
    axes[2].bar(labels, ratios, color=[BLUE2, TEAL, RED], edgecolor="black", linewidth=0.7)
    add_panel(axes[2], "c")
    axes[2].set_yscale("log")
    axes[2].set_ylim(float(np.min(ratios)) / 2.0, float(np.max(ratios)) * 2.0)
    axes[2].set_ylabel("peak SSC / total", labelpad=2)
    axes[2].set_title("Component accounting")
    axes[2].grid(axis="y", color=LIGHT, lw=0.5)
    save_pub(fig, "fig2_forward_api")


def fig3_transport_projection() -> None:
    runtime_path = ROOT / "output/asgard_doc/runtime_benchmark/runtime_breakdown_summary.csv"
    rows = read_csv(runtime_path)
    data = [
        r
        for r in rows
        if r["solver"] == "fullhide_1d" and r["threads"] in {"1", "8"}
        or r["solver"] == "fullhide_2d" and r["threads"] in {"1", "8"}
    ]
    write_rows(DATA_DIR / "fig3_transport_projection_runtime.csv", data)

    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.45), gridspec_kw={"width_ratios": [1.15, 1.05, 1.0]})
    fig.subplots_adjust(left=0.06, right=0.985, bottom=0.30, top=0.82, wspace=0.65)
    ax = axes[0]
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    draw_box(ax, (0.20, 0.75), "1D\nshells", BLUE, width=0.19, height=0.12, fontsize=6.8)
    draw_box(ax, (0.56, 0.75), r"$N_e$" "\n" r"$(\gamma,R)$", TEAL, width=0.20, height=0.12, fontsize=6.8)
    draw_box(ax, (0.86, 0.75), "thin\nEATS", VIOLET, width=0.19, height=0.12, fontsize=6.8)
    draw_box(ax, (0.20, 0.35), "2D\nshell", BLUE, width=0.19, height=0.12, fontsize=6.8)
    draw_box(ax, (0.56, 0.35), r"$N_e$" "\n" r"$(\gamma,\chi,R)$", TEAL, width=0.20, height=0.12, fontsize=6.8)
    draw_box(ax, (0.86, 0.35), r"$\chi$" "\nEATS", VIOLET, width=0.19, height=0.12, fontsize=6.8)
    for y in (0.75, 0.35):
        draw_arrow(ax, (0.31, y), (0.45, y))
        draw_arrow(ax, (0.67, y), (0.76, y))
    ax.text(0.01, 0.95, "a", fontweight="bold", fontsize=9)
    ax.set_title("Transport-projection contract")

    labels = []
    vals = []
    colors = []
    for r in data:
        labels.append(f"{r['medium']} {r['dimension']} {r['threads']}t")
        vals.append(float(r["high_level_total_median_s"]))
        colors.append(BLUE if r["dimension"] == "1d" else TEAL)
    x = np.arange(len(vals))
    axes[1].bar(x, vals, color=colors, edgecolor="black", linewidth=0.6)
    add_panel(axes[1], "b")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=70, ha="right", fontsize=5.5)
    axes[1].set_ylabel("median wall time (s)")
    axes[1].set_title("Runtime benchmark")
    axes[1].grid(axis="y", color=LIGHT, lw=0.5)

    stages = ["dynamics", "electron", "SSC", r"$\gamma\gamma$", "EATS"]
    row = next(r for r in data if r["medium"] == "ism" and r["dimension"] == "2d" and r["threads"] == "1")
    stage_vals = [
        float(row["dynamics_median_s"]),
        float(row["electron_synch_median_s"]),
        float(row["ssc_median_s"]),
        float(row["annihilation_median_s"]),
        float(row["sed_eats_interpolation_median_s"]),
    ]
    axes[2].bar(stages, stage_vals, color=[NEUTRAL, TEAL, BLUE2, RED, VIOLET], edgecolor="black", linewidth=0.6)
    add_panel(axes[2], "c")
    axes[2].set_yscale("log")
    axes[2].set_ylabel("median stage time (s)")
    axes[2].set_title("2D solve profile")
    axes[2].tick_params(axis="x", rotation=45)
    axes[2].grid(axis="y", color=LIGHT, lw=0.5, which="both")
    save_pub(fig, "fig3_transport_projection")


def fig4_hadronic_feedback() -> None:
    process_rows = [
        {"process": "proton synchrotron", "feeds_observer": 1, "feeds_electrons": 0, "feeds_photons": 0},
        {"process": "p-gamma", "feeds_observer": 1, "feeds_electrons": 0, "feeds_photons": 1},
        {"process": "Bethe-Heitler", "feeds_observer": 1, "feeds_electrons": 1, "feeds_photons": 1},
        {"process": "pp", "feeds_observer": 1, "feeds_electrons": 1, "feeds_photons": 0},
        {"process": "gamma-gamma pair/synch", "feeds_observer": 1, "feeds_electrons": 1, "feeds_photons": 1},
        {"process": "neutrino", "feeds_observer": 1, "feeds_electrons": 0, "feeds_photons": 0},
    ]
    write_rows(DATA_DIR / "fig4_hadronic_feedback_contract.csv", process_rows)
    matrix = np.array([[int(r[k]) for k in ("feeds_observer", "feeds_electrons", "feeds_photons")] for r in process_rows])

    fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.55), gridspec_kw={"width_ratios": [1.08, 1.0]})
    fig.subplots_adjust(left=0.06, right=0.98, bottom=0.28, top=0.82, wspace=0.40)
    ax = axes[0]
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    draw_box(ax, (0.15, 0.70), "shock\nenergy", BLUE, width=0.20, height=0.14, fontsize=6.8)
    draw_box(ax, (0.43, 0.70), "protons", TEAL, width=0.19, height=0.12, fontsize=6.8)
    draw_box(ax, (0.70, 0.70), "photons", GOLD, width=0.19, height=0.12, fontsize=6.8)
    draw_box(ax, (0.58, 0.32), r"secondary" "\n" r"$e^\pm$", RED, width=0.22, height=0.14, fontsize=6.8)
    draw_box(ax, (0.86, 0.50), "observer", VIOLET, width=0.20, height=0.12, fontsize=6.8)
    for start, end in [
        ((0.25, 0.70), (0.34, 0.70)),
        ((0.525, 0.70), (0.605, 0.70)),
        ((0.43, 0.64), (0.53, 0.39)),
        ((0.70, 0.64), (0.63, 0.39)),
        ((0.74, 0.64), (0.82, 0.54)),
        ((0.68, 0.36), (0.78, 0.47)),
    ]:
        draw_arrow(ax, start, end)
    ax.text(0.01, 0.95, "a", fontweight="bold", fontsize=9)
    ax.set_title("Joint shell-level feedback")

    ax = axes[1]
    im = ax.imshow(matrix, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    add_panel(ax, "b")
    ax.set_yticks(np.arange(len(process_rows)))
    ax.set_yticklabels([r["process"] for r in process_rows], fontsize=6)
    ax.set_xticks(np.arange(3))
    ax.set_xticklabels(["observer", r"$Q_{e^\pm}$", "photon sink"], rotation=35, ha="right")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(j, i, "yes" if matrix[i, j] else "-", ha="center", va="center", fontsize=6)
    ax.set_title("Implemented feedback roles")
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

    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.35), gridspec_kw={"width_ratios": [1.02, 1.0, 1.05]})
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
    axes[2].text(0.55, 0.07, rf"$E_e/E_{{\rm diss}}={inj/diss:.2f}$", transform=axes[2].transAxes, ha="center")
    add_panel(axes[2], "c")
    axes[2].grid(axis="x", color=LIGHT, lw=0.5, which="both")
    save_pub(fig, "fig5_reverse_shock")


def fig6_runtime_reuse() -> None:
    rows = read_csv(ROOT / "output/asgard_doc/runtime_benchmark/runtime_breakdown_summary.csv")
    selected = [r for r in rows if r["threads"] == "1"]
    write_rows(DATA_DIR / "fig6_runtime_summary.csv", selected)
    metrics = [
        {"metric": "tracked source files", "value": 179},
        {"metric": "tracked source lines", "value": 44260},
        {"metric": "Fortran kernel files", "value": 87},
        {"metric": "Fortran kernel lines", "value": 21386},
    ]
    write_rows(DATA_DIR / "fig6_code_metrics.csv", metrics)

    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.45), gridspec_kw={"width_ratios": [0.9, 1.0, 1.1]})
    fig.subplots_adjust(left=0.08, right=0.985, bottom=0.29, top=0.83, wspace=0.55)
    labels = [f"{r['medium']} {r['dimension']}" for r in selected]
    total = np.array([float(r["high_level_total_median_s"]) for r in selected])
    x = np.arange(len(total))
    axes[0].bar(x, total, color=[BLUE if "1d" in label else TEAL for label in labels], edgecolor="black", linewidth=0.6)
    add_panel(axes[0], "a")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels, rotation=45, ha="right")
    axes[0].set_ylabel("median wall time (s)")
    axes[0].set_title("End-to-end query cost")
    axes[0].grid(axis="y", color=LIGHT, lw=0.5)

    stage_names = ["dynamics", "electron", "SSC", r"$\gamma\gamma$", "EATS", "Python"]
    r = next(row for row in selected if row["medium"] == "ism" and row["dimension"] == "1d")
    vals = np.array([
        float(r["dynamics_median_s"]),
        float(r["electron_synch_median_s"]),
        float(r["ssc_median_s"]),
        float(r["annihilation_median_s"]),
        float(r["sed_eats_interpolation_median_s"]),
        float(r["python_middle_layer_median_s"]),
    ])
    y_stage = np.arange(len(vals))
    axes[1].barh(y_stage, vals, color=[NEUTRAL, TEAL, BLUE2, RED, VIOLET, GOLD], edgecolor="black", linewidth=0.6)
    axes[1].set_yticks(y_stage)
    axes[1].set_yticklabels(stage_names, fontsize=6.2)
    axes[1].set_xscale("log")
    axes[1].set_xlabel("median stage time (s)")
    axes[1].grid(axis="x", color=LIGHT, lw=0.5, which="both")
    add_panel(axes[1], "b")
    axes[1].set_title("1D ISM stage share")

    metric_labels = ["source files", "source lines", "Fortran files", "Fortran lines"]
    axes[2].barh(metric_labels, [m["value"] for m in metrics], color=[NEUTRAL, NEUTRAL, VIOLET, VIOLET])
    add_panel(axes[2], "c")
    axes[2].set_xscale("log")
    axes[2].set_xlabel("count")
    axes[2].set_title("Codebase scale")
    axes[2].grid(axis="x", color=LIGHT, lw=0.5, which="both")
    save_pub(fig, "fig6_runtime_reuse")


def main() -> None:
    ensure_dirs()
    fig1_architecture()
    fig2_forward_api()
    fig3_transport_projection()
    fig4_hadronic_feedback()
    fig5_reverse_shock()
    fig6_runtime_reuse()


if __name__ == "__main__":
    main()
