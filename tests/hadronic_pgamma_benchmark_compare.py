from __future__ import annotations

from pathlib import Path
import sys

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

import hadronic_pgamma_lepton_reference as lep_ref
from hadronic_pgamma_benchmark_report import (
    OUTPUT_DIR,
    assert_benchmark_summary,
    build_error_summary,
    collect_benchmark_payload,
    write_report,
)


def _save_pair(fig: plt.Figure, stem: Path) -> tuple[Path, Path]:
    png = stem.with_suffix(".png")
    pdf = stem.with_suffix(".pdf")
    fig.savefig(png, dpi=220, bbox_inches="tight")
    fig.savefig(pdf, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return png, pdf


def _plot_secondary_phi_regression(payload: dict[str, object], output_dir: Path) -> tuple[Path, Path]:
    lepton_actual = payload["lepton"]["actual"]
    lepton_reference = payload["lepton"]["reference"]
    x_grid = lepton_actual["x_grid"]

    fig, axes = plt.subplots(2, 3, figsize=(14.0, 8.2), constrained_layout=True)
    for idx, channel in enumerate(lep_ref.CHANNEL_ORDER):
        ax = axes.flat[idx]
        phi_ref = lepton_reference["phi"][channel]
        phi_act = lepton_actual["phi"][channel]
        rel = np.abs(phi_act - phi_ref) / np.maximum(np.abs(phi_ref), 1.0e-300)

        ax.loglog(x_grid, np.maximum(x_grid * phi_ref, 1.0e-300), "--", linewidth=1.6, label="reference")
        ax.loglog(x_grid, np.maximum(x_grid * phi_act, 1.0e-300), "-", linewidth=1.3, label="current")
        panel_ymax = float(np.max(np.maximum([x_grid * phi_ref, x_grid * phi_act], 1.0e-300)))
        ax.set_ylim(bottom=panel_ymax * 1.0e-15)
        ax.set_title(f"{channel} (max rel={np.max(rel):.1e})", fontsize=10)
        ax.set_xlabel("x = Esec / Ep")
        ax.set_ylabel("x Phi")
        ax.grid(True, which="both", alpha=0.25, linestyle=":")
        ax.legend(fontsize=8)

    return _save_pair(fig, output_dir / "ka2008_secondary_phi_regression")


def _plot_loss_gamma_regression(payload: dict[str, object], output_dir: Path) -> tuple[Path, Path]:
    loss_actual = payload["loss"]["actual"]
    loss_reference = payload["loss"]["reference"]
    lepton_actual = payload["lepton"]["actual"]
    lepton_reference = payload["lepton"]["reference"]

    gamma_grid = loss_actual["gamma_grid"]
    loss_ref = loss_reference["loss_rate"]
    loss_act = loss_actual["loss_rate"]
    loss_rel = np.abs(loss_act - loss_ref) / np.maximum(np.abs(loss_ref), 1.0e-300)

    e_grid = lepton_actual["spectrum_energy_grid"]
    gamma_spec_ref = lepton_reference["gamma_spectrum"]
    gamma_spec_act = lepton_actual["gamma_spectrum"]
    gamma_spec_rel = np.abs(gamma_spec_act - gamma_spec_ref) / np.maximum(np.abs(gamma_spec_ref), 1.0e-300)

    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.0), constrained_layout=True)
    axes[0, 0].loglog(gamma_grid, loss_ref, "--", label="reference")
    axes[0, 0].loglog(gamma_grid, loss_act, "-", label="current")
    axes[0, 0].set_ylim(bottom=float(np.max(np.maximum([loss_ref, loss_act], 1.0e-300))) * 1.0e-15)
    axes[0, 0].set_title("strict isotropic p-gamma loss rate")
    axes[0, 0].set_xlabel("gamma_p")
    axes[0, 0].set_ylabel("t^{-1}_pgamma")
    axes[0, 0].grid(True, which="both", alpha=0.25, linestyle=":")
    axes[0, 0].legend(fontsize=8)

    axes[1, 0].semilogx(gamma_grid, loss_rel, color="tab:red")
    axes[1, 0].set_xlabel("gamma_p")
    axes[1, 0].set_ylabel("relative error")
    axes[1, 0].grid(True, which="both", alpha=0.25, linestyle=":")
    axes[1, 0].set_title("loss relative error")

    axes[0, 1].loglog(e_grid, gamma_spec_ref, "--", label="reference")
    axes[0, 1].loglog(e_grid, gamma_spec_act, "-", label="current")
    axes[0, 1].set_ylim(bottom=float(np.max(np.maximum([gamma_spec_ref, gamma_spec_act], 1.0e-300))) * 1.0e-15)
    axes[0, 1].set_title("KA2008 gamma spectrum")
    axes[0, 1].set_xlabel("E_gamma [GeV]")
    axes[0, 1].set_ylabel("Q_gamma")
    axes[0, 1].grid(True, which="both", alpha=0.25, linestyle=":")
    axes[0, 1].legend(fontsize=8)

    axes[1, 1].semilogx(e_grid, gamma_spec_rel, color="tab:red")
    axes[1, 1].set_xlabel("E_gamma [GeV]")
    axes[1, 1].set_ylabel("relative error")
    axes[1, 1].grid(True, which="both", alpha=0.25, linestyle=":")
    axes[1, 1].set_title("gamma spectrum relative error")

    return _save_pair(fig, output_dir / "ka2008_loss_gamma_regression")


def _plot_error_summary(summary: dict[str, object], output_dir: Path) -> tuple[Path, Path]:
    labels = [
        "kinematics",
        "loss",
        "lepton phi",
        "lepton spectrum",
        "gamma spectrum",
    ]
    values = [
        summary["kinematics"]["max_rel_error"],
        summary["loss"]["max_rel_error"],
        summary["lepton"]["phi_max_rel_error"],
        summary["lepton"]["spectrum_max_rel_error"],
        summary["lepton"]["gamma_spectrum_max_rel_error"],
    ]

    fig, ax = plt.subplots(figsize=(8.8, 4.4), constrained_layout=True)
    ax.bar(labels, values, color=["#4e79a7", "#f28e2b", "#59a14f", "#e15759", "#76b7b2"])
    ax.set_yscale("log")
    ax.set_ylabel("max relative error")
    ax.set_title("ka2008_reference / hummer_2010_response regression error summary")
    ax.grid(True, which="both", axis="y", alpha=0.25, linestyle=":")
    return _save_pair(fig, output_dir / "ka2008_hummer_error_summary")


def main() -> None:
    output_dir = OUTPUT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    payload = collect_benchmark_payload()
    summary = build_error_summary(payload)
    assert_benchmark_summary(payload, summary)

    figure_paths: list[Path] = []
    figure_paths.extend(_plot_secondary_phi_regression(payload, output_dir))
    figure_paths.extend(_plot_loss_gamma_regression(payload, output_dir))
    figure_paths.extend(_plot_error_summary(summary, output_dir))
    summary_json_path, summary_md_path = write_report(summary, output_dir=output_dir, figure_paths=figure_paths)

    print(f"global_max_rel_error={summary['global_max_rel_error']:.3e}")
    print(f"saved={summary_json_path}")
    print(f"saved={summary_md_path}")
    for path in figure_paths:
        print(f"saved={path}")


if __name__ == "__main__":
    main()
