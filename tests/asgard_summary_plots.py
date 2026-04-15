from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap


ROOT = Path(__file__).resolve().parents[1]
INPUT_JSON = ROOT / "output" / "asgard_doc" / "comprehensive_validation_asgard.json"
OUTPUT_DIR = ROOT / "output" / "asgard_doc"
OUTPUT_PDF = OUTPUT_DIR / "comprehensive_validation_summary.pdf"


GOOD_BAD = LinearSegmentedColormap.from_list("good_bad", ["#2e7d32", "#f9a825", "#c62828"])


def _load_payload() -> dict:
    return json.loads(INPUT_JSON.read_text(encoding="utf-8"))


def _ratio(item: dict) -> float:
    measured = item.get("measured")
    expected = item.get("expected")
    tolerance = item.get("tolerance")
    if measured is None or expected is None or tolerance in (None, 0):
        return np.nan
    return abs(float(measured) - float(expected)) / float(tolerance)


def _split_name(item: dict) -> tuple[str, str]:
    name = item["name"]
    if name.startswith("fwd-"):
        _, medium, phase, qty = name.split("-", 3)
        return f"FS {medium} {phase}", qty
    if name.startswith("rvs-"):
        _, shell, medium, phase, qty = name.split("-", 4)
        return f"RS {shell} {medium} {phase}", qty
    if name.startswith("regime-"):
        _, regime, segment = name.split("-", 2)
        return f"Regime {regime}", segment
    return name, "value"


def _collect_matrix(items: list[dict], qty_order: list[str]) -> tuple[list[str], np.ndarray, dict[tuple[int, int], dict]]:
    row_labels = []
    for item in items:
        label, _ = _split_name(item)
        if label not in row_labels:
            row_labels.append(label)
    rows = {label: i for i, label in enumerate(row_labels)}
    cols = {qty: j for j, qty in enumerate(qty_order)}
    matrix = np.full((len(row_labels), len(qty_order)), np.nan, dtype=float)
    meta: dict[tuple[int, int], dict] = {}
    for item in items:
        label, qty = _split_name(item)
        if qty not in cols:
            continue
        i = rows[label]
        j = cols[qty]
        matrix[i, j] = _ratio(item)
        meta[(i, j)] = item
    return row_labels, matrix, meta


def _draw_heatmap(ax, title: str, items: list[dict], qty_order: list[str]) -> None:
    row_labels, matrix, meta = _collect_matrix(items, qty_order)
    if matrix.size == 0:
        ax.axis("off")
        ax.set_title(title)
        return
    masked = np.ma.masked_invalid(matrix)
    im = ax.imshow(masked, aspect="auto", cmap=GOOD_BAD, vmin=0.0, vmax=3.0)
    ax.set_title(title, fontsize=12, pad=10)
    ax.set_xticks(range(len(qty_order)))
    ax.set_xticklabels(qty_order, rotation=0)
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels)
    ax.set_xlabel("|measured - expected| / tolerance")
    for (i, j), item in meta.items():
        measured = item.get("measured")
        expected = item.get("expected")
        passed = item.get("passed", False)
        if measured is None or expected is None:
            text = "n/a"
        else:
            text = f"{measured:.2f}\n[{expected:.2f}]"
        ax.text(j, i, text, ha="center", va="center", fontsize=8, color="white" if not np.isnan(matrix[i, j]) and matrix[i, j] > 1.0 else "black")
        if not passed:
            ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor="black", linewidth=1.2))
    for x in np.arange(-0.5, len(qty_order), 1.0):
        ax.axvline(x, color="white", linewidth=0.6, alpha=0.7)
    for y in np.arange(-0.5, len(row_labels), 1.0):
        ax.axhline(y, color="white", linewidth=0.6, alpha=0.7)
    return im


def _category_counts(regression: list[dict]) -> dict[str, tuple[int, int]]:
    groups = {
        "FS dynamics": [],
        "RS dynamics": [],
        "FS/RS frequencies": [],
        "Spectral regimes": [],
    }
    for item in regression:
        _, qty = _split_name(item)
        if item["name"].startswith("fwd-") and qty in {"u", "r", "B", "N_p"}:
            groups["FS dynamics"].append(item)
        elif item["name"].startswith("rvs-") and qty in {"u", "r", "B", "N_p"}:
            groups["RS dynamics"].append(item)
        elif qty.startswith("nu_"):
            groups["FS/RS frequencies"].append(item)
        elif item["name"].startswith("regime-"):
            groups["Spectral regimes"].append(item)
    counts = {}
    for key, values in groups.items():
        passed = sum(1 for item in values if item["passed"])
        counts[key] = (passed, len(values))
    return counts


def _draw_overview(fig, payload: dict) -> None:
    regression = payload["regression"]
    counts = _category_counts(regression)
    timings = payload["benchmark"]

    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.3], width_ratios=[1.0, 1.1])
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[1, :])

    ax0.axis("off")
    summary = payload["summary"]
    text = (
        "ASGARD vs ASGARD Summary\n\n"
        f"Benchmark pass: {summary['benchmark_pass']}\n"
        f"Regression pass: {summary['regression_pass']}\n"
        f"Regression total: {summary['regression_total']}\n"
        f"Regression failed: {summary['regression_failed']}\n"
        f"Generated from: {INPUT_JSON.name}"
    )
    ax0.text(0.02, 0.98, text, va="top", ha="left", fontsize=12, family="monospace")

    labels = list(counts.keys())
    ratios = [passed / total if total else 0.0 for passed, total in counts.values()]
    ax1.barh(labels, ratios, color=["#2e7d32" if x >= 0.8 else "#f9a825" if x >= 0.5 else "#c62828" for x in ratios])
    ax1.set_xlim(0.0, 1.0)
    ax1.set_xlabel("Pass fraction")
    ax1.set_title("Category pass rates")
    for i, (label, (passed, total)) in enumerate(counts.items()):
        ax1.text(ratios[i] + 0.02, i, f"{passed}/{total}", va="center", fontsize=10)

    bench_labels = [f"{b['config']['jet']}/{b['config']['medium']}/{b['config']['radiation']}" for b in timings]
    bench_times = [b["timing_ms"] for b in timings]
    ax2.barh(bench_labels, bench_times, color="#1565c0")
    ax2.set_xlabel("Timing [ms]")
    ax2.set_title("Benchmark overview")
    for i, ms in enumerate(bench_times):
        ax2.text(ms * 1.01, i, f"{ms:.1f}", va="center", fontsize=9)

    fig.suptitle("Comprehensive Validation Overview", fontsize=16)


def _draw_dynamics(fig, payload: dict) -> None:
    regression = payload["regression"]
    fwd_dyn = [item for item in regression if item["name"].startswith("fwd-") and _split_name(item)[1] in {"u", "r", "B", "N_p"}]
    rvs_dyn = [item for item in regression if item["name"].startswith("rvs-") and _split_name(item)[1] in {"u", "r", "B", "N_p"}]

    gs = fig.add_gridspec(2, 1, hspace=0.35)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[1, 0])
    im = _draw_heatmap(ax0, f"Forward-shock dynamics ({sum(i['passed'] for i in fwd_dyn)}/{len(fwd_dyn)})", fwd_dyn, ["u", "r", "B", "N_p"])
    _draw_heatmap(ax1, f"Reverse-shock dynamics ({sum(i['passed'] for i in rvs_dyn)}/{len(rvs_dyn)})", rvs_dyn, ["r", "B", "N_p"])
    cbar = fig.colorbar(im, ax=[ax0, ax1], shrink=0.85, pad=0.02)
    cbar.set_label("|measured - expected| / tolerance")
    fig.suptitle("Dynamics Summary", fontsize=16)


def _draw_frequencies(fig, payload: dict) -> None:
    regression = payload["regression"]
    fwd_freq = [item for item in regression if item["name"].startswith("fwd-") and _split_name(item)[1].startswith("nu_")]
    rvs_freq = [item for item in regression if item["name"].startswith("rvs-") and _split_name(item)[1].startswith("nu_")]

    gs = fig.add_gridspec(2, 1, hspace=0.35)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[1, 0])
    im = _draw_heatmap(ax0, f"Forward-shock frequencies ({sum(i['passed'] for i in fwd_freq)}/{len(fwd_freq)})", fwd_freq, ["nu_m", "nu_c", "nu_M"])
    _draw_heatmap(ax1, f"Reverse-shock frequencies ({sum(i['passed'] for i in rvs_freq)}/{len(rvs_freq)})", rvs_freq, ["nu_m", "nu_c", "nu_M"])
    cbar = fig.colorbar(im, ax=[ax0, ax1], shrink=0.85, pad=0.02)
    cbar.set_label("|measured - expected| / tolerance")
    fig.suptitle("Frequency Evolution Summary", fontsize=16)


def _draw_spectrum(fig, payload: dict) -> None:
    spectrum = [item for item in payload["regression"] if item["name"].startswith("regime-")]
    labels = [item["name"].replace("regime-", "") for item in spectrum]
    measured = np.array([np.nan if item["measured"] is None else item["measured"] for item in spectrum], dtype=float)
    expected = np.array([np.nan if item["expected"] is None else item["expected"] for item in spectrum], dtype=float)
    passed = np.array([item["passed"] for item in spectrum], dtype=bool)

    gs = fig.add_gridspec(1, 2, width_ratios=[1.6, 1.0], wspace=0.35)
    ax0 = fig.add_subplot(gs[0, 0])
    y = np.arange(len(labels))
    ax0.scatter(expected, y, marker="x", s=70, color="black", label="expected")
    ax0.scatter(measured, y, s=45, color=np.where(passed, "#2e7d32", "#c62828"), label="measured")
    ax0.set_yticks(y)
    ax0.set_yticklabels(labels, fontsize=9)
    ax0.set_xlabel("Spectral slope")
    ax0.set_title("Measured vs expected")
    ax0.axvline(0.0, color="0.75", linewidth=0.8)
    ax0.legend(loc="lower right")

    ax1 = fig.add_subplot(gs[0, 1])
    ratios = np.array([_ratio(item) for item in spectrum], dtype=float)
    ax1.barh(y, ratios, color=np.where(passed, "#2e7d32", "#c62828"))
    ax1.axvline(1.0, color="black", linestyle="--", linewidth=1.0)
    ax1.set_yticks(y)
    ax1.set_yticklabels([])
    ax1.set_xlabel("|measured - expected| / tolerance")
    ax1.set_title(f"Spectral regimes ({passed.sum()}/{len(passed)})")

    fig.suptitle("Spectral Regime Summary", fontsize=16)


def _draw_benchmark(fig, payload: dict) -> None:
    benchmark = payload["benchmark"]
    labels = [f"{item['config']['jet']}/{item['config']['medium']}" for item in benchmark]
    timings = np.array([item["timing_ms"] for item in benchmark], dtype=float)

    gs = fig.add_gridspec(2, 1, hspace=0.35)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[1, 0])

    ax0.barh(labels, timings, color="#1565c0")
    ax0.set_xlabel("Timing [ms]")
    ax0.set_title("Benchmark timing")
    for i, ms in enumerate(timings):
        ax0.text(ms * 1.01, i, f"{ms:.1f}", va="center", fontsize=9)

    dim_labels = []
    mean_vals = []
    max_vals = []
    for item in benchmark:
        cfg_label = f"{item['config']['jet']}/{item['config']['medium']}"
        for dim, stats in item["convergence"].items():
            dim_labels.append(f"{cfg_label} {dim}")
            mean_vals.append(stats["worst_mean_error"])
            max_vals.append(stats["worst_max_error"])
    y = np.arange(len(dim_labels))
    if dim_labels:
        ax1.barh(y - 0.18, mean_vals, height=0.35, color="#43a047", label="worst mean")
        ax1.barh(y + 0.18, max_vals, height=0.35, color="#ef6c00", label="worst max")
        ax1.axvline(0.05, color="#43a047", linestyle="--", linewidth=1.0)
        ax1.axvline(0.15, color="#ef6c00", linestyle="--", linewidth=1.0)
        ax1.set_yticks(y)
        ax1.set_yticklabels(dim_labels)
        ax1.legend(loc="lower right")
    else:
        ax1.text(0.5, 0.5, "No convergence scans recorded for structured timing-only entries.", ha="center", va="center", transform=ax1.transAxes)
    ax1.set_xlabel("Relative error")
    ax1.set_title("Convergence summary")

    fig.suptitle("Benchmark Overview", fontsize=16)


def _save_png(fig, name: str) -> None:
    path = OUTPUT_DIR / name
    fig.savefig(path, dpi=200, bbox_inches="tight")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    payload = _load_payload()

    with PdfPages(OUTPUT_PDF) as pdf:
        for draw_fn, png_name in [
            (_draw_overview, "validation-summary-overview.png"),
            (_draw_dynamics, "validation-summary-dynamics.png"),
            (_draw_frequencies, "validation-summary-frequencies.png"),
            (_draw_spectrum, "validation-summary-spectrum.png"),
            (_draw_benchmark, "validation-summary-benchmark.png"),
        ]:
            fig = plt.figure(figsize=(12, 9), constrained_layout=True)
            draw_fn(fig, payload)
            pdf.savefig(fig)
            _save_png(fig, png_name)
            plt.close(fig)

    print(OUTPUT_PDF)


if __name__ == "__main__":
    main()
