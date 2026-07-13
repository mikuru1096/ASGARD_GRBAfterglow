"""Shared, functional contracts for reproducible ASGARD benchmarks."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import platform
import resource
import subprocess
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
DATA_ROOT = ROOT / "paper" / "source_data" / "benchmarks"
FIGURE_ROOT = ROOT / "paper" / "figures" / "benchmarks"
PALETTE = {
    "black": "#000000",
    "orange": "#E69F00",
    "sky": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
}


def plot_style() -> dict[str, object]:
    """Return the common print-readable matplotlib style."""
    return {
        "font.family": "sans-serif",
        "font.size": 8.0,
        "axes.labelsize": 8.0,
        "axes.titlesize": 8.0,
        "xtick.labelsize": 7.0,
        "ytick.labelsize": 7.0,
        "legend.fontsize": 7.0,
        "lines.linewidth": 1.25,
        "lines.markersize": 4.0,
        "axes.linewidth": 0.7,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "savefig.bbox": "tight",
    }


def panel_label(axis, label: str) -> None:
    axis.text(0.02, 0.98, label, transform=axis.transAxes, ha="left", va="top", fontweight="bold")


def significant(values, fraction: float) -> np.ndarray:
    """Select finite values above a stated fraction of their physical peak."""
    data = np.asarray(values)
    return np.isfinite(data) & (data > fraction * np.nanmax(data))


def joint_significant(left, right, fraction: float) -> np.ndarray:
    return significant(left, fraction) & significant(right, fraction)


def summary(values) -> dict[str, float | int]:
    """Report robust timing/error statistics without changing the samples."""
    data = np.asarray(values, dtype=float)
    data = data[np.isfinite(data)]
    q25, median, q75, p95 = np.percentile(data, (25, 50, 75, 95))
    return {
        "count": int(data.size),
        "median": float(median),
        "iqr": float(q75 - q25),
        "p50": float(median),
        "p95": float(p95),
    }


def peak_rss() -> float:
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return float(rss / (1024 if sys.platform != "darwin" else 1024**2))


def save_figure(figure, path: Path, *, dpi: int = 600) -> list[Path]:
    """Export an editable PDF master and a 600-dpi PNG."""
    path.parent.mkdir(parents=True, exist_ok=True)
    outputs = [path.with_suffix(".pdf"), path.with_suffix(".png")]
    for output in outputs:
        figure.savefig(output, dpi=dpi if output.suffix == ".png" else None)
    return outputs


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _command(*args: str) -> str:
    return subprocess.run(
        args, cwd=ROOT, check=True, text=True, encoding="utf-8", capture_output=True
    ).stdout.strip()


def environment(mode: str, *, threads: int, grid: object, repeats: int) -> dict[str, object]:
    packages = {}
    for name in ("numpy", "scipy", "matplotlib", "afterglowpy", "jetsimpy", "vegasafterglow", "pyblastafterglowmag"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    return {
        "mode": mode,
        "git_commit": _command("git", "rev-parse", "HEAD"),
        "git_dirty": bool(_command("git", "status", "--porcelain")),
        "python": sys.version.split()[0],
        "compiler": _command("gfortran", "--version").splitlines()[0],
        "platform": platform.platform(),
        "os": platform.system(),
        "cpu": platform.processor() or _command("uname", "-m"),
        "threads": threads,
        "grid": grid,
        "repeats": repeats,
        "packages": packages,
    }


def write_json(path: Path, payload: object) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def manifest(data_root: Path = DATA_ROOT, figure_root: Path = FIGURE_ROOT) -> dict[str, object]:
    target = data_root / "manifest.json"
    files = sorted(
        path
        for root in (data_root, figure_root)
        if root.exists()
        for path in root.rglob("*")
        if path.is_file() and path != target
    )
    return {
        "files": [
            {"path": path.relative_to(ROOT).as_posix(), "bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in files
        ]
    }


def write_manifest(path: Path = DATA_ROOT / "manifest.json") -> Path:
    return write_json(path, manifest())
