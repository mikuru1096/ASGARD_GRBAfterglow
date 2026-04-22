from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


PACKAGE_ROOT = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[1]
OUTPUT_ROOT = ROOT / "output"
ASGARD_DOC_DIR = OUTPUT_ROOT / "asgard_doc"
BENCHMARK_EXP_TAIL_DIR = OUTPUT_ROOT / "benchmark_exp_tail"
DATA_LIGHT_CURVE_DIR = PACKAGE_ROOT / "data_light_curve"
DATA_SPECTRUM_DIR = PACKAGE_ROOT / "data_spectrum"


@dataclass(frozen=True)
class OutputPaths:
    root: Path = ROOT
    package_root: Path = PACKAGE_ROOT
    output_root: Path = OUTPUT_ROOT
    asgard_doc_dir: Path = ASGARD_DOC_DIR
    benchmark_exp_tail_dir: Path = BENCHMARK_EXP_TAIL_DIR
    data_light_curve_dir: Path = DATA_LIGHT_CURVE_DIR
    data_spectrum_dir: Path = DATA_SPECTRUM_DIR


PATHS = OutputPaths()


def ensure_output_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def asgard_doc_path(name: str) -> Path:
    return ASGARD_DOC_DIR / name


def benchmark_exp_tail_path(name: str) -> Path:
    return BENCHMARK_EXP_TAIL_DIR / name


def data_light_curve_path(name: str) -> Path:
    return DATA_LIGHT_CURVE_DIR / name


def data_spectrum_path(name: str) -> Path:
    return DATA_SPECTRUM_DIR / name
