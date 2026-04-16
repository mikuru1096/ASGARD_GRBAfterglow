from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUTPUT_ROOT = ROOT / "output"
ASGARD_DOC_DIR = OUTPUT_ROOT / "asgard_doc"
BENCHMARK_EXP_TAIL_DIR = OUTPUT_ROOT / "benchmark_exp_tail"


@dataclass(frozen=True)
class OutputPaths:
    root: Path = ROOT
    output_root: Path = OUTPUT_ROOT
    asgard_doc_dir: Path = ASGARD_DOC_DIR
    benchmark_exp_tail_dir: Path = BENCHMARK_EXP_TAIL_DIR


PATHS = OutputPaths()


def ensure_output_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def asgard_doc_path(name: str) -> Path:
    return ASGARD_DOC_DIR / name


def benchmark_exp_tail_path(name: str) -> Path:
    return BENCHMARK_EXP_TAIL_DIR / name
