from pathlib import Path
import sys


def ensure_repo_root_on_path() -> None:
    root = str(Path(__file__).resolve().parents[1])
    if root not in sys.path:
        sys.path.insert(0, root)
