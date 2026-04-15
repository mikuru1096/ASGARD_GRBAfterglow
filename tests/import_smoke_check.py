from __future__ import annotations

from pathlib import Path
import importlib
import os
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def main() -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        prev_cwd = Path.cwd()
        os.chdir(tmpdir)
        try:
            hand_my = importlib.import_module("hand_my")
            hand_reverse = importlib.import_module("hand_reverse")
        finally:
            os.chdir(prev_cwd)

        assert hasattr(hand_my, "build_demo_model")
        assert hasattr(hand_reverse, "build_demo_model")
        assert not any(Path(tmpdir).glob("*.pdf"))

    print("PASS: import smoke check succeeded.")


if __name__ == "__main__":
    main()
