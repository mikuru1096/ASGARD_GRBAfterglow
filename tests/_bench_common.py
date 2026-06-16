from __future__ import annotations

from _repo_path import ensure_repo_root_on_path
import time
from typing import Callable


ensure_repo_root_on_path()


def run_case(name: str, fn: Callable[[], object]) -> dict[str, object]:
    print(f"  {name} ...", flush=True)
    t0 = time.perf_counter()
    try:
        payload = fn()
        item = {"name": name, "status": "PASS", "seconds": time.perf_counter() - t0, "payload": payload}
    except NotImplementedError as exc:
        item = {"name": name, "status": "UNSUPPORTED", "seconds": time.perf_counter() - t0, "payload": str(exc)}
    except Exception as exc:
        item = {"name": name, "status": "FAIL", "seconds": time.perf_counter() - t0, "payload": f"{type(exc).__name__}: {exc}"}
    print(f"  {name}: {item['status']} {item['seconds']:.2f}s", flush=True)
    return item
