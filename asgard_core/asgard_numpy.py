from __future__ import annotations

import numpy as np


def trapezoid(y, x=None, axis: int = -1):
    fn = getattr(np, "trapezoid", None)
    if fn is not None:
        return fn(y, x=x, axis=axis)
    return np.trapz(y, x=x, axis=axis)
