import os

if os.name == "nt" and hasattr(os, "add_dll_directory"):
    _mingw_bin = os.environ.get("ASGARD_MINGW_BIN", r"C:\msys64\mingw64\bin")
    if os.path.isdir(_mingw_bin):
        os.add_dll_directory(_mingw_bin)

from src.Constants import constants

__all__ = [
    "constants",
    "Dynamics",
    "Electron",
    "Hadronic",
    "Interpolation",
    "Radiation",
]

__doc__ = "ASGARD Fortran runtime package."
