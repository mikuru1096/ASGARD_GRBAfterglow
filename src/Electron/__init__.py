from src.Electron.FS_electron_charint_1d import fs_electron_charint_1d
from src.Electron.FS_electron_fullhide_1d import fs_electron_fullhide_1d
from src.Electron.FS_electron_fullhide_2d import fs_electron_fullhide_2d
from src.Electron.FS_electron_slc1_1d import fs_electron_slc1_1d
from src.Electron.FS_electron_t2g1_1d import fs_electron_t2g1_1d
from src.Electron.FS_electron_weno5_1d import fs_electron_weno5_1d

fs_electron_charint = fs_electron_charint_1d
fs_electron_fullhide = fs_electron_fullhide_1d
fs_electron_slc1 = fs_electron_slc1_1d
fs_electron_t2g1 = fs_electron_t2g1_1d
fs_electron_weno5 = fs_electron_weno5_1d

__all__ = [
    "fs_electron_charint",
    "fs_electron_charint_1d",
    "fs_electron_fullhide",
    "fs_electron_fullhide_1d",
    "fs_electron_fullhide_2d",
    "fs_electron_slc1",
    "fs_electron_slc1_1d",
    "fs_electron_t2g1",
    "fs_electron_t2g1_1d",
    "fs_electron_weno5",
    "fs_electron_weno5_1d",
]

__doc__ = "Electron solver bindings."
