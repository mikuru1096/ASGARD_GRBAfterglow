try:
    from src.Electron.FS_electron_fullhide_v3 import fs_electron_fullhide
except ImportError:
    try:
        from src.Electron.FS_electron_fullhide_v2 import fs_electron_fullhide
    except ImportError:
        from src.Electron.FS_electron_fullhide import fs_electron_fullhide

try:
    from src.Electron.FS_electron_slc1 import fs_electron_slc1
except ImportError:
    try:
        from src.Electron.FS_electron_slc1_v3 import fs_electron_slc1
    except ImportError:
        from src.Electron.FS_electron_slc1_v2 import fs_electron_slc1

from src.Electron.FS_electron_charint import fs_electron_charint

from src.Electron.FS_electron_weno5 import fs_electron_weno5
try:
    from src.Electron.FS_electron_t2g1_v3 import fs_electron_t2g1
except ImportError:
    try:
        from src.Electron.FS_electron_t2g1_v2 import fs_electron_t2g1
    except ImportError:
        from src.Electron.FS_electron_t2g1 import fs_electron_t2g1
#from src.Electron.FS_electron_t2g2 import fs_electron_t2g2

__all__ = [
           "fs_electron_fullhide",
           "fs_electron_slc1",
           "fs_electron_charint",
           "fs_electron_weno5",
           "fs_electron_t2g1",
#           "fs_electron_t2g2"
          ]
