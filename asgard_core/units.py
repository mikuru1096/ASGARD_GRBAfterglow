from __future__ import annotations

from math import pi

rad = 1.0
deg = pi / 180.0
arcmin = deg / 60.0
arcsec = arcmin / 60.0
mas = 1.0e-3 * arcsec
uas = 1.0e-6 * arcsec

second = 1.0
minute = 60.0 * second
hour = 60.0 * minute
day = 24.0 * hour
year = 365.25 * day

Hz = 1.0
kHz = 1.0e3 * Hz
MHz = 1.0e6 * Hz
GHz = 1.0e9 * Hz

Jy = 1.0
mJy = 1.0e-3 * Jy
uJy = 1.0e-6 * Jy
nJy = 1.0e-9 * Jy

eV = 2.417989242e14 * Hz
keV = 1.0e3 * eV
MeV = 1.0e6 * eV
GeV = 1.0e9 * eV
TeV = 1.0e12 * eV
