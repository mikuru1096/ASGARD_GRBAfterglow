from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from src import constants

import src.Hadronic.hadronic_forward_1d as hadronic_module


HUMMER_OPERATOR = "fortran_operator"
HUMMER_DECAY = "fortran_decay"
ACCEL_BACKEND = "fortran_acceleration"
BH_BACKEND = "fortran_bethe_heitler"
HIC_BACKEND = "fortran_hadronic_ic"
PP_BACKEND = "fortran_pp_shell"
SECONDARY_BACKEND = "fortran_secondary_radiation"
SPECIES_BACKEND = "fortran_species_transport"

MP_GEV = constants.para_m_p_gev


def strict_grid(values: np.ndarray, name: str, *, finite: bool = False) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be a 1d array.")
    if arr.size < 2:
        raise ValueError(f"{name} must contain at least 2 points.")
    if np.any(arr <= 0.0):
        raise ValueError(f"{name} must be strictly positive.")
    if np.any(np.diff(arr) <= 0.0):
        raise ValueError(f"{name} must be strictly increasing.")
    if finite and not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain finite values.")
    return arr


def same_grid(
    values: np.ndarray,
    reference: np.ndarray,
    name: str,
    *,
    finite: bool = False,
    nonnegative: bool = False,
) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.shape != reference.shape:
        raise ValueError(f"{name} must match grid shape.")
    if finite and not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain finite values.")
    if nonnegative and np.any(arr < 0.0):
        raise ValueError(f"{name} must be non-negative.")
    return arr


def loginterp(x_src: np.ndarray, y_src: np.ndarray, x_dst: np.ndarray) -> np.ndarray:
    x = np.asarray(x_src, dtype=float)
    y = np.asarray(y_src, dtype=float)
    target = np.asarray(x_dst, dtype=float)
    out = np.zeros_like(target, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y) & (x > 0.0) & (y > 0.0)
    if np.count_nonzero(valid) < 2:
        return out
    xv = x[valid]
    yv = y[valid]
    inside = np.isfinite(target) & (target >= xv[0]) & (target <= xv[-1])
    out[inside] = np.exp(np.interp(np.log(target[inside]), np.log(xv), np.log(yv)))
    return out


@dataclass(frozen=True, slots=True)
class HummerOutput:
    """Hummer-2010 pγ wrapper output; rate fields are per GeV."""

    gamma_gev: np.ndarray
    gamma_rate: np.ndarray
    neutrino_gev: np.ndarray
    neutrino_rate: np.ndarray
    process_gev: np.ndarray
    process_rate: np.ndarray
    pion_nu: np.ndarray
    muon_nu: np.ndarray
    muon_electron: np.ndarray
    hadron_gev: np.ndarray
    proton_reinj: np.ndarray
    neutron_reinj: np.ndarray
    proton_loss: np.ndarray
    neutron_loss: np.ndarray
    photon_gev: np.ndarray
    photon_loss: np.ndarray


def solve_pgamma(
    proton_gev: np.ndarray,
    proton_density: np.ndarray,
    photon_gev: np.ndarray,
    photon_density: np.ndarray,
    gamma_gev: np.ndarray,
    neutrino_gev: np.ndarray,
    process_gev: np.ndarray,
    neutron_density: np.ndarray | None = None,
) -> HummerOutput:
    ep = strict_grid(proton_gev, "proton_gev")
    protons = same_grid(proton_density, ep, "proton_density", nonnegative=True)
    neutrons = np.zeros_like(protons) if neutron_density is None else same_grid(
        neutron_density,
        ep,
        "neutron_density",
        nonnegative=True,
    )
    eph = strict_grid(photon_gev, "photon_gev")
    photons = same_grid(photon_density, eph, "photon_density", nonnegative=True)
    egamma = strict_grid(gamma_gev, "gamma_gev")
    enu = strict_grid(neutrino_gev, "neutrino_gev")
    eproc = strict_grid(process_gev, "process_gev")

    qpi0, qpip, qpim, qpro, qntr, afpro, afntr, afrad = hadronic_module.pg_operator(
        ep,
        protons,
        eph,
        photons,
        neutrons,
    )
    qpi0 = np.asarray(qpi0, dtype=float)
    qpip = np.asarray(qpip, dtype=float)
    qpim = np.asarray(qpim, dtype=float)

    (
        gamma_rate,
        process_rate,
        _,
        _,
        _,
        _,
        pionnu,
        muonnu,
        muone,
        neutrino_rate,
    ) = hadronic_module.decay_operator(
        ep,
        qpi0,
        qpip,
        qpim,
        egamma,
        enu,
        eproc,
    )

    return HummerOutput(
        gamma_gev=egamma,
        gamma_rate=np.asarray(gamma_rate, dtype=float),
        neutrino_gev=enu,
        neutrino_rate=np.asarray(neutrino_rate, dtype=float),
        process_gev=eproc,
        process_rate=np.asarray(process_rate, dtype=float),
        pion_nu=np.asarray(pionnu, dtype=float),
        muon_nu=np.asarray(muonnu, dtype=float),
        muon_electron=np.asarray(muone, dtype=float),
        hadron_gev=ep,
        proton_reinj=np.asarray(qpro, dtype=float),
        neutron_reinj=np.asarray(qntr, dtype=float),
        proton_loss=np.asarray(afpro, dtype=float),
        neutron_loss=np.asarray(afntr, dtype=float),
        photon_gev=eph,
        photon_loss=np.asarray(afrad, dtype=float),
    )
