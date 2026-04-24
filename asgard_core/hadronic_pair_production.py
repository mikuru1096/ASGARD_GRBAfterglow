from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

try:
    import src.Hadronic.FS_hadronic_1d as hadronic_fortran_module
except ImportError:
    hadronic_fortran_module = None


LIGHT_SPEED_CGS = 2.99792458e10
SIGMA_T_CGS = 6.6524587158e-25
ELECTRON_MASS_GEV = 5.1099895e-4
REG = 1.0e-50

T_SSC_CGS = math.sqrt(SIGMA_T_CGS) / LIGHT_SPEED_CGS
N_SSC_CGS = SIGMA_T_CGS ** (-1.5)
RATE_SSC_CGS = N_SSC_CGS / T_SSC_CGS

_A0V = np.zeros(2101, dtype=float)
_A_A0_HV = np.zeros(2101, dtype=float)
_HAS_FORTRAN_PAIR_PRODUCTION = hadronic_fortran_module is not None and hasattr(
    hadronic_fortran_module, "fs_hadronic_pair_production_shell"
)
PAIR_PRODUCTION_BACKEND = "fortran_pair_production" if _HAS_FORTRAN_PAIR_PRODUCTION else "unavailable"


@dataclass(frozen=True)
class PairProductionOutput:
    photon_energy_gev: np.ndarray
    electron_energy_gev: np.ndarray
    photon_loss_rate: np.ndarray
    pair_injection_rate_per_gev_per_species: np.ndarray
    pair_injection_rate_per_gev_total: np.ndarray
    absorbed_power_gev_per_cm3_s: float
    injected_power_gev_per_cm3_s: float


def solve_pair_production(
    photon_energy_gev: np.ndarray,
    photon_density_per_gev: np.ndarray,
    electron_energy_gev: np.ndarray,
    max_com_energy_factor: int = 138,
) -> PairProductionOutput:
    e_ph = _as_strictly_increasing(photon_energy_gev, "photon_energy_gev")
    n_ph = _as_matching_and_nonnegative(photon_density_per_gev, e_ph, "photon_density_per_gev")
    e_e = _as_strictly_increasing(electron_energy_gev, "electron_energy_gev")
    if int(max_com_energy_factor) < 1:
        raise ValueError("max_com_energy_factor must be >= 1.")

    if not _HAS_FORTRAN_PAIR_PRODUCTION:
        raise RuntimeError("Pair production core must be provided by the Fortran backend.")

    (
        photon_loss_rate,
        pair_rate_per_gev_per_species,
        pair_rate_per_gev_total,
        absorbed_power,
        injected_power,
    ) = hadronic_fortran_module.fs_hadronic_pair_production_shell(
        e_ph,
        n_ph,
        e_e,
        int(max_com_energy_factor),
    )
    return PairProductionOutput(
        photon_energy_gev=e_ph,
        electron_energy_gev=e_e,
        photon_loss_rate=np.asarray(photon_loss_rate, dtype=float),
        pair_injection_rate_per_gev_per_species=np.asarray(pair_rate_per_gev_per_species, dtype=float),
        pair_injection_rate_per_gev_total=np.asarray(pair_rate_per_gev_total, dtype=float),
        absorbed_power_gev_per_cm3_s=float(absorbed_power),
        injected_power_gev_per_cm3_s=float(injected_power),
    )


def _calc_pair_injection_full(
    gm_e: np.ndarray,
    x_ph: np.ndarray,
    pp_ng: np.ndarray,
    dln: float,
    ind_min_energy_pho: int,
    max_com_energy_factor: int,
    afpair_ssc: np.ndarray,
    photons_log_ssc: np.ndarray,
) -> np.ndarray:
    n_lep = gm_e.size
    n_ph = x_ph.size
    epspair_ssc = np.zeros(n_lep, dtype=float)

    for i in range(1, n_lep - 1):
        outer = _outer_pp(gm_e[i], dln, ind_min_energy_pho, n_ph)
        accum = 0.0
        for j in range(outer, n_ph):
            inner = min(_inner_pp(x_ph[j], gm_e[i], dln, ind_min_energy_pho, n_ph), n_ph - 1)
            kmax = min(-j + 2 * ind_min_energy_pho + max_com_energy_factor, n_ph)
            for k in range(inner, kmax):
                accum += _rgg_d1(gm_e[i], x_ph[j], x_ph[k]) * pp_ng[k] * pp_ng[j]
        epspair_ssc[i] = accum * dln * dln * 0.75 * gm_e[i]

    r_alpha_phot = float(np.sum(x_ph[ind_min_energy_pho:] * afpair_ssc[ind_min_energy_pho:] * photons_log_ssc[ind_min_energy_pho:]))
    r_eps_raw = float(np.sum(gm_e * epspair_ssc))
    if r_eps_raw > 1.0e-100:
        epspair_ssc *= 0.5 * r_alpha_phot / r_eps_raw
    epspair_ssc[0] = 0.0
    return epspair_ssc


def _build_photon_loss_kernel(x_ph: np.ndarray) -> np.ndarray:
    n_ph = x_ph.size
    ker = np.zeros((n_ph, n_ph), dtype=float)
    for i in range(n_ph):
        xi = x_ph[i]
        for j in range(n_ph):
            xj = x_ph[j]
            ker[i, j] = 0.375 * _phibar(xi, xj) / (xi * xi * xj * xj)
    return ker


def _phibar(a: float, b: float) -> float:
    s = a * b
    if s <= 1.0:
        return 0.0
    if s <= 1.1:
        return _phibar1(s)
    if s < 5.0:
        return _phibar2(s, -1.0 + 2.0 * s * (1.0 + math.sqrt(1.0 - 1.0 / s)))
    return _phibar3(s)


def _phibar1(s: float) -> float:
    s1 = s - 1.0
    s2 = s1 * math.sqrt(s1)
    return s2 * (1.333333 + 1.2 * s1 - 253.0 * s1 * s1 / 70.0)


def _phibar2(s: float, w: float) -> float:
    v = math.log(w)
    u = 1.0 - 1.0 / s
    return (
        (2.0 - 4.0 * s) * math.sqrt(u)
        + v * (4.0 * math.log(1.0 + w) - 3.0 * v + s * (1.0 + u * u))
        - 3.289868
        + _phisum(w)
    )


def _phibar3(s: float) -> float:
    w = math.log(4.0 * s)
    return (2.0 * s + w) * (w - 2.0) + (w + 1.125) / s - 0.289868


def _phisum(w: float) -> float:
    total = 0.0
    for i in range(1, 15):
        total += _sign_int(i) / (w**i * i * i)
    return -4.0 * total


def _inner_pp(x: float, gm: float, dln: float, ind_min_energy_pho: int, n_ph: int) -> int:
    if x <= 0.5:
        fval = _fpp_m(x, gm)
        if abs(fval) <= 1.0e-300:
            return 0
        arg = gm - x + 0.5 * (fval + 1.0 / fval)
        if arg <= 0.0:
            return 0
        q = int(math.log(arg) / dln + ind_min_energy_pho + 1)
        return min(max(q, 0), n_ph)
    if x < 1.0 and gm < _gm_b(x):
        fval = _fpp_m(x, gm)
        if abs(fval) <= 1.0e-300:
            return 0
        arg = gm - x + 0.5 * (fval + 1.0 / fval)
        if arg <= 0.0:
            return 0
        q = int(math.log(arg) / dln + ind_min_energy_pho + 1)
        return min(max(q, 0), n_ph)
    if x > 1.0 and gm < _gm_b(x):
        fval = _fpp_p(x, gm)
        if abs(fval) <= 1.0e-300:
            return 0
        arg = gm - x + 0.5 * (fval + 1.0 / fval)
        if arg <= 0.0:
            return 0
        q = int(math.log(arg) / dln + ind_min_energy_pho + 1)
        return min(max(q, 0), n_ph)
    return 0


def _outer_pp(gm: float, dln: float, ind_min_energy_pho: int, n_ph: int) -> int:
    q = int(math.log(_x_l(gm)) / dln + ind_min_energy_pho)
    return min(max(q, 0), n_ph - 1)


def _x_l(gm: float) -> float:
    return 0.5 * gm * (1.0 - _beta(gm))


def _gm_b(x: float) -> float:
    return x - (x - 1.0) / (2.0 * x - 1.0)


def _fpp_m(x: float, gp: float) -> float:
    return 2.0 * x - gp * (1.0 - _beta(gp))


def _fpp_p(x: float, gp: float) -> float:
    return 2.0 * x - gp * (1.0 + _beta(gp))


def _rgg_d1(gm: float, x: float, x1: float) -> float:
    gp = x + x1 - gm
    if gp <= 1.0 or gm <= 1.0:
        return 0.0

    if gp > 10.0 and gm > 10.0:
        kl = 0.25 / gp + 0.09375 / gp**3 + 0.25 / gm + 0.09375 / gm**3
        ku = 0.5 * (gp + math.sqrt(gp * gp - 1.0) + gm + math.sqrt(gm * gm - 1.0))
    elif gp > 10.0:
        kl = 0.5 * (0.5 / gp + gm - math.sqrt(gm * gm - 1.0))
        ku = gp + 0.5 * (gm + math.sqrt(gm * gm - 1.0))
    elif gm > 10.0:
        kl = 0.5 * (gp - math.sqrt(gp * gp - 1.0) + 0.5 / gm)
        ku = gm + 0.5 * (gp + math.sqrt(gp * gp - 1.0))
    else:
        kl = 0.5 * (gp - math.sqrt(gp * gp - 1.0) + gm - math.sqrt(gm * gm - 1.0))
        ku = 0.5 * (gp + math.sqrt(gp * gp - 1.0) + gm + math.sqrt(gm * gm - 1.0))

    if x < kl or x1 < kl or x > ku or x1 > ku:
        return 0.0

    tval = x + x1
    x_upper = _xcm_u(gm, x, x1)
    x_lower = _xcm_l(gm, x, x1)
    u = 4.0 * x_upper * x_upper / (tval * tval)
    v = 4.0 * x_lower * x_lower / (tval * tval)

    h_nh_u = ((gm - x) * (gm - x) - 1.0) * x_upper * x_upper / (x * x1)
    h_nh_l = ((gm - x) * (gm - x) - 1.0) * x_lower * x_lower / (x * x1)
    h_ih_u = ((gm - x1) * (gm - x1) - 1.0) * x_upper * x_upper / (x * x1)
    h_ih_l = ((gm - x1) * (gm - x1) - 1.0) * x_lower * x_lower / (x * x1)

    td2a = _td2(gm, x, x1, x_upper)
    td2b = _td2(gm, x, x1, x_lower)
    td2c = _td2(gm, x1, x, x_upper)
    td2d = _td2(gm, x1, x, x_lower)

    if h_nh_u > 1000.0 and h_nh_l > 1000.0:
        denom = math.sqrt((gm - x) * (gm - x) - 1.0)
        td2a = -(x_upper * x_upper + 2.0 * math.log(x_upper)) / denom
        td2b = -(x_lower * x_lower + 2.0 * math.log(x_lower)) / denom
    if h_ih_u > 1000.0 and h_ih_l > 1000.0:
        denom = math.sqrt((gm - x1) * (gm - x1) - 1.0)
        td2c = -(x_upper * x_upper + 2.0 * math.log(x_upper)) / denom
        td2d = -(x_lower * x_lower + 2.0 * math.log(x_lower)) / denom

    if u < 0.01 and v < 0.01:
        return -0.25 * (0.5 * tval * (u - v) + td2a - td2b + td2c - td2d)
    return -0.25 * (tval * (math.sqrt(1.0 - v) - math.sqrt(1.0 - u + 1.0e-9)) + td2a - td2b + td2c - td2d)


def _xcm_u(gm: float, x: float, x1: float) -> float:
    gp = x + x1 - gm
    if gp > 15.0 and gm > 15.0:
        return min(math.sqrt(x * x1), math.sqrt(gp * gm))
    return min(
        math.sqrt(x * x1),
        math.sqrt(0.5 * (gm * gp + 1.0 + math.sqrt(gm * gm - 1.0) * math.sqrt(gp * gp - 1.0))),
    )


def _xcm_l(gm: float, x: float, x1: float) -> float:
    gp = x + x1 - gm
    if gp > 10.0:
        if gm > 10.0:
            return 0.5 * math.sqrt(2.0 + gp / gm + gm / gp)
        if gm < 1.001:
            return 1.001
        return math.sqrt(0.5 * ((gm - math.sqrt(gm * gm - 1.0 + REG)) * (gp - 0.5 / gp) + 1.0))
    if gp < 1.001:
        return 1.001
    if gm > 10.0:
        return math.sqrt(0.5 * ((gp - math.sqrt(gp * gp - 1.0 + REG)) * (gm - 0.5 / gm) + 1.0))
    if gm < 1.001:
        return 1.001
    return math.sqrt(0.5 * (gp * gm + 1.0 - math.sqrt((gp * gp - 1.0) * (gm * gm - 1.0))))


def _td2(gm: float, x: float, x1: float, xcm: float) -> float:
    y = math.sqrt(x * x1)
    y2 = y * y
    q = xcm / y
    h = ((gm - x) * (gm - x) - 1.0) * q * q
    if h > 1000.0:
        return -(
            xcm * xcm + h / (xcm * xcm) + 0.886294 - 0.5 * (gm * (x1 - x) + x * x) / y2 + math.log(h)
        ) * q / math.sqrt(h)
    if h < -0.999:
        return q * (q * q * (y2 - 1.0) * (-1.5707963) + 0.5 * ((x1 - x) / y2 - 6.28319))

    ah = math.sqrt(1.0 + h)
    return q * q * q * (y2 - 1.0) * _a_a0_h(h) - ah / (xcm * y) + 0.5 * q * (
        (1.0 + (gm * (x1 - x) + x * x) / y2 - 2.0 * q * q) / ah - 4.0 * _a0(h)
    )


def _a_a0_h(h: float) -> float:
    if h <= 20.0:
        i = 100 + int(h / 0.01)
        return _A_A0_HV[i]
    return -(1.0 - (0.693147 + 0.5 * math.log(h)) / h) / math.sqrt(h)


def _a0(h: float) -> float:
    if h < 20.0:
        i = 100 + int(h / 0.01)
        return _A0V[i]
    return (0.693147 + 0.5 * math.log(h)) / math.sqrt(h)


def _a0f(h: float) -> float:
    if h > 1.0e-6:
        return math.log(math.sqrt(h) + math.sqrt(1.0 + h)) / math.sqrt(h)
    if h < -1.0e-6:
        return math.asin(math.sqrt(-h) - 1.0e-8) / math.sqrt(-h)
    return 1.0 - h / 6.0


def _a_a0_hf(h: float) -> float:
    if -0.1 <= h <= 0.1:
        return -0.666667 + (0.2 - (0.107143 - 0.0694444 * h) * h) * h
    return (_a0f(h) - math.sqrt(1.0 + h + 1.0e-8)) / h


def _beta(gamma: float) -> float:
    if gamma < 10.0:
        return math.sqrt(1.0 - 1.0 / (gamma * gamma) + REG)
    return 1.0 - 0.5 / (gamma * gamma) * (1.0 + 0.25 / (gamma * gamma))


def _sign_int(i: int) -> int:
    return 1 if (i % 2 == 0) else -1


def _as_strictly_increasing(values: np.ndarray, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be a 1d array.")
    if arr.size < 2:
        raise ValueError(f"{name} must contain at least two points.")
    if np.any(arr <= 0.0):
        raise ValueError(f"{name} must be strictly positive.")
    if np.any(np.diff(arr) <= 0.0):
        raise ValueError(f"{name} must be strictly increasing.")
    return arr


def _as_matching_and_nonnegative(values: np.ndarray, grid: np.ndarray, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.shape != grid.shape:
        raise ValueError(f"{name} must match grid shape.")
    if np.any(arr < 0.0):
        raise ValueError(f"{name} must be non-negative.")
    return arr


def _log_spacing(grid: np.ndarray, name: str) -> float:
    dln = np.diff(np.log(grid))
    if not np.allclose(dln, dln[0], rtol=1.0e-6, atol=1.0e-12):
        raise ValueError(f"{name} must be logarithmically uniform.")
    return float(dln[0])


def _grid_index_offset(e_min: float, ph_min: float, dln: float) -> int:
    ratio = math.log(e_min / ph_min) / dln
    rounded = int(round(ratio))
    if not math.isclose(ratio, float(rounded), rel_tol=1.0e-9, abs_tol=1.0e-9):
        raise ValueError("Grid minima are not aligned on the shared logarithmic lattice.")
    return rounded


for _idx in range(2101):
    _h = -1.0 + 0.01 * _idx
    _A0V[_idx] = _a0f(_h)
    _A_A0_HV[_idx] = _a_a0_hf(_h)
