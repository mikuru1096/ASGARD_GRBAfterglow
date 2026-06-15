from __future__ import annotations

import numpy as np

import numpy as np
import src.Hadronic.hadronic_forward_1d as hadronic_legacy_module
from src import constants


MPI0_MASS_GEV = 0.137
PROTON_MASS_GEV = constants.para_m_p_e * constants.para_erg2ev * 1.0e-9
ETA0 = 2.0 * MPI0_MASS_GEV / PROTON_MASS_GEV + (MPI0_MASS_GEV / PROTON_MASS_GEV) ** 2
R_PI = MPI0_MASS_GEV / PROTON_MASS_GEV
ETA2 = 4.0 * R_PI * (1.0 + R_PI)

_TABLE1_RHO = np.array([
    1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
    3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 20.0, 30.0,
    40.0, 100.0,
], dtype=float)
_TABLE1_S = np.array([
    0.0768, 0.106, 0.182, 0.201, 0.219, 0.216, 0.233, 0.233, 0.248, 0.244,
    0.188, 0.131, 0.120, 0.107, 0.102, 0.0932, 0.0838, 0.0761, 0.107, 0.0928,
    0.0772, 0.0479,
], dtype=float)
_TABLE1_DELTA = np.array([
    0.544, 0.540, 0.750, 0.791, 0.788, 0.831, 0.839, 0.825, 0.805, 0.779,
    1.23, 1.82, 2.05, 2.19, 2.23, 2.29, 2.37, 2.43, 2.27, 2.33,
    2.42, 2.59,
], dtype=float)
_TABLE1_B = np.array([
    2.86e-19, 2.24e-18, 5.61e-18, 1.02e-17, 1.60e-17, 2.23e-17, 3.10e-17, 4.07e-17, 5.30e-17, 6.74e-17,
    1.51e-16, 1.24e-16, 1.37e-16, 1.62e-16, 1.71e-16, 1.78e-16, 1.84e-16, 1.93e-16, 4.74e-16, 7.70e-16,
    1.06e-15, 2.73e-15,
], dtype=float)

_TABLE2_RHO = np.array([
    1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
    3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 30.0, 100.0,
], dtype=float)
_TABLE2_EPLUS = (
    np.array([0.367, 0.282, 0.260, 0.239, 0.224, 0.207, 0.198, 0.193, 0.187, 0.181, 0.122, 0.106, 0.0983, 0.0875, 0.0830, 0.0783, 0.0735, 0.0644, 0.0333, 0.0224], dtype=float),
    np.array([3.12, 2.96, 2.83, 2.76, 2.69, 2.66, 2.62, 2.56, 2.52, 2.49, 2.48, 2.50, 2.46, 2.46, 2.44, 2.44, 2.45, 2.50, 2.77, 2.86], dtype=float),
    np.array([8.09e-19, 7.70e-18, 2.05e-17, 3.66e-17, 5.48e-17, 7.39e-17, 9.52e-17, 1.20e-16, 1.47e-16, 1.75e-16, 3.31e-16, 4.16e-16, 5.57e-16, 6.78e-16, 7.65e-16, 8.52e-16, 9.17e-16, 9.57e-16, 3.07e-15, 1.58e-14], dtype=float),
)
_TABLE2_NUMU_BAR = (
    np.array([0.365, 0.287, 0.250, 0.238, 0.220, 0.206, 0.197, 0.193, 0.187, 0.178, 0.123, 0.106, 0.0944, 0.0829, 0.0801, 0.0752, 0.0680, 0.0615, 0.0361, 0.0228], dtype=float),
    np.array([3.09, 2.96, 2.89, 2.76, 2.71, 2.67, 2.62, 2.56, 2.52, 2.51, 2.48, 2.56, 2.57, 2.58, 2.54, 2.53, 2.56, 2.60, 2.78, 2.88], dtype=float),
    np.array([8.09e-19, 7.70e-18, 1.99e-17, 3.62e-17, 5.39e-17, 7.39e-17, 9.48e-17, 1.20e-16, 1.47e-16, 1.74e-16, 3.38e-16, 5.17e-16, 7.61e-16, 9.57e-16, 1.11e-15, 1.25e-15, 1.36e-15, 1.46e-15, 5.87e-15, 3.10e-14], dtype=float),
)
_TABLE2_NUMU = (
    np.array([0.0, 0.0778, 0.242, 0.377, 0.440, 0.450, 0.461, 0.451, 0.464, 0.446, 0.366, 0.249, 0.204, 0.174, 0.156, 0.140, 0.121, 0.107, 0.0705, 0.0463], dtype=float),
    np.array([0.0, 0.306, 0.792, 1.09, 1.06, 0.953, 0.956, 0.922, 0.912, 0.940, 1.49, 2.03, 2.18, 2.24, 2.28, 2.32, 2.39, 2.46, 2.53, 2.62], dtype=float),
    np.array([1.08e-18, 9.91e-18, 2.47e-17, 4.43e-17, 6.70e-17, 9.04e-17, 1.18e-16, 1.32e-16, 1.77e-16, 2.11e-16, 3.83e-16, 5.09e-16, 7.26e-16, 9.26e-16, 1.07e-15, 1.19e-15, 1.29e-15, 1.40e-15, 5.65e-15, 3.01e-14], dtype=float),
)
_TABLE2_NE = (
    np.array([0.768, 0.569, 0.491, 0.395, 0.31, 0.323, 0.305, 0.285, 0.270, 0.259, 0.158, 0.129, 0.113, 0.0996, 0.0921, 0.0861, 0.0800, 0.0723, 0.0411, 0.0283], dtype=float),
    np.array([2.49, 2.35, 2.41, 2.45, 2.45, 2.43, 2.40, 2.39, 2.37, 2.35, 2.42, 2.46, 2.45, 2.46, 2.46, 2.45, 2.47, 2.51, 2.70, 2.77], dtype=float),
    np.array([9.43e-19, 9.22e-18, 2.35e-17, 4.20e-17, 6.26e-17, 8.57e-17, 1.13e-16, 1.39e-16, 1.70e-16, 2.05e-16, 3.81e-16, 4.74e-16, 6.30e-16, 7.65e-16, 8.61e-16, 9.61e-16, 1.03e-15, 1.10e-15, 3.55e-15, 1.86e-14], dtype=float),
)

_TABLE3_RHO = np.array([3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 30.0, 100.0], dtype=float)
_TABLE3_EMINUS = (
    np.array([0.658, 0.348, 0.286, 0.256, 0.258, 0.220, 0.217, 0.192, 0.125, 0.0507], dtype=float),
    np.array([3.09, 2.81, 2.39, 2.27, 2.13, 2.20, 2.13, 2.19, 2.27, 2.63], dtype=float),
    np.array([6.43e-19, 9.91e-18, 1.24e-16, 2.67e-16, 3.50e-16, 4.03e-16, 4.48e-16, 4.78e-16, 1.64e-15, 4.52e-15], dtype=float),
)
_TABLE3_NUEBAR = (
    np.array([0.985, 0.378, 0.31, 0.327, 0.308, 0.292, 0.260, 0.233, 0.135, 0.0770], dtype=float),
    np.array([2.63, 2.98, 2.31, 2.11, 2.03, 1.98, 2.02, 2.07, 2.24, 2.40], dtype=float),
    np.array([6.61e-19, 9.74e-18, 1.34e-16, 2.91e-16, 3.81e-16, 4.48e-16, 4.83e-16, 5.13e-16, 1.75e-15, 5.48e-15], dtype=float),
)


def photon_density_hz_to_gev(photon_nu_hz: np.ndarray, photon_density_per_hz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Convert `n_nu` to `n_E` using `E = h nu`.

    If `n_nu` is given per Hz, then `n_E(E) = n_nu(nu) / h_GeV`.
    """
    nu = np.asarray(photon_nu_hz, dtype=float)
    density = np.asarray(photon_density_per_hz, dtype=float)
    if nu.ndim != 1 or density.ndim != 1 or nu.shape != density.shape:
        raise ValueError("photon_nu_hz and photon_density_per_hz must be 1d arrays with the same shape.")
    if np.any(nu <= 0.0):
        raise ValueError("photon_nu_hz must be positive.")
    photon_energy_gev, photon_density_per_gev = hadronic_legacy_module.fs_hadronic_photon_density_hz_to_gev(
        nu,
        density,
    )
    return np.asarray(photon_energy_gev, dtype=float), np.asarray(photon_density_per_gev, dtype=float)


def kelner_aharonian_2008_gamma_phi(eta: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Gamma-channel stable-secondary parameterization from Kelner & Aharonian (2008)."""
    eta_arr = np.asarray(eta, dtype=float)
    x_arr = np.asarray(x, dtype=float)
    eta_b, x_b = np.broadcast_arrays(eta_arr, x_arr)
    out = np.zeros_like(eta_b)

    valid_eta = eta_b >= ETA0
    if not np.any(valid_eta):
        return out

    rho = eta_b[valid_eta] / ETA0
    _require_table_range(rho, _TABLE1_RHO[0], _TABLE1_RHO[-1], "gamma")
    b_gamma = _interp_table_linear(_TABLE1_RHO, _TABLE1_B, rho)
    s_gamma = _interp_table_linear(_TABLE1_RHO, _TABLE1_S, rho)
    delta_gamma = _interp_table_linear(_TABLE1_RHO, _TABLE1_DELTA, rho)
    x_minus, x_plus = _x_pm_gamma(eta_b[valid_eta])
    psi = 2.5 + 0.4 * np.log(rho)

    xv = x_b[valid_eta]
    low = xv < x_minus
    mid = (xv >= x_minus) & (xv <= x_plus)

    loc = np.zeros_like(xv)
    if np.any(low):
        loc[low] = b_gamma[low] * np.power(np.log(2.0), psi[low])
    if np.any(mid):
        y = (xv[mid] - x_minus[mid]) / np.maximum(x_plus[mid] - x_minus[mid], 1.0e-300)
        shape = np.log(2.0 / (1.0 + y * y))
        expo = -s_gamma[mid] * np.power(np.log(xv[mid] / x_minus[mid]), delta_gamma[mid])
        loc[mid] = b_gamma[mid] * np.exp(expo) * np.power(shape, psi[mid])
    out[valid_eta] = loc
    return out


def kelner_aharonian_2008_gamma_spectrum(
    gamma_energy_gev: np.ndarray,
    proton_energy_gev: np.ndarray,
    proton_density_per_gev: np.ndarray,
    photon_energy_gev: np.ndarray,
    photon_density_per_gev: np.ndarray,
) -> np.ndarray:
    """Eq. (10) gamma production spectrum for arbitrary proton and photon distributions."""
    e_gamma = np.asarray(gamma_energy_gev, dtype=float)
    e_p = np.asarray(proton_energy_gev, dtype=float)
    f_p = np.asarray(proton_density_per_gev, dtype=float)
    eps = np.asarray(photon_energy_gev, dtype=float)
    f_ph = np.asarray(photon_density_per_gev, dtype=float)
    if e_p.ndim != 1 or f_p.ndim != 1 or e_p.shape != f_p.shape:
        raise ValueError("proton_energy_gev and proton_density_per_gev must be 1d with matching shapes.")
    if eps.ndim != 1 or f_ph.ndim != 1 or eps.shape != f_ph.shape:
        raise ValueError("photon_energy_gev and photon_density_per_gev must be 1d with matching shapes.")
    if np.any(e_gamma <= 0.0) or np.any(e_p <= PROTON_MASS_GEV) or np.any(eps <= 0.0):
        raise ValueError("All energies must be positive and proton energy must exceed rest mass.")

    spectrum = np.zeros_like(e_gamma)
    for i, eg in enumerate(e_gamma):
        outer = np.zeros_like(e_p)
        for j, ep in enumerate(e_p):
            x = eg / ep
            if x <= 0.0:
                continue
            eta = 4.0 * eps * ep / (PROTON_MASS_GEV * PROTON_MASS_GEV)
            phi = kelner_aharonian_2008_gamma_phi(eta, x)
            inner = np.trapezoid(f_ph * phi, eps)
            outer[j] = f_p[j] * inner / ep
        spectrum[i] = np.trapezoid(outer, e_p)
    return spectrum


def kelner_aharonian_2008_secondary_phi(channel: str, eta: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Stable-secondary parameterization for KA2008 photomeson fits."""
    key = _normalize_channel(channel)
    if key == "gamma":
        return kelner_aharonian_2008_gamma_phi(eta, x)
    if key in {"e_plus", "nu_mu_bar", "nu_e"}:
        return _ka2008_table2_phi(key, eta, x)
    if key == "nu_mu":
        return _ka2008_numu_phi(eta, x)
    if key in {"e_minus", "nu_e_bar"}:
        return _ka2008_table3_phi(key, eta, x)
    raise ValueError(f"Unsupported KA2008 secondary channel: {channel!r}")


def kelner_aharonian_2008_secondary_spectrum(
    channel: str,
    secondary_energy_gev: np.ndarray,
    proton_energy_gev: np.ndarray,
    proton_density_per_gev: np.ndarray,
    photon_energy_gev: np.ndarray,
    photon_density_per_gev: np.ndarray,
) -> np.ndarray:
    """General KA2008 stable-secondary spectrum for a named channel."""
    e_sec = np.asarray(secondary_energy_gev, dtype=float)
    e_p = np.asarray(proton_energy_gev, dtype=float)
    f_p = np.asarray(proton_density_per_gev, dtype=float)
    eps = np.asarray(photon_energy_gev, dtype=float)
    f_ph = np.asarray(photon_density_per_gev, dtype=float)
    if e_p.ndim != 1 or f_p.ndim != 1 or e_p.shape != f_p.shape:
        raise ValueError("proton_energy_gev and proton_density_per_gev must be 1d with matching shapes.")
    if eps.ndim != 1 or f_ph.ndim != 1 or eps.shape != f_ph.shape:
        raise ValueError("photon_energy_gev and photon_density_per_gev must be 1d with matching shapes.")
    if np.any(e_sec <= 0.0) or np.any(e_p <= PROTON_MASS_GEV) or np.any(eps <= 0.0):
        raise ValueError("All energies must be positive and proton energy must exceed rest mass.")

    spectrum = np.zeros_like(e_sec)
    for i, es in enumerate(e_sec):
        outer = np.zeros_like(e_p)
        for j, ep in enumerate(e_p):
            x = es / ep
            if x <= 0.0:
                continue
            eta = 4.0 * eps * ep / (PROTON_MASS_GEV * PROTON_MASS_GEV)
            phi = kelner_aharonian_2008_secondary_phi(channel, eta, x)
            inner = np.trapezoid(f_ph * phi, eps)
            outer[j] = f_p[j] * inner / ep
        spectrum[i] = np.trapezoid(outer, e_p)
    return spectrum


def _interp_table_linear(x_tab: np.ndarray, y_tab: np.ndarray, x: np.ndarray) -> np.ndarray:
    x_arr = np.asarray(x, dtype=float)
    return np.interp(x_arr, x_tab, y_tab)


def _require_table_range(x: np.ndarray, xmin: float, xmax: float, label: str) -> None:
    if np.any(x < xmin) or np.any(x > xmax):
        raise ValueError(
            f"KA2008 {label} fit is only tabulated for eta/eta0 in [{xmin}, {xmax}]. "
            "Use a literature-based extension before evaluating outside this domain."
        )


def _normalize_channel(channel: str) -> str:
    key = channel.strip().lower().replace("-", "_")
    aliases = {
        "gamma": "gamma",
        "g": "gamma",
        "eplus": "e_plus",
        "e_plus": "e_plus",
        "positron": "e_plus",
        "nuebar": "nu_e_bar",
        "nu_ebar": "nu_e_bar",
        "nu_e_bar": "nu_e_bar",
        "antinue": "nu_e_bar",
        "numubar": "nu_mu_bar",
        "nu_mubar": "nu_mu_bar",
        "nu_mu_bar": "nu_mu_bar",
        "antinumu": "nu_mu_bar",
        "numu": "nu_mu",
        "nu_mu": "nu_mu",
        "e_minus": "e_minus",
        "eminus": "e_minus",
        "electron": "e_minus",
        "nue": "nu_e",
        "nu_e": "nu_e",
    }
    if key not in aliases:
        raise ValueError(f"Unsupported KA2008 channel: {channel!r}")
    return aliases[key]


def _x_pm_gamma(eta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    eta_arr = np.asarray(eta, dtype=float)
    radicand = np.maximum((eta_arr - R_PI * R_PI - 2.0 * R_PI) * (eta_arr - R_PI * R_PI + 2.0 * R_PI), 0.0)
    root = np.sqrt(radicand)
    pref = 1.0 / (2.0 * (1.0 + eta_arr))
    x_plus = pref * (eta_arr + R_PI * R_PI + root)
    x_minus = pref * (eta_arr + R_PI * R_PI - root)
    return x_minus, x_plus


def _ka2008_table2_phi(channel: str, eta: np.ndarray, x: np.ndarray) -> np.ndarray:
    eta_arr = np.asarray(eta, dtype=float)
    x_arr = np.asarray(x, dtype=float)
    eta_b, x_b = np.broadcast_arrays(eta_arr, x_arr)
    out = np.zeros_like(eta_b)

    valid = eta_b >= ETA0
    if not np.any(valid):
        return out

    rho = eta_b[valid] / ETA0
    _require_table_range(rho, _TABLE2_RHO[0], _TABLE2_RHO[-1], channel)
    s_tab, d_tab, b_tab = _table2_series(channel)
    s_l = _interp_table_linear(_TABLE2_RHO, s_tab, rho)
    d_l = _interp_table_linear(_TABLE2_RHO, d_tab, rho)
    b_l = _interp_table_linear(_TABLE2_RHO, b_tab, rho)
    x_minus, x_plus = _x_table2_bounds(channel, eta_b[valid])
    psi = 2.5 + 1.4 * np.log(rho)
    out[valid] = _build_secondary_phi_piece(x_b[valid], x_minus, x_plus, s_l, d_l, b_l, psi)
    return out


def _ka2008_numu_phi(eta: np.ndarray, x: np.ndarray) -> np.ndarray:
    eta_arr = np.asarray(eta, dtype=float)
    x_arr = np.asarray(x, dtype=float)
    eta_b, x_b = np.broadcast_arrays(eta_arr, x_arr)
    out = np.zeros_like(eta_b)

    valid = eta_b >= ETA0
    if not np.any(valid):
        return out

    rho = eta_b[valid] / ETA0
    _require_table_range(rho, _TABLE2_RHO[0], _TABLE2_RHO[-1], "nu_mu")
    s_l, d_l, b_l = _table2_series("nu_mu")
    s_l = _interp_table_linear(_TABLE2_RHO, s_l, rho)
    d_l = _interp_table_linear(_TABLE2_RHO, d_l, rho)
    b_l = _interp_table_linear(_TABLE2_RHO, b_l, rho)
    x_minus_gamma, x_plus_gamma = _x_pm_gamma(eta_b[valid])
    x_minus = 0.427 * x_minus_gamma
    x_plus = np.where(
        rho < 2.14,
        0.427 * x_plus_gamma,
        np.where(
            rho < 10.0,
            (0.427 + 0.0729 * (rho - 2.14)) * x_plus_gamma,
            x_plus_gamma,
        ),
    )
    psi = 2.5 + 1.4 * np.log(rho)
    out[valid] = _build_secondary_phi_piece(x_b[valid], x_minus, x_plus, s_l, d_l, b_l, psi)
    return out


def _ka2008_table3_phi(channel: str, eta: np.ndarray, x: np.ndarray) -> np.ndarray:
    eta_arr = np.asarray(eta, dtype=float)
    x_arr = np.asarray(x, dtype=float)
    eta_b, x_b = np.broadcast_arrays(eta_arr, x_arr)
    out = np.zeros_like(eta_b)

    valid = eta_b >= ETA2
    if not np.any(valid):
        return out

    rho = eta_b[valid] / ETA0
    _require_table_range(rho, _TABLE3_RHO[0], _TABLE3_RHO[-1], channel)
    s_tab, d_tab, b_tab = _table3_series(channel)
    s_l = _interp_table_linear(_TABLE3_RHO, s_tab, rho)
    d_l = _interp_table_linear(_TABLE3_RHO, d_tab, rho)
    b_l = _interp_table_linear(_TABLE3_RHO, b_tab, rho)
    x_minus, x_plus = _x_two_pion_bounds(eta_b[valid])
    x_minus = 0.5 * x_minus
    psi = 6.0 * (1.0 - np.exp(1.5 * (4.0 - rho)))
    psi = np.where(rho >= 4.0, psi, 0.0)
    out[valid] = _build_secondary_phi_piece(x_b[valid], x_minus, x_plus, s_l, d_l, b_l, psi)
    return out


def _build_secondary_phi_piece(
    x: np.ndarray,
    x_minus: np.ndarray,
    x_plus: np.ndarray,
    s_l: np.ndarray,
    d_l: np.ndarray,
    b_l: np.ndarray,
    psi: np.ndarray,
) -> np.ndarray:
    out = np.zeros_like(x)
    low = x < x_minus
    mid = (x >= x_minus) & (x < x_plus)
    if np.any(low):
        out[low] = b_l[low] * np.power(np.log(2.0), psi[low])
    if np.any(mid):
        y = (x[mid] - x_minus[mid]) / np.maximum(x_plus[mid] - x_minus[mid], 1.0e-300)
        shape = np.log(2.0 / (1.0 + y * y))
        expo = -s_l[mid] * np.power(np.log(x[mid] / x_minus[mid]), d_l[mid])
        out[mid] = b_l[mid] * np.exp(expo) * np.power(shape, psi[mid])
    return out


def _table2_series(channel: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    key = _normalize_channel(channel)
    if key == "e_plus":
        return _TABLE2_EPLUS
    if key == "nu_mu_bar":
        return _TABLE2_NUMU_BAR
    if key == "nu_mu":
        return _TABLE2_NUMU
    if key == "nu_e":
        return _TABLE2_NE
    raise ValueError(f"Unsupported Table 2 channel: {channel!r}")


def _table3_series(channel: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    key = _normalize_channel(channel)
    if key == "e_minus":
        return _TABLE3_EMINUS
    if key == "nu_e_bar":
        return _TABLE3_NUEBAR
    raise ValueError(f"Unsupported Table 3 channel: {channel!r}")


def _x_table2_bounds(channel: str, eta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    key = _normalize_channel(channel)
    x_minus_gamma, x_plus_gamma = _x_pm_gamma(eta)
    if key in {"e_plus", "nu_mu_bar", "nu_e"}:
        return 0.25 * x_minus_gamma, x_plus_gamma
    if key == "nu_mu":
        raise ValueError("nu_mu does not use Table2 bounds helper.")
    raise ValueError(f"Unsupported Table2 channel: {channel!r}")


def _x_two_pion_bounds(eta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    eta_arr = np.asarray(eta, dtype=float)
    radicand = np.maximum(eta_arr * (eta_arr - ETA2), 0.0)
    root = np.sqrt(radicand)
    pref = 1.0 / (2.0 * (1.0 + eta_arr))
    x_plus = pref * (eta_arr - 2.0 * R_PI + root)
    x_minus = pref * (eta_arr - 2.0 * R_PI - root)
    return x_minus, x_plus
