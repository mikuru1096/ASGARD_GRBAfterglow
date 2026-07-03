from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from src import constants

import src.Hadronic.hadronic_forward_1d as hadronic_fortran_module


HUMMER2010_OPERATOR_BACKEND = "fortran_operator"
HUMMER2010_DECAY_BACKEND = "fortran_decay"
ACCELERATION_BACKEND = "fortran_acceleration"
BETHE_HEITLER_BACKEND = "fortran_bethe_heitler"
HADRONIC_IC_BACKEND = "fortran_hadronic_ic"
PP_DELTA_BACKEND = "fortran_pp_shell"
SECONDARY_RADIATION_BACKEND = "fortran_secondary_radiation"
SPECIES_TRANSPORT_BACKEND = "fortran_species_transport"


ERG_PER_GEV = constants.para_gev2erg
GEV_TO_ERG = ERG_PER_GEV
ELECTRON_MASS_GEV = constants.para_m_e_gev
PROTON_MASS_GEV = constants.para_m_p_gev
NEUTRON_MASS_GEV = constants.para_m_n_gev
PI0_MASS_GEV = constants.para_m_pi0_gev
PI_PLUS_MASS_GEV = constants.para_m_pi_charged_gev
PION_CHARGED_MASS_GEV = PI_PLUS_MASS_GEV
MUON_MASS_GEV = constants.para_m_mu_gev
CHARGED_PION_DECAY_S = 2.6033e-8
MUON_DECAY_S = 2.1969811e-6
AM3_C_CGS = constants.para_c
AM3_SIGMA_T_CGS = constants.para_sigmat
AM3_MASS_ELECTRON_GEV = ELECTRON_MASS_GEV
AM3_MASS_PROTON_GEV = PROTON_MASS_GEV
AM3_MASS_PION_CHARGED_GEV = PION_CHARGED_MASS_GEV
AM3_MASS_MUON_GEV = MUON_MASS_GEV


def as_strictly_increasing_grid(values: np.ndarray, name: str, *, require_finite: bool = False) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be a 1d array.")
    if arr.size < 2:
        raise ValueError(f"{name} must contain at least 2 points.")
    if np.any(arr <= 0.0):
        raise ValueError(f"{name} must be strictly positive.")
    if np.any(np.diff(arr) <= 0.0):
        raise ValueError(f"{name} must be strictly increasing.")
    if require_finite and not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain finite values.")
    return arr


def as_matching(values: np.ndarray, reference: np.ndarray, name: str, *, require_finite: bool = False) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.shape != reference.shape:
        raise ValueError(f"{name} must match grid shape.")
    if require_finite and not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain finite values.")
    return arr


def as_matching_nonnegative(values: np.ndarray, reference: np.ndarray, name: str) -> np.ndarray:
    arr = as_matching(values, reference, name)
    if np.any(arr < 0.0):
        raise ValueError(f"{name} must be non-negative.")
    return arr


def photon_density_hz_to_gev(photon_nu_hz: np.ndarray, photon_density_per_hz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    nu = np.asarray(photon_nu_hz, dtype=float)
    density = np.asarray(photon_density_per_hz, dtype=float)
    if nu.ndim != 1 or density.ndim != 1 or nu.shape != density.shape:
        raise ValueError("photon_nu_hz and photon_density_per_hz must be 1d arrays with the same shape.")
    if np.any(nu <= 0.0):
        raise ValueError("photon_nu_hz must be positive.")
    photon_energy_gev = constants.para_h_gev * nu
    photon_density_per_gev = density / constants.para_h_gev
    return photon_energy_gev, photon_density_per_gev


def positive_loglog_interp(x_src: np.ndarray, y_src: np.ndarray, x_dst: np.ndarray) -> np.ndarray:
    x = np.asarray(x_src, dtype=float)
    y = np.asarray(y_src, dtype=float)
    target = np.asarray(x_dst, dtype=float)
    out = np.zeros_like(target, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y) & (x > 0.0) & (y > 0.0)
    if np.count_nonzero(valid) < 2:
        return out
    xv = x[valid]
    yv = y[valid]
    in_range = np.isfinite(target) & (target >= xv[0]) & (target <= xv[-1])
    out[in_range] = np.exp(np.interp(np.log(target[in_range]), np.log(xv), np.log(yv)))
    return out


@dataclass(frozen=True)
class Hummer2010Output:
    gamma_energy_gev: np.ndarray
    gamma_rate_per_gev: np.ndarray
    neutrino_energy_gev: np.ndarray
    neutrino_rate_per_gev: np.ndarray
    process_energy_gev: np.ndarray
    process_rate_per_gev: np.ndarray
    pion_plus_source_rate_per_gev: np.ndarray
    pion_minus_source_rate_per_gev: np.ndarray
    muon_plus_right_source_rate_per_gev: np.ndarray
    muon_plus_left_source_rate_per_gev: np.ndarray
    muon_minus_left_source_rate_per_gev: np.ndarray
    muon_minus_right_source_rate_per_gev: np.ndarray
    prompt_pion_neutrino_rate_per_gev: np.ndarray
    muon_neutrino_rate_per_gev: np.ndarray
    muon_electron_rate_per_gev: np.ndarray
    hadron_energy_gev: np.ndarray
    proton_reinjection_rate_per_gev: np.ndarray
    neutron_reinjection_rate_per_gev: np.ndarray
    proton_loss_rate: np.ndarray
    neutron_loss_rate: np.ndarray
    photon_loss_energy_gev: np.ndarray
    photon_loss_rate: np.ndarray


def solve_hummer2010_pgamma(
    proton_energy_gev: np.ndarray,
    proton_density_per_gev: np.ndarray,
    photon_energy_gev: np.ndarray,
    photon_density_per_gev: np.ndarray,
    gamma_energy_gev: np.ndarray,
    neutrino_energy_gev: np.ndarray,
    process_energy_gev: np.ndarray,
    neutron_density_per_gev: np.ndarray | None = None,
) -> Hummer2010Output:
    e_p = as_strictly_increasing_grid(proton_energy_gev, "proton_energy_gev")
    n_p = as_matching_nonnegative(proton_density_per_gev, e_p, "proton_density_per_gev")
    n_n = np.zeros_like(n_p) if neutron_density_per_gev is None else as_matching_nonnegative(
        neutron_density_per_gev, e_p, "neutron_density_per_gev"
    )
    e_gam_t = as_strictly_increasing_grid(photon_energy_gev, "photon_energy_gev")
    n_gam_t = as_matching_nonnegative(photon_density_per_gev, e_gam_t, "photon_density_per_gev")
    e_gamma = as_strictly_increasing_grid(gamma_energy_gev, "gamma_energy_gev")
    e_nu = as_strictly_increasing_grid(neutrino_energy_gev, "neutrino_energy_gev")
    e_proc = as_strictly_increasing_grid(process_energy_gev, "process_energy_gev")

    (
        qpi0_rate,
        qpip_rate,
        qpim_rate,
        qpro_reinj_rate,
        qntr_reinj_rate,
        afpg_pro,
        afpg_ntr,
        afpg_rad,
    ) = hadronic_fortran_module.pg_operator(
        e_p,
        n_p,
        e_gam_t,
        n_gam_t,
        n_n,
    )
    qpi0_rate = np.asarray(qpi0_rate, dtype=float)
    qpip_rate = np.asarray(qpip_rate, dtype=float)
    qpim_rate = np.asarray(qpim_rate, dtype=float)

    (
        gamma_rate,
        process_rate,
        mu_plus_r_rate,
        mu_plus_l_rate,
        mu_minus_l_rate,
        mu_minus_r_rate,
        prompt_pion_neutrino_rate,
        muon_neutrino_rate,
        muon_electron_rate,
        neutrino_rate,
    ) = hadronic_fortran_module.decay_operator(
        e_p,
        qpi0_rate,
        qpip_rate,
        qpim_rate,
        e_gamma,
        e_nu,
        e_proc,
    )

    return Hummer2010Output(
        gamma_energy_gev=e_gamma,
        gamma_rate_per_gev=np.asarray(gamma_rate, dtype=float),
        neutrino_energy_gev=e_nu,
        neutrino_rate_per_gev=np.asarray(neutrino_rate, dtype=float),
        process_energy_gev=e_proc,
        process_rate_per_gev=np.asarray(process_rate, dtype=float),
        pion_plus_source_rate_per_gev=qpip_rate,
        pion_minus_source_rate_per_gev=qpim_rate,
        muon_plus_right_source_rate_per_gev=np.asarray(mu_plus_r_rate, dtype=float),
        muon_plus_left_source_rate_per_gev=np.asarray(mu_plus_l_rate, dtype=float),
        muon_minus_left_source_rate_per_gev=np.asarray(mu_minus_l_rate, dtype=float),
        muon_minus_right_source_rate_per_gev=np.asarray(mu_minus_r_rate, dtype=float),
        prompt_pion_neutrino_rate_per_gev=np.asarray(prompt_pion_neutrino_rate, dtype=float),
        muon_neutrino_rate_per_gev=np.asarray(muon_neutrino_rate, dtype=float),
        muon_electron_rate_per_gev=np.asarray(muon_electron_rate, dtype=float),
        hadron_energy_gev=e_p,
        proton_reinjection_rate_per_gev=np.asarray(qpro_reinj_rate, dtype=float),
        neutron_reinjection_rate_per_gev=np.asarray(qntr_reinj_rate, dtype=float),
        proton_loss_rate=np.asarray(afpg_pro, dtype=float),
        neutron_loss_rate=np.asarray(afpg_ntr, dtype=float),
        photon_loss_energy_gev=e_gam_t,
        photon_loss_rate=np.asarray(afpg_rad, dtype=float),
    )


@dataclass(frozen=True)
class BetheHeitlerOutput:
    electron_energy_gev: np.ndarray
    pair_rate_per_gev: np.ndarray
    proton_loss_rate: np.ndarray
    photon_loss_rate: np.ndarray


def solve_bethe_heitler(
    proton_energy_gev: np.ndarray,
    proton_density_per_gev: np.ndarray,
    photon_energy_gev: np.ndarray,
    photon_density_per_gev: np.ndarray,
    electron_energy_gev: np.ndarray,
) -> BetheHeitlerOutput:
    e_p = as_strictly_increasing_grid(proton_energy_gev, "proton_energy_gev")
    n_p = as_matching(proton_density_per_gev, e_p, "proton_density_per_gev")
    e_ph = as_strictly_increasing_grid(photon_energy_gev, "photon_energy_gev")
    n_ph = as_matching(photon_density_per_gev, e_ph, "photon_density_per_gev")
    e_e = as_strictly_increasing_grid(electron_energy_gev, "electron_energy_gev")

    pair_rate_per_gev, proton_loss_rate, photon_loss_rate = hadronic_fortran_module.bethe_heitler(
        e_p,
        n_p,
        e_ph,
        n_ph,
        e_e,
    )
    pair_rate_per_gev = np.asarray(pair_rate_per_gev, dtype=float)
    proton_loss_rate = np.asarray(proton_loss_rate, dtype=float)
    photon_loss_rate = np.asarray(photon_loss_rate, dtype=float)

    return BetheHeitlerOutput(
        electron_energy_gev=e_e,
        pair_rate_per_gev=pair_rate_per_gev,
        proton_loss_rate=proton_loss_rate,
        photon_loss_rate=photon_loss_rate,
    )


@dataclass(frozen=True)
class HadronicInverseComptonKernel:
    """AM3-compatible delta-kernel indices for hadronic IC emission."""

    dln_energy: float
    ind_min_energy_pho_hadgrid: int
    delta_e_p: np.ndarray
    jmax_p: np.ndarray
    delta_e_pi: np.ndarray
    jmax_pi: np.ndarray
    delta_e_mu: np.ndarray
    jmax_mu: np.ndarray


@dataclass(frozen=True)
class HadronicInverseComptonOutput:
    """Hadronic IC emissivity output (proton/pion/muon channels)."""

    photon_energy_gev: np.ndarray
    epsilon_p_ic: np.ndarray
    epsilon_pi_ic: np.ndarray
    epsilon_mu_ic: np.ndarray
    coeff_p_cgs: float
    coeff_pi_cgs: float
    coeff_mu_cgs: float
    kernel: HadronicInverseComptonKernel


def initialize_hadronic_inverse_compton_kernel(
    hadron_energy_gev: np.ndarray,
    photon_energy_gev: np.ndarray,
    ind_min_energy_pho_hadgrid: int,
) -> HadronicInverseComptonKernel:
    """Initialize AM3 HadIC small-step kernels.

    This follows AM3 `HadIC::initialize_kernels` exactly:
    gamma = E_had / m,
    deltaE = int(log(gamma^2)/dlnE),
    jmax1 = int(log(0.5*m/gamma)/dlnE) + ind_min_energy_pho_hadgrid,
    jmax = min(nbins_pho, jmax1 + deltaE).
    """

    e_had = as_strictly_increasing_grid(hadron_energy_gev, "hadron_energy_gev", require_finite=True)
    e_ph = as_strictly_increasing_grid(photon_energy_gev, "photon_energy_gev", require_finite=True)
    zeros_ph = np.zeros_like(e_ph)
    zeros_had = np.zeros_like(e_had)
    (
        _eps_p,
        _eps_pi,
        _eps_mu,
        _coeff_p,
        _coeff_pi,
        _coeff_mu,
        dln_energy,
        delta_e_p,
        jmax_p,
        delta_e_pi,
        jmax_pi,
        delta_e_mu,
        jmax_mu,
    ) = hadronic_fortran_module.hadronic_ic(
        e_had,
        e_ph,
        zeros_ph,
        zeros_had,
        zeros_had,
        zeros_had,
        zeros_had,
        zeros_had,
        zeros_had,
        zeros_had,
        int(ind_min_energy_pho_hadgrid),
    )

    return HadronicInverseComptonKernel(
        dln_energy=float(dln_energy),
        ind_min_energy_pho_hadgrid=int(ind_min_energy_pho_hadgrid),
        delta_e_p=np.asarray(delta_e_p, dtype=np.int64),
        jmax_p=np.asarray(jmax_p, dtype=np.int64),
        delta_e_pi=np.asarray(delta_e_pi, dtype=np.int64),
        jmax_pi=np.asarray(jmax_pi, dtype=np.int64),
        delta_e_mu=np.asarray(delta_e_mu, dtype=np.int64),
        jmax_mu=np.asarray(jmax_mu, dtype=np.int64),
    )


def solve_hadronic_inverse_compton(
    hadron_energy_gev: np.ndarray,
    photon_energy_gev: np.ndarray,
    photons_on_had_grid_per_gev: np.ndarray,
    protons_per_gev: np.ndarray,
    pion_plus_per_gev: np.ndarray,
    pion_minus_per_gev: np.ndarray,
    muon_minus_left_per_gev: np.ndarray,
    muon_minus_right_per_gev: np.ndarray,
    muon_plus_left_per_gev: np.ndarray,
    muon_plus_right_per_gev: np.ndarray,
    ind_min_energy_pho_hadgrid: int = 0,
    kernel: HadronicInverseComptonKernel | None = None,
) -> HadronicInverseComptonOutput:
    """Compute hadronic-grid IC emission in AM3 delta approximation.

    Input/Output units follow AM3's HadIC path:
    - energies in GeV
    - distributions as number density per GeV [cm^-3 GeV^-1]
    - output epsilon_* arrays are AM3-style emissivity terms with
      coeff = c * sigma_T * (m_e / m)^2 and 1 dlnE integration factor.
    """

    e_had = as_strictly_increasing_grid(hadron_energy_gev, "hadron_energy_gev", require_finite=True)
    e_ph = as_strictly_increasing_grid(photon_energy_gev, "photon_energy_gev", require_finite=True)
    photons = as_matching(photons_on_had_grid_per_gev, e_ph, "photons_on_had_grid_per_gev", require_finite=True)
    protons = as_matching(protons_per_gev, e_had, "protons_per_gev", require_finite=True)
    pion_plus = as_matching(pion_plus_per_gev, e_had, "pion_plus_per_gev", require_finite=True)
    pion_minus = as_matching(pion_minus_per_gev, e_had, "pion_minus_per_gev", require_finite=True)
    muon_minus_left = as_matching(muon_minus_left_per_gev, e_had, "muon_minus_left_per_gev", require_finite=True)
    muon_minus_right = as_matching(muon_minus_right_per_gev, e_had, "muon_minus_right_per_gev", require_finite=True)
    muon_plus_left = as_matching(muon_plus_left_per_gev, e_had, "muon_plus_left_per_gev", require_finite=True)
    muon_plus_right = as_matching(muon_plus_right_per_gev, e_had, "muon_plus_right_per_gev", require_finite=True)

    (
        epsilon_p,
        epsilon_pi,
        epsilon_mu,
        coeff_p,
        coeff_pi,
        coeff_mu,
        dln_energy,
        delta_e_p,
        jmax_p,
        delta_e_pi,
        jmax_pi,
        delta_e_mu,
        jmax_mu,
    ) = hadronic_fortran_module.hadronic_ic(
        e_had,
        e_ph,
        photons,
        protons,
        pion_plus,
        pion_minus,
        muon_minus_left,
        muon_minus_right,
        muon_plus_left,
        muon_plus_right,
        int(ind_min_energy_pho_hadgrid),
    )
    if kernel is None:
        kernel = HadronicInverseComptonKernel(
            dln_energy=float(dln_energy),
            ind_min_energy_pho_hadgrid=int(ind_min_energy_pho_hadgrid),
            delta_e_p=np.asarray(delta_e_p, dtype=np.int64),
            jmax_p=np.asarray(jmax_p, dtype=np.int64),
            delta_e_pi=np.asarray(delta_e_pi, dtype=np.int64),
            jmax_pi=np.asarray(jmax_pi, dtype=np.int64),
            delta_e_mu=np.asarray(delta_e_mu, dtype=np.int64),
            jmax_mu=np.asarray(jmax_mu, dtype=np.int64),
        )

    coeff_p_target = _am3_ic_coeff(AM3_MASS_PROTON_GEV)
    coeff_pi_target = _am3_ic_coeff(AM3_MASS_PION_CHARGED_GEV)
    coeff_mu_target = _am3_ic_coeff(AM3_MASS_MUON_GEV)
    epsilon_p = np.asarray(epsilon_p, dtype=float) * (coeff_p_target / float(coeff_p))
    epsilon_pi = np.asarray(epsilon_pi, dtype=float) * (coeff_pi_target / float(coeff_pi))
    epsilon_mu = np.asarray(epsilon_mu, dtype=float) * (coeff_mu_target / float(coeff_mu))

    return HadronicInverseComptonOutput(
        photon_energy_gev=e_ph,
        epsilon_p_ic=epsilon_p,
        epsilon_pi_ic=epsilon_pi,
        epsilon_mu_ic=epsilon_mu,
        coeff_p_cgs=float(coeff_p_target),
        coeff_pi_cgs=float(coeff_pi_target),
        coeff_mu_cgs=float(coeff_mu_target),
        kernel=kernel,
    )


def _am3_ic_coeff(mass_gev: float) -> float:
    mass_ratio = mass_gev / AM3_MASS_ELECTRON_GEV
    return AM3_C_CGS * AM3_SIGMA_T_CGS / mass_ratio / mass_ratio


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
    e_ph = as_strictly_increasing_grid(photon_energy_gev, "photon_energy_gev")
    n_ph = as_matching_nonnegative(photon_density_per_gev, e_ph, "photon_density_per_gev")
    e_e = as_strictly_increasing_grid(electron_energy_gev, "electron_energy_gev")
    if int(max_com_energy_factor) < 1:
        raise ValueError("max_com_energy_factor must be >= 1.")

    (
        photon_loss_rate,
        pair_rate_per_gev_per_species,
        pair_rate_per_gev_total,
        absorbed_power,
        injected_power,
    ) = hadronic_fortran_module.pair_production(
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


@dataclass(frozen=True)
class HadronicPPOutput:
    gamma_energy_gev: np.ndarray
    gamma_rate_per_gev: np.ndarray
    neutrino_energy_gev: np.ndarray
    neutrino_rate_per_gev: np.ndarray
    pair_energy_gev: np.ndarray
    pair_rate_per_gev: np.ndarray
    proton_energy_gev: np.ndarray
    proton_loss_rate: np.ndarray


def solve_pp_shell(
    proton_energy_gev: np.ndarray,
    proton_density_per_gev: np.ndarray,
    target_proton_density_cm3: float,
    gamma_energy_gev: np.ndarray,
    neutrino_energy_gev: np.ndarray,
    pair_energy_gev: np.ndarray,
    *,
    kappa_inelastic: float = 0.5,
    pion_energy_fraction: float = 0.17,
    neutral_pion_fraction: float = 1.0 / 3.0,
) -> HadronicPPOutput:
    """Minimal pp microphysics backend in delta-functional approximation.

    Implemented relations:
    1) Inelastic cross section (Kelner et al. 2006, Eq. 79 form):
       sigma_inel(Tp) = (34.3 + 1.88 L + 0.25 L^2) * [1 - (Tp_th/Tp)^4]^2 mb, L = ln(Tp/Tp_th).
    2) Collision frequency:
       t_pp^{-1}(Ep) = n_H * c * sigma_inel(Ep).
    3) Proton loss term (continuous approximation):
       dN_p/dt|loss = -kappa_inelastic * t_pp^{-1} * N_p.
    4) Secondary source terms (delta mapping E_s = x_s * E_p):
       Q_s(E_s) = (m_s / x_s) * t_pp^{-1}(E_p) * N_p(E_p), E_p = E_s / x_s.

    Validity:
    - Physical threshold enforced by sigma_inel(Tp < Tp_th) = 0.
    - Delta approximation is intended as a minimal backend for broad hadronic sweeps;
      accuracy degrades near threshold and for detailed spectral features.
    """

    e_p = as_strictly_increasing_grid(proton_energy_gev, "proton_energy_gev")
    n_p = as_matching(proton_density_per_gev, e_p, "proton_density_per_gev")
    e_gamma = as_strictly_increasing_grid(gamma_energy_gev, "gamma_energy_gev")
    e_nu = as_strictly_increasing_grid(neutrino_energy_gev, "neutrino_energy_gev")
    e_pair = as_strictly_increasing_grid(pair_energy_gev, "pair_energy_gev")
    if target_proton_density_cm3 < 0.0:
        raise ValueError("target_proton_density_cm3 must be non-negative.")
    if not (0.0 < kappa_inelastic <= 1.0):
        raise ValueError("kappa_inelastic must be in (0, 1].")
    if not (0.0 < pion_energy_fraction < 1.0):
        raise ValueError("pion_energy_fraction must be in (0, 1).")
    if not (0.0 <= neutral_pion_fraction <= 1.0):
        raise ValueError("neutral_pion_fraction must be in [0, 1].")

    gamma_rate, neutrino_rate, pair_rate, proton_loss_rate = hadronic_fortran_module.pp_shell(
        e_p,
        n_p,
        float(target_proton_density_cm3),
        e_gamma,
        e_nu,
        e_pair,
        float(kappa_inelastic),
        float(pion_energy_fraction),
        float(neutral_pion_fraction),
    )
    gamma_rate = np.asarray(gamma_rate, dtype=float)
    neutrino_rate = np.asarray(neutrino_rate, dtype=float)
    pair_rate = np.asarray(pair_rate, dtype=float)
    proton_loss_rate = np.asarray(proton_loss_rate, dtype=float)

    return HadronicPPOutput(
        gamma_energy_gev=e_gamma,
        gamma_rate_per_gev=gamma_rate,
        neutrino_energy_gev=e_nu,
        neutrino_rate_per_gev=neutrino_rate,
        pair_energy_gev=e_pair,
        pair_rate_per_gev=pair_rate,
        proton_energy_gev=e_p,
        proton_loss_rate=proton_loss_rate,
    )
