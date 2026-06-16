from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from src import constants
from asgard_core.hadronic_validation import as_matching_nonnegative, as_strictly_increasing_grid

import src.Hadronic.hadronic_forward_1d as hadronic_fortran_module

NEUTRON_MASS_GEV = constants.para_m_n_gev
PI_PLUS_MASS_GEV = constants.para_m_pi_charged_gev
MUON_MASS_GEV = constants.para_m_mu_gev

SPECIES_TRANSPORT_BACKEND = "fortran_species_transport"


@dataclass(frozen=True)
class NeutronDistribution:
    gamma: np.ndarray
    density_per_gamma: np.ndarray


@dataclass(frozen=True)
class ChargedPionDistribution:
    gamma: np.ndarray
    plus_density_per_gamma: np.ndarray
    minus_density_per_gamma: np.ndarray


@dataclass(frozen=True)
class ChargedMuonDistribution:
    gamma: np.ndarray
    minus_left_density_per_gamma: np.ndarray
    minus_right_density_per_gamma: np.ndarray
    plus_left_density_per_gamma: np.ndarray
    plus_right_density_per_gamma: np.ndarray


@dataclass(frozen=True)
class HadronicSpeciesState:
    neutron: NeutronDistribution
    charged_pion: ChargedPionDistribution
    charged_muon: ChargedMuonDistribution


@dataclass(frozen=True)
class HadronicSpeciesSources:
    neutron_per_gamma_s: np.ndarray
    charged_pion_plus_per_gamma_s: np.ndarray
    charged_pion_minus_per_gamma_s: np.ndarray
    charged_muon_minus_left_per_gamma_s: np.ndarray
    charged_muon_minus_right_per_gamma_s: np.ndarray
    charged_muon_plus_left_per_gamma_s: np.ndarray
    charged_muon_plus_right_per_gamma_s: np.ndarray


def spherical_divergence_rate(radius_cm: float, expansion_speed_cm_s: float) -> float:
    """Comoving expansion divergence for spherical flow: div(v)=3 v_exp/R."""
    if radius_cm <= 0.0:
        raise ValueError("radius_cm must be positive.")
    if expansion_speed_cm_s < 0.0:
        raise ValueError("expansion_speed_cm_s must be non-negative.")
    return 3.0 * expansion_speed_cm_s / radius_cm


def advance_species_transport_explicit(
    state: HadronicSpeciesState,
    sources: HadronicSpeciesSources,
    dt_s: float,
    b_field_g: float,
    divergence_rate_s_inv: float,
) -> HadronicSpeciesState:
    gamma = _shared_gamma_grid(state)
    _validate_sources(gamma, sources)
    (
        neutron_next,
        pion_plus_next,
        pion_minus_next,
        muon_minus_left_next,
        muon_minus_right_next,
        muon_plus_left_next,
        muon_plus_right_next,
    ) = hadronic_fortran_module.fs_hadronic_species_transport_shell(
        gamma,
        float(dt_s),
        float(b_field_g),
        float(divergence_rate_s_inv),
        np.asarray(state.neutron.density_per_gamma, dtype=float),
        np.asarray(state.charged_pion.plus_density_per_gamma, dtype=float),
        np.asarray(state.charged_pion.minus_density_per_gamma, dtype=float),
        np.asarray(state.charged_muon.minus_left_density_per_gamma, dtype=float),
        np.asarray(state.charged_muon.minus_right_density_per_gamma, dtype=float),
        np.asarray(state.charged_muon.plus_left_density_per_gamma, dtype=float),
        np.asarray(state.charged_muon.plus_right_density_per_gamma, dtype=float),
        np.asarray(sources.neutron_per_gamma_s, dtype=float),
        np.asarray(sources.charged_pion_plus_per_gamma_s, dtype=float),
        np.asarray(sources.charged_pion_minus_per_gamma_s, dtype=float),
        np.asarray(sources.charged_muon_minus_left_per_gamma_s, dtype=float),
        np.asarray(sources.charged_muon_minus_right_per_gamma_s, dtype=float),
        np.asarray(sources.charged_muon_plus_left_per_gamma_s, dtype=float),
        np.asarray(sources.charged_muon_plus_right_per_gamma_s, dtype=float),
    )

    return HadronicSpeciesState(
        neutron=NeutronDistribution(gamma=gamma, density_per_gamma=neutron_next),
        charged_pion=ChargedPionDistribution(
            gamma=gamma,
            plus_density_per_gamma=pion_plus_next,
            minus_density_per_gamma=pion_minus_next,
        ),
        charged_muon=ChargedMuonDistribution(
            gamma=gamma,
            minus_left_density_per_gamma=muon_minus_left_next,
            minus_right_density_per_gamma=muon_minus_right_next,
            plus_left_density_per_gamma=muon_plus_left_next,
            plus_right_density_per_gamma=muon_plus_right_next,
        ),
    )


def _shared_gamma_grid(state: HadronicSpeciesState) -> np.ndarray:
    gamma_n = as_strictly_increasing_grid(state.neutron.gamma, "neutron.gamma")
    gamma_pi = as_strictly_increasing_grid(state.charged_pion.gamma, "charged_pion.gamma")
    gamma_mu = as_strictly_increasing_grid(state.charged_muon.gamma, "charged_muon.gamma")
    if not (gamma_n.shape == gamma_pi.shape == gamma_mu.shape):
        raise ValueError("All species gamma grids must share the same shape.")
    if not (np.allclose(gamma_n, gamma_pi) and np.allclose(gamma_n, gamma_mu)):
        raise ValueError("All species must share an identical gamma grid.")
    as_matching_nonnegative(state.neutron.density_per_gamma, gamma_n, "neutron.density_per_gamma")
    as_matching_nonnegative(state.charged_pion.plus_density_per_gamma, gamma_n, "charged_pion.plus_density_per_gamma")
    as_matching_nonnegative(state.charged_pion.minus_density_per_gamma, gamma_n, "charged_pion.minus_density_per_gamma")
    as_matching_nonnegative(state.charged_muon.minus_left_density_per_gamma, gamma_n, "charged_muon.minus_left_density_per_gamma")
    as_matching_nonnegative(state.charged_muon.minus_right_density_per_gamma, gamma_n, "charged_muon.minus_right_density_per_gamma")
    as_matching_nonnegative(state.charged_muon.plus_left_density_per_gamma, gamma_n, "charged_muon.plus_left_density_per_gamma")
    as_matching_nonnegative(state.charged_muon.plus_right_density_per_gamma, gamma_n, "charged_muon.plus_right_density_per_gamma")
    return gamma_n


def _validate_sources(gamma: np.ndarray, sources: HadronicSpeciesSources) -> None:
    as_matching_nonnegative(sources.neutron_per_gamma_s, gamma, "sources.neutron_per_gamma_s")
    as_matching_nonnegative(sources.charged_pion_plus_per_gamma_s, gamma, "sources.charged_pion_plus_per_gamma_s")
    as_matching_nonnegative(sources.charged_pion_minus_per_gamma_s, gamma, "sources.charged_pion_minus_per_gamma_s")
    as_matching_nonnegative(sources.charged_muon_minus_left_per_gamma_s, gamma, "sources.charged_muon_minus_left_per_gamma_s")
    as_matching_nonnegative(sources.charged_muon_minus_right_per_gamma_s, gamma, "sources.charged_muon_minus_right_per_gamma_s")
    as_matching_nonnegative(sources.charged_muon_plus_left_per_gamma_s, gamma, "sources.charged_muon_plus_left_per_gamma_s")
    as_matching_nonnegative(sources.charged_muon_plus_right_per_gamma_s, gamma, "sources.charged_muon_plus_right_per_gamma_s")
