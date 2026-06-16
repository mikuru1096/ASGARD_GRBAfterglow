"""反向激波强子过程：Fortran 核心驱动 + 最简 Python 包装。"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from src import constants

import src.Hadronic.hadronic_reverse_1d as _hadronic_reverse_module


@dataclass(frozen=True)
class ReverseShockHadronicSolution:
    gam_p: np.ndarray
    d_n_gam_p: np.ndarray
    l_had_syn_spec: np.ndarray
    seed_had_syn: np.ndarray


def solve_rs_hadronic_core(
    r_tobs_s: np.ndarray,
    r_gamma: np.ndarray,
    radius_cm: np.ndarray,
    rs_swept_mass_g: np.ndarray,
    rs_b_field_g: np.ndarray,
    v_seed_hz: np.ndarray,
    num_gam_p: int,
    epsilon_p: float,
    include_proton_synch: bool = True,
) -> ReverseShockHadronicSolution:
    tobs = np.asarray(r_tobs_s, dtype=float)
    gamma = np.asarray(r_gamma, dtype=float)
    radius = np.asarray(radius_cm, dtype=float)
    swept = np.asarray(rs_swept_mass_g, dtype=float)
    b_field = np.asarray(rs_b_field_g, dtype=float)
    v_seed = np.asarray(v_seed_hz, dtype=float)
    num_nu = int(v_seed.size)
    num_r = int(radius.size)
    np_gam_p = int(num_gam_p)

    shell_energy = float(epsilon_p) * swept * np.maximum(gamma - 1.0, 0.0) * constants.para_c * constants.para_c

    gam_p, dN_gam_p, P_had_syn, Seed_had_syn = _hadronic_reverse_module.fs_hadronic_reverse_1d(
        tobs, gamma, radius, shell_energy, b_field, v_seed,
        1 if include_proton_synch else 0,
        np_gam_p,
    )
    return ReverseShockHadronicSolution(
        gam_p=np.asarray(gam_p, dtype=float),
        d_n_gam_p=np.asarray(dN_gam_p, dtype=float).reshape(np_gam_p, num_r),
        l_had_syn_spec=np.asarray(P_had_syn, dtype=float).reshape(num_nu, num_r),
        seed_had_syn=np.asarray(Seed_had_syn, dtype=float).reshape(num_nu, num_r),
    )
