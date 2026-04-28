"""电磁对级联：Fortran 核迭代驱动。"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from src import constants

try:
    import src.Hadronic.FS_hadronic_1d as hadronic_fortran_module
except ImportError:
    hadronic_fortran_module = None

_HAS_CASCADE = hadronic_fortran_module is not None and hasattr(
    hadronic_fortran_module, "fs_hadronic_pair_cascade_step"
)


@dataclass(frozen=True)
class PairCascadeOutput:
    photon_energy_gev: np.ndarray
    electron_energy_gev: np.ndarray
    frequency_hz: np.ndarray
    photon_loss_rate_local_s_inv: np.ndarray
    tau_pair_path: np.ndarray
    shell_path_time_s: float
    pair_injection_rate_per_gev_total_local: np.ndarray
    pair_syn_luminosity_hz: np.ndarray
    pair_ic_luminosity_hz: np.ndarray
    pair_total_luminosity_hz: np.ndarray
    absorbed_power_gev_per_cm3_s: float
    injected_power_gev_per_cm3_s: float
    cascade_iterations: int = 1
    cascade_converged: bool = True


def shell_path_time_seconds(radius_cm: float, gamma_bulk: float) -> float:
    return float(radius_cm) / (12.0 * max(float(gamma_bulk), 1.0) * constants.para_c)


def compute_iterative_pair_cascade(
    photon_energy_gev: np.ndarray,
    photon_density_per_gev: np.ndarray,
    electron_energy_gev: np.ndarray,
    radius_cm: float,
    gamma_bulk: float,
    b_field_g: float,
    max_iterations: int = 5,
    convergence_rtol: float = 1e-3,
) -> PairCascadeOutput:
    """迭代电磁对级联：Fortran 核 + Python 收敛控制。"""
    if not _HAS_CASCADE:
        raise RuntimeError("Pair cascade Fortran kernel not available.")

    e_ph = np.asarray(photon_energy_gev, dtype=float)
    n_ph_primary = np.asarray(photon_density_per_gev, dtype=float)
    n_ph_primary = np.maximum(n_ph_primary, 0.0)
    e_e = np.asarray(electron_energy_gev, dtype=float)
    n_ph = n_ph_primary.copy()
    path_time = shell_path_time_seconds(radius_cm, gamma_bulk)

    cascade_accum = np.zeros_like(n_ph)
    total_absorbed = 0.0
    freq_hz = e_ph * constants.para_gev2hz
    converged = False
    final_iter = 1

    for iteration in range(max_iterations):
        cascade_syn, absorbed = hadronic_fortran_module.fs_hadronic_pair_cascade_step(
            e_ph, n_ph, e_e,
            max(float(b_field_g), 1e-30),
            path_time,
        )
        cascade_syn = np.maximum(np.asarray(cascade_syn, dtype=float), 0.0)
        total_absorbed += float(absorbed)

        cascade_new = cascade_syn / max(e_ph * constants.para_h_gev, 1e-60)
        rel_change = 0.0
        if iteration > 0 and np.max(cascade_accum) > 0.0:
            rel_change = float(np.max(np.abs(cascade_new - cascade_accum)) / np.max(cascade_accum))
        cascade_accum = cascade_new

        if rel_change < convergence_rtol and iteration > 0:
            converged = True
            final_iter = iteration + 1
            break

        n_ph = n_ph_primary + cascade_accum
        final_iter = iteration + 1

    # 构建向后兼容的输出
    photon_loss = np.zeros_like(n_ph)
    tau_pair = np.zeros_like(n_ph)

    return PairCascadeOutput(
        photon_energy_gev=e_ph,
        electron_energy_gev=e_e,
        frequency_hz=freq_hz,
        photon_loss_rate_local_s_inv=photon_loss,
        tau_pair_path=tau_pair,
        shell_path_time_s=path_time,
        pair_injection_rate_per_gev_total_local=np.zeros_like(e_e),
        pair_syn_luminosity_hz=cascade_accum,
        pair_ic_luminosity_hz=np.zeros_like(cascade_accum),
        pair_total_luminosity_hz=cascade_accum,
        absorbed_power_gev_per_cm3_s=total_absorbed,
        injected_power_gev_per_cm3_s=0.0,
        cascade_iterations=final_iter,
        cascade_converged=converged,
    )
