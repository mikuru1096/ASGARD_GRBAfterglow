from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from asgard_models import FitResult
from asgard_observables import LIGHT_CURVE_PLOT_SPECS

FORWARD_CHARACTERISTIC_SPECS = (
    ("nu_m", "#FF4500", r"FS $\nu_m$", "-"),
    ("nu_c", "#1E90FF", r"FS $\nu_c$", "-"),
    ("nu_a", "#228B22", r"FS $\nu_a$", "-"),
)
REVERSE_CHARACTERISTIC_SPECS = (
    ("rs_nu_m", "#FF8C00", r"RS $\nu_m$", "--"),
    ("rs_nu_c", "#00BFFF", r"RS $\nu_c$", "--"),
    ("rs_nu_a", "#32CD32", r"RS $\nu_a$", "--"),
)


def plot_light_curve(
    result: FitResult,
    outfile: str = "Radiation_Lightcurves.pdf",
    show: bool = False,
) -> plt.Figure:
    t_obs = result.t_obs_s
    flux = result.bands_flux

    fig, axes = plt.subplots(figsize=(10, 8))
    for index, (_, scale, color, label) in enumerate(LIGHT_CURVE_PLOT_SPECS):
        axes.loglog(t_obs, flux[index] * scale, color=color, linewidth=2.5, label=label)

    axes.set_xlim([1.0e2, 1.0e8])
    axes.set_ylim([1.0e-22, 1.0e-4])
    axes.set_xlabel("Time [s]", fontsize=14)
    axes.set_ylabel("Flux erg/cm$^2$/s or $\\mu$Jy", fontsize=14)
    axes.set_title("LCs", fontsize=18)
    axes.tick_params(which="both", direction="in", bottom=True, top=True, left=True, right=True, labelsize=14)
    axes.legend(bbox_to_anchor=(1, 1), loc="upper right", ncol=1, borderaxespad=0.5, fancybox=True, fontsize=10)

    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    if show:
        plt.show()
    return fig


def plot_spectrum(
    result: FitResult,
    times_s: list[float],
    quantity: str = "nufnu",
    outfile: str = "Radiation_Spectra.pdf",
    show: bool = False,
) -> plt.Figure:
    if result.spectrum_fnu is None or result.spectrum_freq_hz is None:
        raise ValueError("FitResult does not contain spectrum data. Enable spectrum_output in FitConfig first.")
    if quantity not in ("nufnu", "fnu"):
        raise ValueError("quantity must be 'nufnu' or 'fnu'.")

    frequencies = result.spectrum_freq_hz
    spectrum = result.spectrum_fnu
    t_obs = result.t_obs_s

    fig, ax = plt.subplots(figsize=(9, 7))
    for t_req in times_s:
        idx = int(np.abs(t_obs - t_req).argmin())
        y = spectrum[:, idx] * frequencies if quantity == "nufnu" else spectrum[:, idx]
        ax.loglog(frequencies, y, linewidth=2.0, label=f"t={t_obs[idx]:.3e} s")

    ax.set_xlabel("Frequency [Hz]", fontsize=13)
    ax.set_ylabel(
        r"$\nu F_\nu$ [erg cm$^{-2}$ s$^{-1}$]"
        if quantity == "nufnu"
        else r"$F_\nu$ [erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$]",
        fontsize=13,
    )
    ax.set_title("Observed Spectra", fontsize=16)
    ax.tick_params(which="both", direction="in", top=True, right=True, labelsize=12)
    ax.legend(loc="best", fontsize=10)

    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    if show:
        plt.show()
    return fig


def plot_characteristic_frequencies(
    result: FitResult,
    include_reverse: bool = True,
    outfile: str = "Characteristic_Frequencies.pdf",
    show: bool = False,
) -> plt.Figure:
    t_obs = result.characteristic_time_s

    fig, ax = plt.subplots(figsize=(9, 7))
    _plot_characteristic_group(ax, t_obs, result, FORWARD_CHARACTERISTIC_SPECS)

    if include_reverse:
        if result.rs_nu_m is None or result.rs_nu_c is None or result.rs_nu_a is None:
            raise ValueError("FitResult does not contain reverse-shock characteristic frequencies.")
        _plot_characteristic_group(ax, t_obs, result, REVERSE_CHARACTERISTIC_SPECS)

    ax.set_xlabel("Time [s]", fontsize=13)
    ax.set_ylabel("Frequency [Hz]", fontsize=13)
    ax.set_title("Characteristic Frequencies", fontsize=16)
    ax.tick_params(which="both", direction="in", top=True, right=True, labelsize=12)
    ax.legend(loc="best", fontsize=10)

    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    if show:
        plt.show()
    return fig


def _plot_characteristic_group(
    ax: plt.Axes,
    t_obs_s: np.ndarray,
    result: FitResult,
    specs: tuple[tuple[str, str, str, str], ...],
) -> None:
    for attribute, color, label, linestyle in specs:
        _plot_frequency_series(ax, t_obs_s, getattr(result, attribute), color, label, linestyle)


def _plot_frequency_series(
    ax: plt.Axes,
    t_obs_s: np.ndarray,
    frequency_hz: np.ndarray,
    color: str,
    label: str,
    linestyle: str,
) -> None:
    if t_obs_s.shape != frequency_hz.shape:
        raise ValueError(
            f"Characteristic-frequency time axis shape {t_obs_s.shape} does not match frequency shape {frequency_hz.shape}."
        )

    mask = frequency_hz > 0.0
    ax.loglog(t_obs_s[mask], frequency_hz[mask], color=color, linestyle=linestyle, linewidth=2.0, label=label)
