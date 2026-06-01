from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from asgard_core.asgard_config import FitResult
from asgard_core.asgard_observables import LIGHT_CURVE_PLOT_SPECS


def _normalize_spectrum_quantity(quantity: str) -> str:
    kind = quantity.lower()
    if kind == "sed":
        return "nufnu"
    if kind not in ("nufnu", "fnu"):
        raise ValueError("quantity must be 'sed', 'nufnu', or 'fnu'.")
    return kind


def _to_2d_matrix(values: np.ndarray, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 2:
        raise ValueError(f"{name} must be a 2D array, got shape {arr.shape}.")
    return arr


def plot_flux_matrix(
    times_s: np.ndarray,
    freq_hz: np.ndarray,
    flux_matrix: np.ndarray,
    labels: list[str] | None = None,
    scales: list[float] | None = None,
    outfile: str = "Flux_Matrix.pdf",
    show: bool = False,
    xlabel: str = "Time [s]",
    ylabel: str = "Flux [erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$]",
    title: str = "Flux Curves",
    linewidth: float = 2.0,
):
    times = np.asarray(times_s, dtype=float)
    freq = np.asarray(freq_hz, dtype=float)
    matrix = _to_2d_matrix(flux_matrix, "flux_matrix")
    if matrix.shape[0] != freq.size:
        raise ValueError(f"flux_matrix shape {matrix.shape} does not match len(freq_hz)={freq.size}.")
    if matrix.shape[1] != times.size:
        raise ValueError(f"flux_matrix shape {matrix.shape} does not match len(times_s)={times.size}.")

    n = matrix.shape[0]
    if labels is None:
        labels = [f"{f:.3g} Hz" for f in freq]
    if scales is None:
        scales = [1.0] * n
    if len(scales) != n:
        raise ValueError(f"scales length {len(scales)} does not match n_curves {n}.")

    fig, ax = plt.subplots(figsize=(10, 8))
    for i in range(n):
        ax.loglog(times, matrix[i] * scales[i], color=f"C{i}", linewidth=linewidth, label=labels[i])

    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_title(title, fontsize=18)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.tick_params(which="both", direction="in", bottom=True, top=True, left=True, right=True, labelsize=14)
    if n > 0:
        ax.legend(bbox_to_anchor=(1, 1), loc="upper right", ncol=1, borderaxespad=0.5, fancybox=True, fontsize=10)

    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    if show:
        plt.show()
    return fig


def plot_sed_grid(
    times_s: np.ndarray,
    freqs_hz: np.ndarray,
    sed_matrix: np.ndarray,
    labels: list[str] | None = None,
    outfile: str = "SED.pdf",
    show: bool = False,
    quantity: str = "nufnu",
    xlabel: str = "Frequency [Hz]",
    ylabel: str | None = None,
    title: str = "Observed Spectra",
    linewidth: float = 2.0,
):
    quantity = _normalize_spectrum_quantity(quantity)

    times = np.asarray(times_s, dtype=float)
    freqs = np.asarray(freqs_hz, dtype=float)
    sed = _to_2d_matrix(sed_matrix, "sed_matrix")
    if sed.shape[0] != freqs.size or sed.shape[1] != times.size:
        raise ValueError(
            f"sed_matrix shape {sed.shape} should be ({freqs.size}, {times.size})."
        )

    if labels is None:
        labels = [f"t={times[i]:.3e} s" for i in range(len(times))]

    if ylabel is None:
        ylabel = r"$\nu F_\nu$ [erg cm$^{-2}$ s$^{-1}$]" if quantity == "nufnu" else r"$F_\nu$ [erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$]"

    fig, ax = plt.subplots(figsize=(9, 7))
    for i in range(times.size):
        spec = sed[:, i] * freqs if quantity == "nufnu" else sed[:, i]
        ax.loglog(freqs, spec, linewidth=linewidth, label=labels[i])

    ax.set_xlabel(xlabel, fontsize=13)
    ax.set_ylabel(ylabel, fontsize=13)
    ax.set_title(title, fontsize=16)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.tick_params(which="both", direction="in", top=True, right=True, labelsize=12)
    ax.legend(loc="best", fontsize=10, ncol=2)

    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    if show:
        plt.show()
    return fig


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
    quantity = _normalize_spectrum_quantity(quantity)

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
