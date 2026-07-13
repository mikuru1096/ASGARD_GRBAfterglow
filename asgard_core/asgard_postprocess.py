from __future__ import annotations

from pathlib import Path

import numpy as np
from extinction import fitzpatrick99 as f99
from scipy import interpolate

from asgard_core.asgard_config import RuntimeConfig, SimulationSetup
from src import Interpolation, constants

package_root = Path(__file__).resolve().parent
lightdir = package_root / "data_light_curve"
fitbands = ("xrt", "optr", "optz", "opti", "optg", "9GHz", "5.5GHz", "3GHz")
pointbands = fitbands[1:]
fitfreqs = np.array([2.7e17, 4.63e14, 3.39e14, 4.01e14, 6.42e14, 9e9, 5.5e9, 3e9], dtype=float)
xrtbins = 8
obsfreqs = np.concatenate(
    (
        np.logspace(np.log10(0.5 * constants.para_kev2hz), np.log10(10.0 * constants.para_kev2hz), xrtbins),
        fitfreqs[1:],
    ),
    axis=0,
)
zeropoint = np.array([0.0, -48.6, -48.6, -48.6, -48.6, 0.0, 0.0, 0.0], dtype=float)


def band_freqs() -> tuple[int, np.ndarray]:
    return xrtbins, obsfreqs.copy()


def combine_flux(flux: np.ndarray, frequency: np.ndarray, numxrt: int) -> np.ndarray:
    nband = len(pointbands)
    xrtflux = np.trapezoid(flux[:numxrt].T, frequency[:numxrt], axis=1).reshape(1, -1)
    poinflux = flux[numxrt : numxrt + nband] * 1e29
    return np.vstack([xrtflux, poinflux])


def observe_flux(
    setup: SimulationSetup,
    chartime: np.ndarray,
    gamma: np.ndarray,
    radius_cm: np.ndarray,
    sourceflux: np.ndarray,
    frequency: np.ndarray,
    config: RuntimeConfig,
) -> np.ndarray:
    sourcebatch = np.asfortranarray(np.asarray(sourceflux, dtype=float)[:, :, None])
    return observe_flux_batch(setup, chartime, gamma, radius_cm, sourcebatch, frequency, config)[:, :, 0]


def observe_flux_batch(
    setup: SimulationSetup,
    chartime: np.ndarray,
    gamma: np.ndarray,
    radius_cm: np.ndarray,
    sourceflux: np.ndarray,
    frequency: np.ndarray,
    config: RuntimeConfig,
) -> np.ndarray:
    frequency = np.asarray(frequency, dtype=float)
    order = np.argsort(frequency)
    sortedfreq = frequency[order]
    geometry_kernel = str(config.geometry_kernel).lower()
    if geometry_kernel == "sed_adaptive_theta":
        kernel = Interpolation.sed_adaptive_theta_batch
        args = (
            float(config.projection_adaptive_rtol),
            int(config.projection_adaptive_max_depth),
        )
    else:
        kernel = Interpolation.sed_interpolation_batch
        args = ()
    sortedflux = kernel(
        setup.boundary,
        chartime,
        gamma,
        radius_cm,
        np.asfortranarray(sourceflux),
        setup.seed_frequency_hz,
        sortedfreq,
        setup.observer_time_s,
        *args,
        config.eats_num_theta,
        config.eats_num_phi,
        config.num_threads,
    )
    if np.array_equal(order, np.arange(order.shape[0])):
        return sortedflux

    flux = np.empty_like(sortedflux)
    flux[order, :, :] = sortedflux
    return flux


def observe_pairbatch(
    setup: SimulationSetup,
    chartime: np.ndarray,
    gamma: np.ndarray,
    radius_cm: np.ndarray,
    sourceflux: np.ndarray,
    pairfreq: np.ndarray,
    pairstart: np.ndarray,
    config: RuntimeConfig,
) -> np.ndarray:
    return Interpolation.sed_pairbatch(
        setup.boundary,
        chartime,
        gamma,
        radius_cm,
        np.asfortranarray(sourceflux),
        setup.seed_frequency_hz,
        np.asarray(pairfreq, dtype=float),
        setup.observer_time_s,
        np.asarray(pairstart, dtype=np.int32),
        config.eats_num_theta,
        config.eats_num_phi,
        config.num_threads,
    )


def light_chi(
    bandflux: np.ndarray,
    obstime: np.ndarray,
    config: RuntimeConfig,
) -> float:
    chi2 = _chi2sum(
        list(fitbands),
        bandflux[: len(fitbands)],
        obstime,
        fitfreqs,
        config.rv,
        config.ebv,
        zeropoint,
        config.f_sys,
    )
    return float(chi2)


def _chi2sum(
    bands,
    curves,
    serial,
    frequency,
    rv,
    ebv,
    zeropt,
    fsys,
):
    if not lightdir.exists():
        raise FileNotFoundError(f"Data directory does not exist: {lightdir}")

    chi2 = 0.0
    interp = _modelinterp(curves, serial)

    for data_file in lightdir.glob("*"):
        if not data_file.is_file():
            continue
        name = data_file.stem
        if name not in bands:
            continue
        table = _loadtable(data_file)
        times, obsflux, obserr = _parsedata(table)
        band = bands.index(name)
        if name.startswith("opt"):
            obsflux, obserr = _optext(
                obsflux,
                obserr,
                frequency[band],
                rv,
                ebv,
                zeropt[band],
            )
        times = times * 86400
        _checkrange(times, serial)
        modflux = interp[band](times)
        sigma = np.where(modflux > obsflux, obserr[0], obserr[1]) if obserr.ndim == 2 else obserr
        variance = sigma**2 if fsys <= 0 else sigma**2 + (obsflux * fsys) ** 2
        chi2 += np.sum((modflux - obsflux) ** 2 / variance)
    return chi2


def _loadtable(data_file: Path) -> np.ndarray:
    table = np.loadtxt(data_file)
    if table.ndim == 1:
        table = table.reshape(1, -1)
    return table


def _modelinterp(curves: np.ndarray, serial: np.ndarray) -> list[interpolate.interp1d]:
    return [
        interpolate.interp1d(
            serial,
            curves[index, :],
            kind="linear",
            bounds_error=True,
        )
        for index in range(len(curves))
    ]


def _checkrange(samples: np.ndarray, serial: np.ndarray) -> None:
    if np.min(samples) < serial[0] or np.max(samples) > serial[-1]:
        raise ValueError("The model curve cannot fully cover the data range.")


def _parsedata(table):
    ncols = table.shape[1] if table.ndim > 1 else table.size
    if ncols == 6:
        times, flux = table[:, 0], table[:, 3]
        error = np.vstack((table[:, 4], np.abs(table[:, 5])))
    elif ncols == 4:
        times = 0.5 * (table[:, 0] + table[:, 1]) if table[0, 0] < table[0, 1] else table[:, 0]
        flux, error = table[:, 2], table[:, 3]
    elif ncols == 3:
        times, flux, error = table[:, 0], table[:, 1], table[:, 2]
    else:
        raise ValueError(f"The observation data should be 3, 4, or 6 columns. Currently, there are {ncols} columns.")
    return times, flux, error


def _magflux(magnitude, magerr, zeropt):
    flux = 1e29 * 10 ** (0.4 * (zeropt - magnitude))
    error = 0.4 * np.log(10.0) * flux * magerr
    return flux, error


def _optext(
    magnitude,
    magerr,
    frequency,
    rv,
    ebv,
    zeropt,
):
    av = rv * ebv
    wave = np.array([2.997e10 / frequency * 1e8], dtype=float)
    magnitude = magnitude - f99(wave, av, rv)
    return _magflux(magnitude, magerr, zeropt)
