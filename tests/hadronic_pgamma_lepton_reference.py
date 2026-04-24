from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_core.hadronic_pgamma import (
    ETA0,
    kelner_aharonian_2008_gamma_spectrum,
    kelner_aharonian_2008_secondary_phi,
    kelner_aharonian_2008_secondary_spectrum,
)


CHANNEL_ORDER = ("e_plus", "nu_mu_bar", "nu_e", "nu_mu", "e_minus", "nu_e_bar")
CHANNEL_ETA_RATIO = {
    "e_plus": 2.0,
    "nu_mu_bar": 2.0,
    "nu_e": 2.0,
    "nu_mu": 5.0,
    "e_minus": 4.0,
    "nu_e_bar": 4.0,
}
PHI_X_GRID = np.array([1.0e-5, 3.0e-5, 1.0e-4, 3.0e-4, 1.0e-3, 3.0e-3, 1.0e-2, 3.0e-2, 1.0e-1, 3.0e-1], dtype=float)
PHI_REFERENCE = {
    "e_plus": np.array(
        [
            4.9049859221896135e-17,
            4.9049859221896135e-17,
            4.9049859221896135e-17,
            4.9049859221896135e-17,
            4.9049859221896135e-17,
            4.9049859221896135e-17,
            4.9019001568822240e-17,
            3.6260479554835298e-17,
            7.0220256372622712e-18,
            6.8160775072070075e-21,
        ],
        dtype=float,
    ),
    "nu_mu_bar": np.array(
        [
            4.8769574312056725e-17,
            4.8769574312056725e-17,
            4.8769574312056725e-17,
            4.8769574312056725e-17,
            4.8769574312056725e-17,
            4.8769574312056725e-17,
            4.8740647583937517e-17,
            3.6186816063754010e-17,
            6.9726916033999495e-18,
            6.5479935721349962e-21,
        ],
        dtype=float,
    ),
    "nu_e": np.array(
        [
            5.7458406517078330e-17,
            5.7458406517078330e-17,
            5.7458406517078330e-17,
            5.7458406517078330e-17,
            5.7458406517078330e-17,
            5.7458406517078330e-17,
            5.7389126525523122e-17,
            3.7961986305749784e-17,
            5.3687992358259347e-18,
            3.5061119857865257e-21,
        ],
        dtype=float,
    ),
    "nu_mu": np.array(
        [
            1.2715871208356813e-16,
            1.2715871208356813e-16,
            1.2715871208356813e-16,
            1.2715871208356813e-16,
            1.2715871208356813e-16,
            1.2715871208356813e-16,
            1.2061369025498529e-16,
            6.8264284678171995e-17,
            1.1407743315447103e-17,
            9.6996941667708782e-21,
        ],
        dtype=float,
    ),
    "e_minus": np.array(
        [
            9.9100000000000000e-18,
            9.9100000000000000e-18,
            9.9100000000000000e-18,
            9.9100000000000000e-18,
            9.9100000000000000e-18,
            9.9100000000000000e-18,
            9.9100000000000000e-18,
            7.4032950411115620e-18,
            5.1190062269103780e-19,
            7.5655604227106810e-22,
        ],
        dtype=float,
    ),
    "nu_e_bar": np.array(
        [
            9.7400000000000006e-18,
            9.7400000000000006e-18,
            9.7400000000000006e-18,
            9.7400000000000006e-18,
            9.7400000000000006e-18,
            9.7400000000000006e-18,
            9.7400000000000006e-18,
            7.1195997496213267e-18,
            2.4966026500385035e-19,
            3.3621054793488804e-23,
        ],
        dtype=float,
    ),
}

SPECTRUM_ENERGY_GRID = np.array(
    [
        3.1622776601683793e01,
        5.2853885930792470e01,
        8.8339278146640580e01,
        1.4764908816142600e02,
        2.4677871148904760e02,
        4.1246263829013566e02,
        6.8938453790739360e02,
        1.1522280977398216e03,
        1.9258186341850783e03,
        3.2187875118212390e03,
        5.3798384034436920e03,
        8.9917899646601330e03,
        1.5028757502606059e04,
        2.5118864315095820e04,
    ],
    dtype=float,
)
SPECTRUM_REFERENCE = {
    "e_plus": np.array(
        [
            7.1982928991404430e-54,
            7.1763947167655360e-54,
            7.0554518878234260e-54,
            6.7106980577369710e-54,
            6.0357867765879310e-54,
            5.0453004471400290e-54,
            3.8628099108590356e-54,
            2.6703466320385890e-54,
            1.6384792977962480e-54,
            8.6857625428972280e-55,
            3.7712641844580550e-55,
            1.2014280428291396e-55,
            2.3703231753300952e-56,
            2.5033357131821380e-57,
        ],
        dtype=float,
    ),
    "nu_mu_bar": np.array(
        [
            1.2288908522909247e-53,
            1.2246875183897874e-53,
            1.2019016165184193e-53,
            1.1386770486361922e-53,
            1.0169553720356749e-53,
            8.3932848205423420e-54,
            6.2908228368291670e-54,
            4.2087622587639140e-54,
            2.4643369411692726e-54,
            1.2272358862992362e-54,
            4.9348505174051480e-55,
            1.4490668011922176e-55,
            2.6714017977489128e-56,
            2.6374938314442993e-57,
        ],
        dtype=float,
    ),
    "nu_e": np.array(
        [
            8.2714757416982840e-54,
            8.2403771338749860e-54,
            8.0765890803536580e-54,
            7.6254199892185330e-54,
            6.7643634573169370e-54,
            5.5319737161540240e-54,
            4.1052964116884845e-54,
            2.7223099673138460e-54,
            1.5844038739959798e-54,
            7.8762062992377100e-55,
            3.1753844906308353e-55,
            9.3660281875747190e-56,
            1.7326448086124229e-56,
            1.7371455206381457e-57,
        ],
        dtype=float,
    ),
    "nu_mu": np.array(
        [
            1.1814197563774684e-53,
            1.1812224385668820e-53,
            1.1734995297165103e-53,
            1.1362332439849290e-53,
            1.0418880986672666e-53,
            8.7638014968501060e-54,
            6.6191676848544680e-54,
            4.4085591531707394e-54,
            2.5366107622192760e-54,
            1.2223406660336456e-54,
            4.6669610850567590e-55,
            1.2939243889986000e-55,
            2.3659614584363346e-56,
            2.4191415535843250e-57,
        ],
        dtype=float,
    ),
    "e_minus": np.array(
        [
            4.2212393772056680e-54,
            4.2211910689841880e-54,
            4.1956168346449824e-54,
            4.0358704401853110e-54,
            3.6365835588820640e-54,
            2.9724292952794440e-54,
            2.1540104673658465e-54,
            1.3628519219255365e-54,
            7.4095265580002700e-55,
            3.3845827578727625e-55,
            1.2469769397452772e-55,
            3.4425354960477750e-56,
            6.4567605248552180e-57,
            7.0869677543587640e-58,
        ],
        dtype=float,
    ),
    "nu_e_bar": np.array(
        [
            4.5635454532812345e-54,
            4.5634743209861420e-54,
            4.5305116077773260e-54,
            4.3334746286958750e-54,
            3.8561410286118525e-54,
            3.0868520928807380e-54,
            2.1805095137167143e-54,
            1.3451138267454154e-54,
            7.1630299044113300e-55,
            3.2315839731157557e-55,
            1.1883896124073048e-55,
            3.3063415431404980e-56,
            6.2907044980648290e-57,
            7.1000480086764130e-58,
        ],
        dtype=float,
    ),
}
GAMMA_SPECTRUM_REFERENCE = np.array(
    [
        4.5694420038138990e-54,
        4.5694420038138990e-54,
        4.5694420038138990e-54,
        4.5621252271843844e-54,
        4.4821155599690133e-54,
        4.2399165061077647e-54,
        3.7796951485144440e-54,
        3.1232773788398700e-54,
        2.3635468758374300e-54,
        1.6101901958862575e-54,
        9.5236228336299250e-55,
        4.4768459558686850e-55,
        1.3865761368557154e-55,
        2.4194450979389490e-56,
    ],
    dtype=float,
)

LEPTON_PHI_TOLERANCE = 5.0e-11
LEPTON_SPECTRUM_TOLERANCE = 5.0e-11


def _relative_error(actual: np.ndarray, reference: np.ndarray) -> np.ndarray:
    denom = np.maximum(np.abs(reference), 1.0e-300)
    return np.abs(actual - reference) / denom


def _secondary_input() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    photon_energy_gev = np.logspace(-5.0, -4.3, 256)
    photon_density_gev = 1.0e-12 * (photon_energy_gev / photon_energy_gev[0]) ** (-1.5)
    proton_energy_gev = np.logspace(4.5, 5.05, 96)
    proton_density_gev = 1.0e-20 * (proton_energy_gev / proton_energy_gev[0]) ** (-2.2)
    return photon_energy_gev, photon_density_gev, proton_energy_gev, proton_density_gev


def evaluate_lepton_reference() -> dict[str, object]:
    photon_energy_gev, photon_density_gev, proton_energy_gev, proton_density_gev = _secondary_input()

    phi_actual: dict[str, np.ndarray] = {}
    spectrum_actual: dict[str, np.ndarray] = {}
    phi_metrics: dict[str, float] = {}
    spectrum_metrics: dict[str, float] = {}

    for channel in CHANNEL_ORDER:
        eta = np.full_like(PHI_X_GRID, CHANNEL_ETA_RATIO[channel] * ETA0)
        phi = kelner_aharonian_2008_secondary_phi(channel, eta, PHI_X_GRID)
        phi_actual[channel] = phi
        phi_metrics[channel] = float(np.max(_relative_error(phi, PHI_REFERENCE[channel])))

        spectrum = kelner_aharonian_2008_secondary_spectrum(
            channel,
            SPECTRUM_ENERGY_GRID,
            proton_energy_gev,
            proton_density_gev,
            photon_energy_gev,
            photon_density_gev,
        )
        spectrum_actual[channel] = spectrum
        spectrum_metrics[channel] = float(np.max(_relative_error(spectrum, SPECTRUM_REFERENCE[channel])))

    gamma_spectrum = kelner_aharonian_2008_gamma_spectrum(
        SPECTRUM_ENERGY_GRID,
        proton_energy_gev,
        proton_density_gev,
        photon_energy_gev,
        photon_density_gev,
    )
    gamma_rel = _relative_error(gamma_spectrum, GAMMA_SPECTRUM_REFERENCE)

    metrics = {
        "phi_max_rel_error": float(max(phi_metrics.values())),
        "phi_mean_rel_error": float(np.mean(list(phi_metrics.values()))),
        "spectrum_max_rel_error": float(max(spectrum_metrics.values())),
        "spectrum_mean_rel_error": float(np.mean(list(spectrum_metrics.values()))),
        "gamma_spectrum_max_rel_error": float(np.max(gamma_rel)),
        "channel_phi_max_rel": phi_metrics,
        "channel_spectrum_max_rel": spectrum_metrics,
    }
    return {
        "reference": {
            "x_grid": PHI_X_GRID,
            "phi": PHI_REFERENCE,
            "spectrum_energy_grid": SPECTRUM_ENERGY_GRID,
            "spectrum": SPECTRUM_REFERENCE,
            "gamma_spectrum": GAMMA_SPECTRUM_REFERENCE,
            "channel_eta_ratio": CHANNEL_ETA_RATIO,
        },
        "actual": {
            "x_grid": PHI_X_GRID,
            "phi": phi_actual,
            "spectrum_energy_grid": SPECTRUM_ENERGY_GRID,
            "spectrum": spectrum_actual,
            "gamma_spectrum": gamma_spectrum,
        },
        "metrics": metrics,
    }


def assert_lepton_reference(payload: dict[str, object]) -> None:
    metrics = payload["metrics"]
    if float(metrics["phi_max_rel_error"]) > LEPTON_PHI_TOLERANCE:
        raise AssertionError(
            f"KA2008 lepton phi regression max_rel_error={metrics['phi_max_rel_error']:.3e} exceeds {LEPTON_PHI_TOLERANCE:.3e}"
        )
    if float(metrics["spectrum_max_rel_error"]) > LEPTON_SPECTRUM_TOLERANCE:
        raise AssertionError(
            f"KA2008 lepton spectrum regression max_rel_error={metrics['spectrum_max_rel_error']:.3e} exceeds {LEPTON_SPECTRUM_TOLERANCE:.3e}"
        )
    if float(metrics["gamma_spectrum_max_rel_error"]) > LEPTON_SPECTRUM_TOLERANCE:
        raise AssertionError(
            f"KA2008 gamma spectrum regression max_rel_error={metrics['gamma_spectrum_max_rel_error']:.3e} exceeds {LEPTON_SPECTRUM_TOLERANCE:.3e}"
        )

    for channel in CHANNEL_ORDER:
        phi = payload["actual"]["phi"][channel]
        spectrum = payload["actual"]["spectrum"][channel]
        if not np.all(np.isfinite(phi)) or not np.all(phi >= 0.0):
            raise AssertionError(f"KA2008 channel {channel} returned invalid phi values.")
        if not np.all(np.isfinite(spectrum)) or not np.all(spectrum >= 0.0):
            raise AssertionError(f"KA2008 channel {channel} returned invalid spectrum values.")


def main() -> None:
    payload = evaluate_lepton_reference()
    metrics = payload["metrics"]
    print(f"lepton_phi_max_rel_error={metrics['phi_max_rel_error']:.3e}")
    print(f"lepton_spectrum_max_rel_error={metrics['spectrum_max_rel_error']:.3e}")
    print(f"gamma_spectrum_max_rel_error={metrics['gamma_spectrum_max_rel_error']:.3e}")
    assert_lepton_reference(payload)


if __name__ == "__main__":
    main()
