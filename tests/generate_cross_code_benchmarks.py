"""Generate formal source-data tables for cross-code benchmarks.

This script records the fixed comparison environment and writes measured
forward-shock synchrotron fluxes where a local adapter is already defined.
"""

from __future__ import annotations

import csv
import importlib.metadata as meta
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import h5py
import numpy as np


root = Path(__file__).resolve().parents[1]
if str(root) not in sys.path:
    sys.path.insert(0, str(root))

from asgard_core import Model, Observer, ObserverGrid, UniformMedium, gaussian_jet
from tests.public_api_builders import (
    hadronic,
    numerics,
    observer_grid,
    radiation,
    reverse_shock,
    solver_options,
    top_hat_model,
)

data = root / "paper" / "source_data"
data.mkdir(parents=True, exist_ok=True)


def version(name: str) -> str:
    return meta.version(name)


def git_head() -> str:
    return subprocess.check_output(
        ["git", "rev-parse", "--short", "HEAD"],
        cwd=root,
        text=True,
        encoding="utf-8",
    ).strip()


times = np.logspace(0.0, 8.0, 17, dtype=float)
freqs = np.array(
    [1.0e8, 1.0e9, 1.0e10, 1.0e12, 1.0e14, 5.0e14, 1.0e17, 1.0e18, 2.4e22, 2.4e23, 2.4e26],
    dtype=float,
)
rs_times = np.logspace(0.0, 7.0, 15, dtype=float)
rs_freqs = np.array([1.0e8, 1.0e9, 1.0e12, 5.0e14, 1.0e18, 2.4e23, 2.4e26], dtype=float)
geom_times = np.logspace(1.0, 8.0, 15, dtype=float)
geom_freqs = np.array([1.0e8, 1.0e9, 5.0e14, 1.0e18, 2.4e23], dtype=float)
mjy_cgs = 1.0e-26
c_cgs = 2.99792458e10
m_e_g = 9.1093837015e-28
m_p_g = 1.67262192369e-24
mec2_erg = m_e_g * c_cgs**2
mpc2_erg = m_p_g * c_cgs**2


def run_asgard() -> tuple[np.ndarray, float]:
    model = top_hat_model(
        numerics=numerics(
            num_radius=192,
            num_observer_time=48,
            num_electron_gamma=128,
            num_photon_frequency=96,
            num_threads=1,
        ),
        observer_grid=observer_grid(time_min_s=float(times[0]), time_max_s=float(times[-1])),
    )
    start = time.perf_counter()
    flux = np.asarray(model.flux_density_grid(times, freqs).total, dtype=float)
    return flux, time.perf_counter() - start


def run_afterglowpy() -> tuple[np.ndarray, float]:
    import afterglowpy as grb

    params = dict(
        jetType=grb.jet.TopHat,
        specType=grb.jet.SimpleSpec,
        thetaObs=0.0,
        E0=1.0e52,
        thetaCore=0.1,
        thetaWing=0.1,
        b=0.0,
        L0=0.0,
        q=0.0,
        ts=0.0,
        n0=1.0,
        p=2.3,
        epsilon_e=0.1,
        epsilon_B=1.0e-3,
        xi_N=0.1,
        d_L=1.0e26,
        z=0.1,
        g0=300.0,
        spread=False,
    )
    pair_t = np.repeat(times, freqs.size)
    pair_nu = np.tile(freqs, times.size)
    start = time.perf_counter()
    flux = np.asarray(grb.fluxDensity(pair_t, pair_nu, **params), dtype=float)
    return flux.reshape(times.size, freqs.size).T * mjy_cgs, time.perf_counter() - start


def run_jetsimpy() -> tuple[np.ndarray, float]:
    import jetsimpy

    mpc_cm = 3.0856775814913673e24
    params = dict(
        Eiso=1.0e52,
        lf=300.0,
        theta_c=0.1,
        n0=1.0,
        k=0.0,
        A=0.0,
        eps_e=0.1,
        eps_b=1.0e-3,
        p=2.3,
        theta_v=0.0,
        d=1.0e26 / mpc_cm,
        z=0.1,
        b=0.0,
    )
    start = time.perf_counter()
    cols = [
        np.asarray(
            jetsimpy.FluxDensity_tophat(
                times,
                float(nu),
                params,
                tmin=1.0,
                tmax=1.0e9,
                spread=False,
                cal_level=2,
                rtol=5.0e-4,
                max_iter=220,
            ),
            dtype=float,
        )
        for nu in freqs
    ]
    return np.vstack(cols) * mjy_cgs, time.perf_counter() - start


def run_vegas() -> tuple[np.ndarray, float]:
    from VegasAfterglow import ISM, Model, Observer, Radiation, TophatJet

    medium = ISM(n_ism=1.0)
    jet = TophatJet(theta_c=0.1, E_iso=1.0e52, Gamma0=300.0, duration=10.0)
    observer = Observer(lumi_dist=1.0e26, z=0.1, theta_obs=0.0)
    rad = Radiation(eps_e=0.1, eps_B=1.0e-3, p=2.3, xi_e=0.1, ssc=False, kn=True)
    model = Model(
        jet=jet,
        medium=medium,
        observer=observer,
        fwd_rad=rad,
        rvs_rad=None,
        resolutions=(0.035, 0.070, 24),
    )
    start = time.perf_counter()
    flux = np.asarray(model.flux_density_grid(times, freqs).total, dtype=float)
    return flux, time.perf_counter() - start


def pyblast_path() -> Path:
    src = root / "_external" / "PyBlastAfterglowMag"
    link = Path("/tmp/asgard_pyblast_repo")
    if not src.exists():
        raise FileNotFoundError(f"PyBlastAfterglowMag source tree is missing: {src}")
    if link.is_symlink():
        link.unlink()
    elif link.exists():
        raise FileExistsError(f"PyBlastAfterglow run path exists and is not a symlink: {link}")
    os.symlink(src, link, target_is_directory=True)
    exe = link / "src" / "pba.out"
    if not exe.exists():
        raise FileNotFoundError(f"PyBlastAfterglow executable is missing: {exe}")
    return link


def pyblast_lc(work: Path, nu_grid: np.ndarray, obs_grid: np.ndarray) -> np.ndarray:
    with h5py.File(work / "lc.h5", "r") as handle:
        h5_nu = np.asarray(handle["freqs"], dtype=float)
        h5_obs = np.asarray(handle["times"], dtype=float)
        h5_flux = np.zeros_like(h5_nu, dtype=float)
        for name in handle:
            if name.startswith("shell="):
                h5_flux += np.asarray(handle[name]["fluxdens"], dtype=float)
    flux = np.empty((nu_grid.size, obs_grid.size), dtype=float)
    for i, nu in enumerate(nu_grid):
        for j, obs in enumerate(obs_grid):
            mask = np.isclose(h5_nu, nu, rtol=1.0e-12, atol=0.0) & np.isclose(h5_obs, obs, rtol=1.0e-12, atol=0.0)
            if np.count_nonzero(mask) != 1:
                raise RuntimeError(f"PyBlastAfterglow lc.h5 has no unique cell for nu={nu} t={obs}")
            flux[i, j] = float(h5_flux[mask][0]) * mjy_cgs
    return flux


def run_pyblast(kind: str, nu_grid: np.ndarray, obs_grid: np.ndarray) -> tuple[np.ndarray, float]:
    from PyBlastAfterglowMag.id_analytic import JetStruct
    from PyBlastAfterglowMag.parfile_tools import create_parfile

    pbaroot = pyblast_path()
    work = Path("/tmp") / f"asgard_pyblast_{kind}"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    if kind == "tophat":
        struct = {
            "struct": "tophat",
            "Eiso_c": 1.0e52,
            "Gamma0c": 300.0,
            "M0c": -1.0,
            "theta_c": 0.1,
            "theta_w": 0.1,
        }
        layers = 1
        theta_obs = 0.0
    elif kind == "gaussian":
        struct = {
            "struct": "gaussian",
            "Eiso_c": 1.0e52,
            "Gamma0c": 120.0,
            "M0c": -1.0,
            "theta_c": 0.08,
            "theta_w": 0.24,
        }
        layers = 10
        theta_obs = 0.12
    else:
        raise ValueError(f"unknown PyBlastAfterglow benchmark kind: {kind}")
    ids = JetStruct(n_layers_pw=layers, n_layers_a=layers)
    for mode, name in [("piece-wise", "id_pw.h5"), ("adaptive", "id_a.h5")]:
        iddict, idpars = ids.get_1D_id(pars=struct, type=mode)
        ids.save_1d_id(id_dict=iddict, id_pars=idpars, outfpath=str(work / name))
    main = {
        "d_l": 1.0e26,
        "z": 0.1,
        "n_ism": 1.0,
        "theta_obs": theta_obs,
        "rtol": 1.0e-3,
        "lc_freqs": "array " + " ".join(f"{float(v):.16e}" for v in nu_grid),
        "lc_times": "array " + " ".join(f"{float(v):.16e}" for v in obs_grid),
        "tb0": 1.0,
        "tb1": 1.0e9,
        "ntb": 256,
    }
    grb = {
        "eps_e_fs": 0.1,
        "eps_b_fs": 1.0e-3,
        "p_fs": 2.3,
        "do_lc": "yes",
        "save_spec": "no",
        "do_skymap": "no",
        "do_rs": "no",
        "do_rs_radiation": "no",
        "method_synchrotron_fs": "Joh06",
        "method_ele_fs": "analytic",
        "method_gamma_min_fs": "useU_e",
        "method_spread": "None",
        "method_comp_mode": "observFlux",
        "fname_ejecta_id": "id_a.h5",
        "method_eats": "adaptive",
        "nsublayers": 35,
        "ebl_tbl_fpath": str(pbaroot / "data" / "EBL" / "Franceschini18" / "table.h5"),
    }
    create_parfile(str(work) + "/", {"main": main, "grb": grb})
    with (work / "parfile.par").open("a", encoding="utf-8") as handle:
        handle.write(
            "\n# ------------------------ Magnetar -------------------------\n\n"
            "* Parameters\n\n\n* Settings\n\n"
            "run_magnetar = no\nload_magnetar = no\nsave_magnetar = no\n\n"
            "# --------------------------- END ---------------------------\n\n"
        )
    start = time.perf_counter()
    subprocess.check_call([str(pbaroot / "src" / "pba.out"), str(work) + "/", "parfile.par", "err"])
    wall = time.perf_counter() - start
    return pyblast_lc(work, nu_grid, obs_grid), wall


def run_asgard_rs_ssc() -> tuple[dict[str, np.ndarray], float]:
    model = top_hat_model(
        fwd_rad=radiation(include_ssc=True, include_kn_correction=True),
        rvs_rad=radiation(include_ssc=True, include_kn_correction=True),
        reverse_shock=reverse_shock(enabled=True, include_ssc=True),
        numerics=numerics(
            num_radius=192,
            num_observer_time=48,
            num_electron_gamma=128,
            num_photon_frequency=96,
            num_threads=1,
        ),
        observer_grid=observer_grid(time_min_s=float(rs_times[0]), time_max_s=float(rs_times[-1])),
    )
    start = time.perf_counter()
    res = model.flux_density_grid(rs_times, rs_freqs)
    wall = time.perf_counter() - start
    return {
        "fs_sync": np.asarray(res.fwd.sync, dtype=float),
        "fs_ssc": np.asarray(res.fwd.ssc, dtype=float),
        "rs_sync": np.asarray(res.rev.sync, dtype=float),
        "rs_ssc": np.asarray(res.rev.ssc, dtype=float),
    }, wall


def run_vegas_rs_ssc() -> tuple[dict[str, np.ndarray], float]:
    from VegasAfterglow import ISM, Model, Observer, Radiation, TophatJet

    medium = ISM(n_ism=1.0)
    jet = TophatJet(theta_c=0.1, E_iso=1.0e52, Gamma0=300.0, duration=10.0)
    observer = Observer(lumi_dist=1.0e26, z=0.1, theta_obs=0.0)
    rad = Radiation(eps_e=0.1, eps_B=1.0e-3, p=2.3, xi_e=0.1, ssc=True, kn=True)
    model = Model(
        jet=jet,
        medium=medium,
        observer=observer,
        fwd_rad=rad,
        rvs_rad=rad,
        resolutions=(0.035, 0.070, 24),
    )
    start = time.perf_counter()
    res = model.flux_density_grid(rs_times, rs_freqs)
    wall = time.perf_counter() - start
    return {
        "fs_sync": np.asarray(res.fwd.sync, dtype=float),
        "fs_ssc": np.asarray(res.fwd.ssc, dtype=float),
        "rs_sync": np.asarray(res.rvs.sync, dtype=float),
        "rs_ssc": np.asarray(res.rvs.ssc, dtype=float),
    }, wall


def run_asgard_gaussian() -> tuple[np.ndarray, float]:
    model = Model(
        jet=gaussian_jet(1.0e52, 120.0, 0.08, 0.24, None, None, False),
        medium=UniformMedium(1.0),
        observer=Observer(
            z=0.1,
            viewing_angle_rad=0.12,
            viewing_azimuth_rad=0.0,
            luminosity_distance_cm=1.0e26,
        ),
        fwd_rad=radiation(include_ssc=False),
        rvs_rad=None,
        numerics=numerics(
            structured_num_theta=32,
            structured_num_phi=64,
            num_radius=96,
            eats_num_theta=24,
            eats_num_phi=24,
            num_observer_time=40,
            num_electron_gamma=96,
            num_photon_frequency=64,
            num_threads=1,
        ),
        observer_grid=ObserverGrid(float(geom_times[0]), float(geom_times[-1])),
        solver_options=solver_options(),
        reverse_shock=reverse_shock(),
        hadronic=hadronic(),
    )
    start = time.perf_counter()
    flux = np.asarray(model.flux_density_grid(geom_times, geom_freqs).total, dtype=float)
    return flux, time.perf_counter() - start


def run_afterglowpy_gaussian() -> tuple[np.ndarray, float]:
    import afterglowpy as grb

    params = dict(
        jetType=grb.jet.Gaussian,
        specType=grb.jet.SimpleSpec,
        thetaObs=0.12,
        E0=1.0e52,
        thetaCore=0.08,
        thetaWing=0.24,
        b=0.0,
        L0=0.0,
        q=0.0,
        ts=0.0,
        n0=1.0,
        p=2.3,
        epsilon_e=0.1,
        epsilon_B=1.0e-3,
        xi_N=0.1,
        d_L=1.0e26,
        z=0.1,
        g0=120.0,
        spread=False,
    )
    pair_t = np.repeat(geom_times, geom_freqs.size)
    pair_nu = np.tile(geom_freqs, geom_times.size)
    start = time.perf_counter()
    flux = np.asarray(grb.fluxDensity(pair_t, pair_nu, **params), dtype=float)
    return flux.reshape(geom_times.size, geom_freqs.size).T * mjy_cgs, time.perf_counter() - start


def run_jetsimpy_gaussian() -> tuple[np.ndarray, float]:
    import jetsimpy

    mpc_cm = 3.0856775814913673e24
    params = dict(
        Eiso=1.0e52,
        lf=120.0,
        theta_c=0.08,
        theta_v=0.12,
        n0=1.0,
        k=0.0,
        A=0.0,
        eps_e=0.1,
        eps_b=1.0e-3,
        p=2.3,
        d=1.0e26 / mpc_cm,
        z=0.1,
    )
    start = time.perf_counter()
    cols = [
        np.asarray(
            jetsimpy.FluxDensity_gaussian(
                geom_times,
                float(nu),
                params,
                tmin=1.0,
                tmax=1.0e9,
                spread=False,
                cal_level=2,
                rtol=5.0e-4,
                max_iter=220,
            ),
            dtype=float,
        )
        for nu in geom_freqs
    ]
    return np.vstack(cols) * mjy_cgs, time.perf_counter() - start


def run_vegas_gaussian() -> tuple[np.ndarray, float]:
    from VegasAfterglow import GaussianJet, ISM, Model as VegasModel, Observer as VegasObserver, Radiation as VegasRadiation

    medium = ISM(n_ism=1.0)
    jet = GaussianJet(theta_c=0.08, E_iso=1.0e52, Gamma0=120.0)
    observer = VegasObserver(lumi_dist=1.0e26, z=0.1, theta_obs=0.12)
    rad = VegasRadiation(eps_e=0.1, eps_B=1.0e-3, p=2.3, xi_e=0.1, ssc=False, kn=True)
    model = VegasModel(
        jet=jet,
        medium=medium,
        observer=observer,
        fwd_rad=rad,
        rvs_rad=None,
        resolutions=(0.035, 0.070, 24),
    )
    start = time.perf_counter()
    flux = np.asarray(model.flux_density_grid(geom_times, geom_freqs).total, dtype=float)
    return flux, time.perf_counter() - start


def electron_budget(gamma: np.ndarray, dnde: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    dgamma = np.gradient(gamma)
    number = np.sum(dnde * dgamma[:, None], axis=0)
    energy = np.sum((gamma[:, None] - 1.0) * dnde * dgamma[:, None], axis=0) * mec2_erg
    return number, energy


def proton_budget(gamma: np.ndarray, dndp: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    dgamma = np.gradient(gamma)
    number = np.sum(dndp * dgamma[:, None], axis=0)
    energy = np.sum((gamma[:, None] - 1.0) * dndp * dgamma[:, None], axis=0) * mpc2_erg
    return number, energy


def run_asgard_complex_state() -> tuple[list[dict[str, str]], float]:
    start = time.perf_counter()
    base_model = top_hat_model(
        fwd_rad=radiation(include_ssc=False),
        numerics=numerics(
            num_radius=96,
            eats_num_theta=8,
            num_observer_time=40,
            num_electron_gamma=96,
            num_photon_frequency=64,
            num_threads=1,
        ),
    )
    base = base_model.details(1.0e4, 2.0e4).fwd
    chi_model = top_hat_model(
        fwd_rad=radiation(include_ssc=False),
        numerics=numerics(
            num_radius=64,
            eats_num_theta=12,
            num_observer_time=32,
            num_electron_gamma=64,
            num_photon_frequency=64,
            downstream_num_chi=24,
            num_threads=1,
        ),
        solver_options=solver_options(
            electron_solver="fullhide_2d",
            geometry_projection="chi_eats_2d",
        ),
    )
    chi = chi_model.details(1.0e4, 2.0e4).fwd
    hadronic_model = top_hat_model(
        fwd_rad=radiation(
            proton_energy_fraction=0.01,
            include_pgamma=False,
            bethe_heitler=False,
            pair_production=False,
        ),
        hadronic=hadronic(
            enabled=True,
            num_proton_gamma=128,
            num_neutrino_frequency=64,
            pgamma_scheme="disabled",
            pair_cascade_iterations=1,
        ),
        numerics=numerics(
            num_radius=96,
            eats_num_theta=12,
            num_observer_time=32,
            num_electron_gamma=64,
            num_photon_frequency=64,
            num_threads=1,
        ),
    )
    had = hadronic_model.details(1.0e4, 2.0e4).fwd
    pair_model = top_hat_model(
        fwd_rad=radiation(
            proton_energy_fraction=0.01,
            include_pgamma=True,
            bethe_heitler=True,
            pair_production=True,
            pgamma_scheme="hummer_2010_response",
        ),
        hadronic=hadronic(
            enabled=True,
            solver="am3_1d",
            num_proton_gamma=96,
            num_neutrino_frequency=64,
            pgamma_scheme="hummer_2010_response",
            pair_cascade_iterations=2,
        ),
        solver_options=solver_options(
            electron_photon_coupling="joint",
            ssc_cooling_mode="numeric_ic_kn",
        ),
        numerics=numerics(
            num_radius=64,
            eats_num_theta=12,
            num_observer_time=32,
            num_electron_gamma=64,
            num_photon_frequency=64,
            num_threads=1,
            electron_adaptive_substeps=False,
        ),
    )
    pair = pair_model.details(1.0e4, 2.0e4).fwd
    wall = time.perf_counter() - start

    base_radius = np.asarray(base.radius, dtype=float)
    base_gamma = np.asarray(base.Gamma, dtype=float)
    base_b = np.asarray(base.B_comv, dtype=float)
    gamma_e = np.asarray(base.gamma_e, dtype=float)
    dnde = np.asarray(base.dN_dgamma_e, dtype=float)
    electron_number, electron_energy = electron_budget(gamma_e, dnde)

    had_radius = np.asarray(had.radius, dtype=float)
    gamma_p = np.asarray(had.gamma_p, dtype=float)
    dndp = np.asarray(had.dN_dgamma_p, dtype=float)
    proton_number, proton_energy = proton_budget(gamma_p, dndp)
    hadronic_peak = np.nanmax(np.asarray(had.l_had_syn_spec, dtype=float), axis=0)
    pair_radius = np.asarray(pair.radius, dtype=float)
    pair_number, pair_energy = electron_budget(
        np.asarray(pair.gamma_e, dtype=float),
        np.asarray(pair.dN_dgamma_e_bh, dtype=float),
    )
    pg_peak = np.nanmax(np.asarray(pair.l_had_pg_gamma_spec, dtype=float), axis=0)
    bh_peak = np.nanmax(np.asarray(pair.l_had_bethe_heitler_spec, dtype=float), axis=0)
    tau_pg_peak = np.nanmax(np.asarray(pair.tau_pg, dtype=float), axis=0)
    tau_bh_peak = np.nanmax(np.asarray(pair.tau_bh, dtype=float), axis=0)
    pg_survival_min = np.nanmin(np.asarray(pair.pg_photon_survival, dtype=float), axis=0)

    chi_radius = np.asarray(chi.chi_radius_cm, dtype=float)
    chi_grid = np.asarray(chi.chi_grid, dtype=float)
    chi_weight = np.asarray(chi.chi_dvolume_weight, dtype=float)
    tau_chi = np.nanmax(np.asarray(chi.tau_syn_chi, dtype=float), axis=0)
    seed_chi = np.nanmax(np.asarray(chi.seed_syn_chi, dtype=float), axis=0)

    item = codes["ASGARD"]
    rows: list[dict[str, str]] = []

    def add_row(
        metric_group: str,
        value_name: str,
        value: float | str,
        radius: float | str,
        chi_value: float | str = "",
        boundary: str = "ASGARD-only radius-ordered diagnostic; no cross-code ratio",
    ) -> None:
        rows.append(
            {
                "scenario": "complex_state_dashboard",
                "code": "ASGARD",
                "code_version": item["version"],
                "source": item["source"],
                "physical_assumptions": (
                    "ASGARD-only radius-ordered state dashboard; uniform medium; "
                    "proton-synch hadronic rows use legacy_1d; activated p-gamma, "
                    "Bethe-Heitler, and pair-feedback rows use am3_1d with joint coupling; "
                    "chi projection uses the fullhide_2d electron branch"
                ),
                "coordinate": "radius and local state coordinate",
                "frequency_hz": "",
                "observer_time": "",
                "value_name": value_name,
                "value": f"{value:.8e}" if isinstance(value, float) else value,
                "asgard_ratio": "",
                "wall_time": f"{wall:.8e}",
                "boundary": boundary,
                "radius_cm": f"{radius:.8e}" if isinstance(radius, float) else radius,
                "chi": f"{chi_value:.8e}" if isinstance(chi_value, float) else chi_value,
                "metric_group": metric_group,
            }
        )

    for i, radius in enumerate(base_radius):
        add_row("full_state_chain", "Gamma", float(base_gamma[i]), float(radius))
        add_row("full_state_chain", "B_comoving_G", float(base_b[i]), float(radius))
        add_row("full_state_chain", "electron_number", float(electron_number[i]), float(radius))
        add_row("full_state_chain", "electron_energy_erg", float(electron_energy[i]), float(radius))

    for i, radius in enumerate(had_radius):
        add_row("hadronic_state", "proton_number", float(proton_number[i]), float(radius))
        add_row("hadronic_state", "proton_energy_erg", float(proton_energy[i]), float(radius))
        add_row("hadronic_state", "hadronic_synch_peak_cgs", float(hadronic_peak[i]), float(radius))

    for i, radius in enumerate(pair_radius):
        add_row("pair_feedback_state", "bh_pair_number", float(pair_number[i]), float(radius))
        add_row("pair_feedback_state", "bh_pair_energy_erg", float(pair_energy[i]), float(radius))
        add_row("pair_feedback_state", "pgamma_peak_cgs", float(pg_peak[i]), float(radius))
        add_row("pair_feedback_state", "bethe_heitler_peak_cgs", float(bh_peak[i]), float(radius))
        add_row("pair_feedback_state", "max_tau_pg", float(tau_pg_peak[i]), float(radius))
        add_row("pair_feedback_state", "max_tau_bh", float(tau_bh_peak[i]), float(radius))
        add_row("pair_feedback_state", "min_pg_photon_survival", float(pg_survival_min[i]), float(radius))

    for j, radius_col in enumerate(chi_radius.T):
        for a, chi_value in enumerate(chi_grid):
            radius = float(radius_col[a])
            add_row("chi_projection", "chi_dvolume_weight", float(chi_weight[a, j]), radius, float(chi_value))
            add_row("chi_projection", "max_tau_syn_chi", float(tau_chi[a, j]), radius, float(chi_value))
            add_row("chi_projection", "max_seed_syn_chi", float(seed_chi[a, j]), radius, float(chi_value))

    add_row(
        "pair_feedback_boundary",
        "chi_resolved_hadronic_feedback",
        "",
        "",
        "",
        "activated pair-feedback rows are one-dimensional and shell-averaged; chi-resolved hadronic feedback is unsupported",
    )
    return rows, wall


def finite_ratio(value: float, ref: float) -> str:
    if np.isfinite(value) and np.isfinite(ref) and ref > 0.0:
        return f"{value / ref:.8e}"
    return ""


codes = {
    "ASGARD": {
        "version": version("asgard-grb-afterglow"),
        "source": f"local git {git_head()}",
    },
    "afterglowpy": {
        "version": version("afterglowpy"),
        "source": "PyPI afterglowpy==0.8.1",
    },
    "jetsimpy": {
        "version": version("jetsimpy"),
        "source": "git 53d07610830f2247d8c41364a2b469491fa22eb2",
    },
    "VegasAfterglow": {
        "version": version("VegasAfterglow"),
        "source": "PyPI VegasAfterglow==2.0.5",
    },
    "PyBlastAfterglowMag": {
        "version": version("PyBlastAfterglowMag"),
        "source": "git b0a39d170930908a398de56089fde8cfab16883d",
    },
}


base = [
    "scenario",
    "code",
    "code_version",
    "source",
    "physical_assumptions",
    "coordinate",
    "frequency_hz",
    "observer_time",
    "value_name",
    "value",
    "asgard_ratio",
    "wall_time",
    "boundary",
]


def write_csv(name: str, rows: list[dict[str, str]], fields: list[str] = base) -> None:
    with (data / name).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


asgard_flux, asgard_wall = run_asgard()
measured = {
    "ASGARD": (
        asgard_flux,
        asgard_wall,
        "broad-band broad-time FS synchrotron grid; ratios are assumption-audit diagnostics",
    ),
    "afterglowpy": (
        *run_afterglowpy(),
        "afterglowpy returns mJy; values converted to cgs; ratios are assumption-audit diagnostics",
    ),
    "jetsimpy": (
        *run_jetsimpy(),
        "jetsimpy returns mJy; values converted to cgs; ratios are assumption-audit diagnostics",
    ),
    "VegasAfterglow": (
        *run_vegas(),
        "broad-band broad-time FS synchrotron grid; ratios are assumption-audit diagnostics",
    ),
    "PyBlastAfterglowMag": (
        *run_pyblast("tophat", freqs, times),
        "PyBlastAfterglowMag parfile plus repo-local pba.out; mJy values converted to cgs; ratios are assumption-audit diagnostics",
    ),
}

fs_rows = []
for code, item in codes.items():
    if code in measured:
        flux, wall, note = measured[code]
        for i, nu in enumerate(freqs):
            for j, obs in enumerate(times):
                fs_rows.append(
                    {
                        "scenario": "fs_synch_common_minimum",
                        "code": code,
                        "code_version": item["version"],
                        "source": item["source"],
                        "physical_assumptions": (
                            "broad-band top-hat on-axis forward-shock synchrotron; uniform medium; "
                            "E_iso=1e52 erg; Gamma0=300; theta_j=0.1 rad; n=1 cm^-3; "
                            "epsilon_e=0.1; epsilon_B=1e-3; p=2.3; xi_e=0.1; z=0.1; "
                            "d_L=1e26 cm; no spreading; no SSC; no reverse shock"
                        ),
                        "coordinate": "observer time and frequency",
                        "frequency_hz": f"{nu:.8e}",
                        "observer_time": f"{obs:.8e}",
                        "value_name": "flux_density_cgs",
                        "value": f"{float(flux[i, j]):.8e}",
                        "asgard_ratio": "1.00000000e+00"
                        if code == "ASGARD"
                        else finite_ratio(float(flux[i, j]), float(asgard_flux[i, j])),
                        "wall_time": f"{wall:.8e}",
                        "boundary": note,
                    }
                )

asgard_rs, asgard_rs_wall = run_asgard_rs_ssc()
vegas_rs, vegas_rs_wall = run_vegas_rs_ssc()

fig2_rows = []
for code, item, payload, wall in [
    ("ASGARD", codes["ASGARD"], asgard_rs, asgard_rs_wall),
    ("VegasAfterglow", codes["VegasAfterglow"], vegas_rs, vegas_rs_wall),
]:
    for component, arr in payload.items():
        ref = asgard_rs[component]
        for i, nu in enumerate(rs_freqs):
            for j, obs in enumerate(rs_times):
                value = float(arr[i, j])
                fig2_rows.append(
                    {
                        "scenario": "rs_ssc_matched_tophat",
                        "code": code,
                        "code_version": item["version"],
                        "source": item["source"],
                        "physical_assumptions": (
                            "top-hat on-axis reverse-shock and SSC benchmark; uniform medium; "
                            "E_iso=1e52 erg; Gamma0=300; theta_j=0.1 rad; shell_duration=10 s; "
                            "n=1 cm^-3; epsilon_e=0.1; epsilon_B=1e-3; p=2.3; xi_e=0.1; "
                            "z=0.1; d_L=1e26 cm"
                        ),
                        "coordinate": "observer time and frequency",
                        "frequency_hz": f"{nu:.8e}",
                        "observer_time": f"{obs:.8e}",
                        "value_name": f"{component}_flux_density_cgs",
                        "value": f"{value:.8e}",
                        "asgard_ratio": "1.00000000e+00"
                        if code == "ASGARD"
                        else finite_ratio(value, float(ref[i, j])),
                        "wall_time": f"{wall:.8e}",
                        "boundary": "ASGARD and VegasAfterglow matched RS/SSC subset only",
                    }
                )

for code, item, boundary in [
    (
        "afterglowpy",
        codes["afterglowpy"],
        "reverse shock and SSC are unsupported in the matched benchmark",
    ),
    (
        "jetsimpy",
        codes["jetsimpy"],
        "reverse shock and SSC are unsupported in the public branch",
    ),
    (
        "PyBlastAfterglowMag",
        codes["PyBlastAfterglowMag"],
        "RS/SSC subset is not part of the current PyBlastAfterglowMag matched benchmark adapter",
    ),
]:
    fig2_rows.append(
        {
            "scenario": "rs_ssc_matched_tophat",
            "code": code,
            "code_version": item["version"],
            "source": item["source"],
            "physical_assumptions": "top-hat on-axis reverse-shock and SSC benchmark",
            "coordinate": "observer time and frequency",
            "frequency_hz": "",
            "observer_time": "",
            "value_name": "component_flux_density_cgs",
            "value": "",
            "asgard_ratio": "",
            "wall_time": "",
            "boundary": boundary,
        }
    )


asgard_geom, asgard_geom_wall = run_asgard_gaussian()
geom_measured = {
    "ASGARD": (
        asgard_geom,
        asgard_geom_wall,
        "structured/off-axis FS synchrotron subset; Gaussian wing truncated at theta_max=0.24 rad",
    ),
    "afterglowpy": (
        *run_afterglowpy_gaussian(),
        "afterglowpy returns mJy; values converted to cgs; Gaussian wing set by thetaWing=0.24 rad",
    ),
    "jetsimpy": (
        *run_jetsimpy_gaussian(),
        "jetsimpy returns mJy; values converted to cgs; public Gaussian adapter has no explicit theta_max field",
    ),
    "VegasAfterglow": (
        *run_vegas_gaussian(),
        "VegasAfterglow public Gaussian adapter has no explicit theta_max field",
    ),
    "PyBlastAfterglowMag": (
        *run_pyblast("gaussian", geom_freqs, geom_times),
        "PyBlastAfterglowMag parfile plus repo-local pba.out; mJy values converted to cgs; Gaussian wing set by theta_w=0.24 rad",
    ),
}

for code, item in codes.items():
    if code in geom_measured:
        flux, wall, note = geom_measured[code]
        for i, nu in enumerate(geom_freqs):
            for j, obs in enumerate(geom_times):
                value = float(flux[i, j])
                fig2_rows.append(
                    {
                        "scenario": "gaussian_off_axis_fs_synch",
                        "code": code,
                        "code_version": item["version"],
                        "source": item["source"],
                        "physical_assumptions": (
                            "Gaussian structured jet; off-axis observer; FS synchrotron only; "
                            "uniform medium; E_iso=1e52 erg; Gamma0=120; theta_c=0.08 rad; "
                            "theta_obs=0.12 rad; n=1 cm^-3; epsilon_e=0.1; "
                            "epsilon_B=1e-3; p=2.3; xi_e=0.1; z=0.1; d_L=1e26 cm; no spreading"
                        ),
                        "coordinate": "observer time and frequency",
                        "frequency_hz": f"{nu:.8e}",
                        "observer_time": f"{obs:.8e}",
                        "value_name": "fs_synch_flux_density_cgs",
                        "value": f"{value:.8e}",
                        "asgard_ratio": "1.00000000e+00"
                        if code == "ASGARD"
                        else finite_ratio(value, float(asgard_geom[i, j])),
                        "wall_time": f"{wall:.8e}",
                        "boundary": note,
                    }
                )

fig3_rows, _complex_wall = run_asgard_complex_state()
for code, item, boundary in [
    (
        "afterglowpy",
        codes["afterglowpy"],
        "unsupported for ASGARD hadronic, pair, chi-resolved, and radius-ordered state dashboard",
    ),
    (
        "jetsimpy",
        codes["jetsimpy"],
        "unsupported for ASGARD hadronic, pair, chi-resolved, and radius-ordered state dashboard",
    ),
    (
        "VegasAfterglow",
        codes["VegasAfterglow"],
        "unsupported for ASGARD hadronic, pair, chi-resolved, and radius-ordered state dashboard",
    ),
    (
        "PyBlastAfterglowMag",
        codes["PyBlastAfterglowMag"],
        "unsupported for ASGARD hadronic, pair, chi-resolved, and radius-ordered state dashboard",
    ),
]:
    fig3_rows.append(
        {
            "scenario": "complex_state_dashboard",
            "code": code,
            "code_version": item["version"],
            "source": item["source"],
            "physical_assumptions": "ASGARD-only local-state diagnostic",
            "coordinate": "radius and local state coordinate",
            "frequency_hz": "",
            "observer_time": "",
            "value_name": "state diagnostic",
            "value": "",
            "asgard_ratio": "",
            "wall_time": "",
            "boundary": boundary,
            "radius_cm": "",
            "chi": "",
            "metric_group": "unsupported_boundary",
        }
    )

fig3_fields = base + ["radius_cm", "chi", "metric_group"]


table_fields = base + [
    "physical_channels",
    "geometry_media",
    "exported_states",
    "benchmark_evidence",
    "explicit_boundaries",
]

cap_rows = [
    {
        **{
            "scenario": "cross_code_capability",
            "code": "ASGARD",
            "code_version": codes["ASGARD"]["version"],
            "source": codes["ASGARD"]["source"],
            "physical_assumptions": "radius-ordered local state before observer projection",
            "coordinate": "radius, frequency, and observer time",
            "frequency_hz": "",
            "observer_time": "",
            "value_name": "capability entry",
            "value": "",
            "asgard_ratio": "",
            "wall_time": "",
            "boundary": "formal ratios available only for explicitly matched subsets",
        },
        "physical_channels": (
            "FS synchrotron; FS SSC or IC; optional RS; optional shell-averaged "
            "hadronic and pair feedback; absorption; observer assembly"
        ),
        "geometry_media": "uniform medium; k=2 stellar wind; tabulated medium; structured components",
        "exported_states": "dynamics, particles, photons, components, and diagnostics",
        "benchmark_evidence": (
            "paper figures 1--5 and figA1; FigB1 broad-band FS rows; "
            "FigB2 ASGARD/Vegas RS+SSC subset and Gaussian off-axis FS subset; "
            "FigB3 ASGARD-only state dashboard"
        ),
        "explicit_boundaries": (
            "no chi-resolved hadronic transport; no IC-mediated electromagnetic cascade; "
            "no jet-spreading backend; no custom Medium dispatch; no wind k != 2; "
            "no active sky-image or polarization backend"
        ),
    },
    {
        **{
            "scenario": "cross_code_capability",
            "code": "afterglowpy",
            "code_version": codes["afterglowpy"]["version"],
            "source": codes["afterglowpy"]["source"],
            "physical_assumptions": "public forward-shock synchrotron model",
            "coordinate": "observer time and frequency",
            "frequency_hz": "",
            "observer_time": "",
            "value_name": "capability entry",
            "value": "",
            "asgard_ratio": "",
            "wall_time": "",
            "boundary": "formal ratios available only for explicitly matched subsets",
        },
        "physical_channels": "FS synchrotron",
        "geometry_media": "structured jets and off-axis observers",
        "exported_states": "flux products; not ASGARD radius-ordered local state",
        "benchmark_evidence": "installed package; FigB1 broad-band FS rows; FigB2 Gaussian geometry subset",
        "explicit_boundaries": "not used for ASGARD RS, SSC, hadronic, or pair state claims",
    },
    {
        **{
            "scenario": "cross_code_capability",
            "code": "jetsimpy",
            "code_version": codes["jetsimpy"]["version"],
            "source": codes["jetsimpy"]["source"],
            "physical_assumptions": "public hydrodynamic afterglow model",
            "coordinate": "observer time and frequency",
            "frequency_hz": "",
            "observer_time": "",
            "value_name": "capability entry",
            "value": "",
            "asgard_ratio": "",
            "wall_time": "",
            "boundary": "formal ratios available only for explicitly matched subsets",
        },
        "physical_channels": "FS synchrotron",
        "geometry_media": "structured jets and off-axis observers",
        "exported_states": "flux products; not ASGARD radius-ordered local state",
        "benchmark_evidence": "installed source commit; FigB1 broad-band FS rows; FigB2 Gaussian geometry subset",
        "explicit_boundaries": "reverse shock and energy injection unsupported in public branch",
    },
    {
        **{
            "scenario": "cross_code_capability",
            "code": "VegasAfterglow",
            "code_version": codes["VegasAfterglow"]["version"],
            "source": codes["VegasAfterglow"]["source"],
            "physical_assumptions": "public GRB afterglow model",
            "coordinate": "observer time and frequency",
            "frequency_hz": "",
            "observer_time": "",
            "value_name": "capability entry",
            "value": "",
            "asgard_ratio": "",
            "wall_time": "",
            "boundary": "formal ratios available only for explicitly matched subsets",
        },
        "physical_channels": "public afterglow channels exposed by package API",
        "geometry_media": "public geometry and media controls exposed by package API",
        "exported_states": "flux products; not ASGARD radius-ordered local state",
        "benchmark_evidence": "installed package; FigB1 broad-band FS rows; FigB2 RS+SSC and Gaussian geometry subsets",
        "explicit_boundaries": "not used as an ASGARD full-dynamics or RS light-curve target",
    },
    {
        **{
            "scenario": "cross_code_capability",
            "code": "PyBlastAfterglowMag",
            "code_version": codes["PyBlastAfterglowMag"]["version"],
            "source": codes["PyBlastAfterglowMag"]["source"],
            "physical_assumptions": "public afterglow package interface",
            "coordinate": "observer time and frequency",
            "frequency_hz": "",
            "observer_time": "",
            "value_name": "capability entry",
            "value": "",
            "asgard_ratio": "",
            "wall_time": "",
            "boundary": "formal ratios available only for explicitly matched subsets",
        },
        "physical_channels": "public afterglow channels exposed by package API",
        "geometry_media": "public geometry and media controls exposed by package API",
        "exported_states": "package products; not ASGARD radius-ordered local state",
        "benchmark_evidence": "installed source commit; repo-local pba.out; FigB1 broad-band FS rows; FigB2 Gaussian geometry subset",
        "explicit_boundaries": "not used for ASGARD hadronic, pair, or chi-resolved state claims",
    },
]


write_csv("figB1_cross_code_fs.csv", fs_rows)
write_csv("figB2_rs_ssc_geometry.csv", fig2_rows)
write_csv("figB3_asgard_complex_state.csv", fig3_rows, fig3_fields)
write_csv("table_cross_code_capabilities.csv", cap_rows, table_fields)
