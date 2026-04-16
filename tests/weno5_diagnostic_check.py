from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_paths import asgard_doc_path
from asgard_presets import build_baseline_config
from asgard_setup import build_simulation_setup
from asgard_runtime import solve_dynamics, solve_electron
import src.Electron.FS_electron_weno5 as weno5_module
from src import constants


OUTPUT_JSON = asgard_doc_path("weno5_diagnostic.json")


def _ambient_density(radius_cm: float, config) -> float:
    if config.a_star > 0.0:
        density = config.a_star * 3.0e35 / radius_cm**2
        if density <= config.d_ne / 4.0:
            density = config.d_ne
        if radius_cm < config.r0:
            density = config.a_star * 3.0e35 / config.r0**2
        return float(density)
    return float(
        config.d_ne
        * (1.0 + (config.f_jump - 1.0) * np.exp(-(np.log10(radius_cm) - np.log10(config.r_tr)) ** 2 / (2.0 * config.f_wide**2)))
    )


def _gamma_m(config, gamma_loc: float, gamma_max: float) -> float:
    temp_gam = config.epsilon_e / config.f_e * constants.para_m_p_div_m_e * (gamma_loc - 1.0)
    if config.p < 2.05 and config.p >= 2.0:
        return 0.05 / 1.05 * temp_gam + 1.0
    if config.p < 2.0 and config.p > 1.0:
        return ((2.0 - config.p) / (config.p - 1.0) * temp_gam * gamma_max ** (config.p - 2.0)) ** (1.0 / (config.p - 1.0)) + 1.0
    return (config.p - 2.0) / (config.p - 1.0) * temp_gam + 1.0


def _diagnose_case(num_points: int) -> dict:
    config = build_baseline_config(
        electron_solver="weno5",
        num_threads=1,
        num_gam_e=num_points,
        num_nu=num_points,
        num_r=160,
        include_forward_ssc=True,
    )
    setup = build_simulation_setup(config)
    dynamics = solve_dynamics(setup.boundary, config)
    electron = solve_electron(setup.boundary, dynamics, setup.seed_frequency_hz, config)

    max_dn = np.max(electron.d_n_gam_e, axis=0)
    max_ps = np.max(electron.l_syn_spec, axis=0)
    bad = np.where((~np.isfinite(max_dn)) | (~np.isfinite(max_ps)) | (max_dn > 1.0e50) | (max_ps > 1.0e50))[0]
    first_bad = None if bad.size == 0 else int(bad[0])

    shell = min(80, config.num_r - 2)
    radius_cm = float(dynamics.radius[shell - 1])
    gamma_loc = float(0.5 * (dynamics.r_gamma[shell] + dynamics.r_gamma[shell - 1]))
    density = _ambient_density(radius_cm, config)
    magnetic_field = 0.39 * np.sqrt(config.epsilon_b * density * gamma_loc * max(gamma_loc - 1.0, 0.0))
    gamma_max = 3.0 * constants.para_m_energy / np.sqrt(8.0 * magnetic_field * constants.para_e**3)
    gamma_m = _gamma_m(config, gamma_loc, gamma_max)
    gamma_c = 7.7e8 * (1.0 + config.z) / (gamma_loc * magnetic_field * magnetic_field * dynamics.r_tobs[shell])
    beta = np.sqrt(1.0 - gamma_loc**-2)
    gamma_max_arr = np.array(gamma_max, dtype=float)
    dEl = weno5_module.get_y.get_forward_cooling(
        config.index_y,
        config.epsilon_e,
        config.epsilon_b,
        config.p,
        magnetic_field,
        gamma_m,
        gamma_c,
        gamma_max_arr,
        radius_cm,
        gamma_loc,
        beta,
        density,
        config.num_threads,
        electron.gam_e,
        setup.seed_frequency_hz,
        electron.l_syn_spec[:, shell],
        electron.seed_syn[:, shell],
    )
    dEl1 = (dEl + 1.0 / radius_cm) / np.log(10.0)
    ddd = float(dynamics.radius[shell] - dynamics.radius[shell - 1])
    f_r = (1.35e-19) / beta / gamma_loc * magnetic_field**2 / np.pi
    ddr = 0.1 / (f_r * gamma_max + 1.333 / (dynamics.radius[shell] + dynamics.radius[shell - 1]))
    base_l1 = int(ddd / ddr) + 10
    ddr = ddd / base_l1
    dx = float(np.log10(electron.gam_e[1] / electron.gam_e[0]))
    cfl_geom = ddr / dx
    if np.all(np.isnan(dEl1)):
        max_speed = float("nan")
        cfl_eff = float("nan")
        required_l1 = None
    else:
        max_speed = float(np.nanmax(np.abs(dEl1)))
        cfl_eff = float(cfl_geom * max_speed)
        required_l1 = int(np.ceil(base_l1 * max(cfl_eff, 1.0))) if np.isfinite(cfl_eff) else None
        stabilized_l1 = int(np.ceil(base_l1 * max(cfl_eff / 0.8, 1.0))) if np.isfinite(cfl_eff) else None
    return {
        "num_points": num_points,
        "first_bad_shell": first_bad,
        "bands_flux_max": float(np.nanmax(electron.l_syn_spec)),
        "d_n_gam_e_max": float(np.nanmax(electron.d_n_gam_e)),
        "nu_a_max": float(np.nanmax(electron.nu_a)),
        "base_l1_shell80": base_l1,
        "cfl_geom_shell80": cfl_geom,
        "max_speed_shell80": max_speed,
        "cfl_eff_shell80": cfl_eff,
        "required_l1_cfl1_shell80": required_l1,
        "stabilized_l1_cfl08_shell80": stabilized_l1,
    }


def main() -> None:
    payload = {
        "cases": [_diagnose_case(n) for n in (61, 121, 241)],
    }
    OUTPUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_JSON.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(OUTPUT_JSON)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
