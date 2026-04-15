from __future__ import annotations

from asgard_models import FitConfig


LEGACY_CONFIG_MAPPING = {
    "Num_threads": "num_threads",
    "index_dyn": "index_dyn",
    "index_Y": "index_y",
    "index_syn_intger": "index_syn_integr",
    "Num_gam_e": "num_gam_e",
    "Num_nu": "num_nu",
    "Num_R": "num_r",
    "Num_theta": "num_theta",
    "Num_phi": "num_phi",
    "z": "z",
    "Eta_0": "eta_0",
    "Epsilon_e": "epsilon_e",
    "Epsilon_b": "epsilon_b",
    "p": "p",
    "OpeningAngle_jet": "opening_angle_jet",
    "theta_v": "theta_v",
    "f_e": "f_e",
    "E_iso": "e_iso",
    "dNe": "d_ne",
    "A_star": "a_star",
    "R0": "r0",
    "Ebv": "ebv",
    "Rv": "rv",
    "Lyman_Ar": "lyman_ar",
    "f_sys": "f_sys",
    "weno5": "weno5",
    "reverse": "reverse",
    "plot_LC": "plot_lc",
    "E_inj_t1": "e_inj_t1",
    "E_inj_t2": "e_inj_t2",
    "L_inj_0": "l_inj_0",
    "q_inj": "q_inj",
    "R_tr": "r_tr",
    "f_jump": "f_jump",
    "f_wide": "f_wide",
}

LEGACY_IGNORE_KEYS = {"plot_syn_curve", "plot_spectrum", "do_plot_spec"}

LEGACY_REVERSE_MAPPING = {
    "Delta_t": "delta_t_s",
    "Epsilon_e_r": "epsilon_e",
    "Epsilon_b_r": "epsilon_b",
    "p_r": "p",
    "f_e_r": "f_e",
}


def legacy_kwargs_to_config(kwargs: dict) -> FitConfig:
    config = FitConfig()
    for key, value in kwargs.items():
        if key in LEGACY_IGNORE_KEYS:
            continue
        if key in LEGACY_REVERSE_MAPPING:
            setattr(config.reverse_shock, LEGACY_REVERSE_MAPPING[key], value)
            config.reverse_shock.enabled = True
            continue
        if key not in LEGACY_CONFIG_MAPPING:
            raise KeyError(f"Unsupported legacy key: {key}")
        setattr(config, LEGACY_CONFIG_MAPPING[key], value)
    return config
