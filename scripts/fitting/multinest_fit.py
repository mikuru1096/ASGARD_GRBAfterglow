#!/usr/bin/env python
from __future__ import absolute_import, unicode_literals, print_function
from pymultinest.solve import solve
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_core.asgard_fit import compile_problem, eval_loglike
from asgard_core.asgard_models import FitConfig

import json


parameters = [
    "$log_{10}n_0$",
    "$log_{10}A_{\\star}$",
    "$log_{10}E_{\\rm k,iso}$",
    "$p$",
    "$log_{10}\\Gamma$",
    "$log_{10}\\epsilon_{e}$",
    "$log_{10}\\epsilon_{B}$",
    "$log_{10}\\theta_{j}$",
    "$\\theta_v/\\theta_j$",
    "$log_{10}\\xi_e$",
    "$E(B-V)$",
    "$Ly\\alpha(A_r)$",
]
n_params = len(parameters)
prefix = "chains/3-"
COMPILED_PROBLEM = None
BASE_CONFIG = FitConfig(
    num_threads=8,
    num_r=500,
    num_theta=300,
    num_phi=50,
    num_tobs=200,
    t_obs_min_log10=2.0,
    t_obs_max_log10=8.0,
    z=4.59,
    d_ne=1.0,
    a_star=1.0,
    eta_0=100.0,
    epsilon_e=1.0e-2,
    epsilon_b=1.0e-5,
    p=2.3,
    opening_angle_jet=0.1,
    f_e=0.1,
    e_iso=1.0e55,
    ebv=0.0,
    lyman_ar=0.0,
    f_sys=-1.0,
)


# probability function, taken from the eggbox problem.
def prior_transform(cube):
    # the argument, cube, consists of values from 0 to 1
    # we have to convert them to physical scales

    params = cube.copy()
    params[0] = -3  + cube[0]*6 #n_0
    params[1] = -2  + cube[1]*3 #A_star
    params[2] = 54.  + cube[2]*2 #E_iso
    params[3] = 2.01 + cube[3]*0.99 #p
    params[4] = 1.4  + cube[4]*2.1 #Gamma_0
    params[5] = -3   + cube[5]*2.99 #Epsilon_e
    params[6] = -8 + cube[6]*7.99 #Epsilon_B
    params[7] = -2.5 * cube[7] #log10(theta_j)
    params[8] =  0  + cube[8]*5 #theta_v/theta_j
    params[9] = -3 * cube[9] #f_e
    params[10] = 1 * cube[10] #Ebv
    params[11] = 8 * cube[11] #Lyman_Ar
 #   params[12] = -3 * cube[12] #relative system error
    return params
    
    
def myloglike(params):
    n_0, A_star, E_iso, p, Eta_0, Epsilon_e, Epsilon_b, theta_j, theta_v, f_e, Ebv, Lyman_Ar = params
    
    if Epsilon_e < Epsilon_b or Eta_0 < -theta_j:
        return -1e308

    E_iso = 10 ** E_iso
    Eta_0 = 10 ** Eta_0
    n_0 = 10 ** n_0
    A_star = 10 ** A_star
    theta_j = 10 ** theta_j
    theta_v = theta_v * theta_j
    Epsilon_e = 10 ** Epsilon_e
    Epsilon_b = 10 ** Epsilon_b
    f_e = 10 ** f_e
    f_sys = -1 #10 ** f_sys

    z = 4.59
    
    config = _config_from_params(
        n_0=n_0,
        A_star=A_star,
        E_iso=E_iso,
        p=p,
        Eta_0=Eta_0,
        Epsilon_e=Epsilon_e,
        Epsilon_b=Epsilon_b,
        theta_j=theta_j,
        theta_v=theta_v,
        f_e=f_e,
        Ebv=Ebv,
        Lyman_Ar=Lyman_Ar,
        f_sys=f_sys,
        z=z,
    )
    if config is None:
        return -1e308
    return eval_loglike(_get_compiled_problem(), config)


def _config_from_params(
    *,
    n_0: float,
    A_star: float,
    E_iso: float,
    p: float,
    Eta_0: float,
    Epsilon_e: float,
    Epsilon_b: float,
    theta_j: float,
    theta_v: float,
    f_e: float,
    Ebv: float,
    Lyman_Ar: float,
    f_sys: float,
    z: float,
):
    if Epsilon_e < Epsilon_b or Eta_0 < -theta_j:
        return None
    config = FitConfig(**vars(BASE_CONFIG))
    config.d_ne = 10.0**n_0
    config.a_star = 10.0**A_star
    config.e_iso = 10.0**E_iso
    config.p = float(p)
    config.eta_0 = 10.0**Eta_0
    config.epsilon_e = 10.0**Epsilon_e
    config.epsilon_b = 10.0**Epsilon_b
    config.opening_angle_jet = 10.0**theta_j
    config.theta_v = theta_v * config.opening_angle_jet
    config.f_e = 10.0**f_e
    config.ebv = float(Ebv)
    config.lyman_ar = float(Lyman_Ar)
    config.f_sys = float(f_sys)
    config.z = float(z)
    return config


def _get_compiled_problem():
    global COMPILED_PROBLEM
    if COMPILED_PROBLEM is None:
        COMPILED_PROBLEM = compile_problem(BASE_CONFIG)
    return COMPILED_PROBLEM


def main() -> None:
    try:
        os.mkdir("chains")
    except OSError:
        pass

    with open(f"{prefix}params.json", "w") as f:
        json.dump(parameters, f, indent=2)

    result = solve(
        LogLikelihood=myloglike,
        Prior=prior_transform,
        n_dims=n_params,
        outputfiles_basename=prefix,
        verbose=True,
    )

    print()
    print("evidence: %(logZ).1f +- %(logZerr).1f" % result)
    print()
    print("parameter values:")
    for name, col in zip(parameters, result["samples"].transpose()):
        print("%15s : %.3f +- %.3f" % (name, col.mean(), col.std()))


if __name__ == "__main__":
    main()
