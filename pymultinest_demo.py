#!/usr/bin/env python
from __future__ import absolute_import, unicode_literals, print_function
from pymultinest.solve import solve
import os
from asgard_inference import build_inference_config, evaluate_fit_loglike

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
    
    config = build_inference_config(
        d_ne=n_0,
        a_star=A_star,
        z=z,
        ebv=Ebv,
        lyman_ar=Lyman_Ar,
        f_sys=f_sys,
        theta_v=theta_v,
        num_phi=50,
        num_threads=8,
        num_r=500,
        num_theta=300,
        eta_0=Eta_0,
        epsilon_e=Epsilon_e,
        epsilon_b=Epsilon_b,
        p=p,
        opening_angle_jet=theta_j,
        f_e=f_e,
        e_iso=E_iso,
    )

    return evaluate_fit_loglike(config)


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
