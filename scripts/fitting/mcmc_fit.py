import emcee
import corner
from schwimmbad import MPIPool
import numpy as np
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_core.asgard_fit import compile_problem, eval_loglike
from asgard_core.asgard_config import FitConfig


COMPILED_PROBLEM = None
BASE_CONFIG = FitConfig(
    num_threads=20,
    num_r=100,
    num_theta=300,
    num_phi=50,
    num_tobs=128,
    t_obs_min_log10=2.0,
    t_obs_max_log10=8.0,
    z=1.0,
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


def lnprior(theta):  # the log-prior probability function(system error): theta is the parameter spaces
    (E_iso, p, Eta_0, dNe, OpeningAngle_jet, Epsilon_e, Epsilon_b, f_e) = theta
    
    if (50. < E_iso < 55. and 
       2.0 < p < 2.9 and
       1.5 < Eta_0 < 3. and
       -4. < dNe < 3. and
       -3. < OpeningAngle_jet < 0. and
       -4 < Epsilon_e < -0.1 and
       -8. < Epsilon_b < -0.1 and
       -3. < f_e < 0. and
       Epsilon_b < Epsilon_e):

        return 0.0
    return -np.inf

def lnlike(theta):
    (E_iso, p, Eta_0, dNe, OpeningAngle_jet, Epsilon_e, Epsilon_b, f_e) = theta
    config = _config_from_theta(theta)
    if config is None:
        return -1e308
    return eval_loglike(_get_compiled_problem(), config)


def _config_from_theta(theta):
    E_iso, p, Eta_0, dNe, OpeningAngle_jet, Epsilon_e, Epsilon_b, f_e = theta
    if Epsilon_e < Epsilon_b or Eta_0 < -OpeningAngle_jet:
        return None
    config = FitConfig(**vars(BASE_CONFIG))
    config.e_iso = 10.0**E_iso
    config.p = float(p)
    config.eta_0 = 10.0**Eta_0
    config.d_ne = 10.0**dNe
    config.opening_angle_jet = 10.0**OpeningAngle_jet
    config.epsilon_e = 10.0**Epsilon_e
    config.epsilon_b = 10.0**Epsilon_b
    config.f_e = 10.0**f_e
    config.theta_v = 0.0
    return config


def _get_compiled_problem():
    global COMPILED_PROBLEM
    if COMPILED_PROBLEM is None:
        COMPILED_PROBLEM = compile_problem(BASE_CONFIG)
    return COMPILED_PROBLEM


def lnprob(theta):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta)

    
def main() -> None:
    ndim, nwalkers, nsteps = 8, 50, 4000
    filename = "chain_fe_z=1.h5"
    backend = emcee.backends.HDFBackend(filename)
    try:
        print("Initial size: {0}".format(backend.iteration))
    except Exception:
        print("new backend created.")

    with MPIPool() as pool:
        sampler = emcee.EnsembleSampler(
            nwalkers,
            ndim,
            lnprob,
            backend=backend,
            moves=[
                (emcee.moves.DEMove(), 0.7),
                (emcee.moves.DESnookerMove(), 0.2),
                (emcee.moves.StretchMove(), 0.1),
            ],
            pool=pool,
        )

        print("Running MCMC...")

        sampler.run_mcmc(None, 0 * (nsteps - backend.iteration), progress=True)

    print("Done.")

    parameters = [
        "$log_{10}E_{\\rm k,iso}$",
        "$p$",
        "$log_{10}\\Gamma$",
        "$log_{10}n_0$",
        "$log_{10}\\theta_{\\rm j}$",
        "$log_{10}\\epsilon_{e}$",
        "$log_{10}\\epsilon_{B}$",
        "$log_{10}f_{e}$",
    ]

    burnin = np.int32(0.7 * (nwalkers * nsteps))
    samples = sampler.flatchain[burnin:, ]

    E, p, Gamma, n0, thetaj, ee, eb, fe = map(
        lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
        zip(*np.percentile(samples, [16, 50, 84], axis=0)),
    )

    results = [E, p, Gamma, n0, thetaj, ee, eb, fe]

    np.savetxt("readm.txt", results, fmt="%s")
    # Plot 1 sigma error onto the data.
    fig = corner.corner(
        samples,
        labels=parameters,
        truths=[E[0], p[0], Gamma[0], n0[0], thetaj[0], ee[0], eb[0], fe[0]],
        quantiles=[0.16, 0.5, 0.84],
        plot_datapoints=False,
        show_titles=True,
        title_fmt=".2f",
        smooth=True,
        smooth1d=True,
        truth_color="darkorange",
        bins=30,
        range=[0.99] * 8,
        title_kwargs={"fontsize": 24},
        label_kwargs={"fontsize": 18},
    )

    fig.savefig("Ryde-triangle.pdf", dpi=300, bbox_inches="tight")

    print("Mean acceptance fraction:", np.mean(sampler.acceptance_fraction))

    try:
        print("Autocorrelation time:", sampler.acor)
    except ImportError:
        print("You can install acor: http://github.com/dfm/acor")


if __name__ == "__main__":
    main()
