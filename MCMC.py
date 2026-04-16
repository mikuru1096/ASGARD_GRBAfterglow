import emcee
import corner
from schwimmbad import MPIPool
import numpy as np

from asgard_inference import (
    build_log_inference_config,
    compile_inference_problem,
    evaluate_compiled_loglike,
    populate_log_inference_config,
)


COMPILED_PROBLEM = None
COMPILED_CONFIG = None


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
    problem, config = _get_compiled_context()
    populate_log_inference_config(
        config,
        f_e_log10=f_e,
        e_iso_log10=E_iso,
        d_ne_log10=dNe,
        eta_0_log10=Eta_0,
        p=p,
        opening_angle_jet_log10=OpeningAngle_jet,
        epsilon_e_log10=Epsilon_e,
        epsilon_b_log10=Epsilon_b,
        z=1.0,
        theta_v=0.0,
        num_threads=20,
        num_r=100,
    )
    return evaluate_compiled_loglike(problem, config)


def _get_compiled_context():
    global COMPILED_PROBLEM, COMPILED_CONFIG
    if COMPILED_PROBLEM is None or COMPILED_CONFIG is None:
        COMPILED_CONFIG = build_log_inference_config(
            f_e_log10=-1.0,
            e_iso_log10=53.0,
            d_ne_log10=-1.0,
            eta_0_log10=2.0,
            p=2.5,
            opening_angle_jet_log10=-1.0,
            epsilon_e_log10=-1.0,
            epsilon_b_log10=-3.0,
            z=1.0,
            theta_v=0.0,
            num_threads=20,
            num_r=100,
        )
        COMPILED_PROBLEM = compile_inference_problem(COMPILED_CONFIG)
    return COMPILED_PROBLEM, COMPILED_CONFIG


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
