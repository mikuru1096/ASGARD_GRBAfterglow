import emcee
import corner
from matplotlib.pyplot import cm, savefig, show
import sys
from multiprocessing import Pool
from schwimmbad import MPIPool
import numpy as np
import matplotlib.pyplot as plt
from src import *
from mergered import fit
from astropy.cosmology import FlatLambdaCDM
from astropy import units
from matplotlib import animation
import gc
import os


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
    
    E_iso = 10 ** E_iso
    Eta_0 = 10 ** Eta_0
    dNe = 10 ** dNe
    OpeningAngle_jet = 10 ** OpeningAngle_jet
    Epsilon_e = 10 ** Epsilon_e
    Epsilon_b = 10 ** Epsilon_b
    f_e = 10 ** f_e
    
    Num_threads = 20
    
    parameters={'f_e': f_e,
                'E_iso': E_iso,
                'dNe': dNe,
                'Eta_0': Eta_0,
                'p': p,
                'OpeningAngle_jet': OpeningAngle_jet,
                'Epsilon_e': Epsilon_e,
                'Epsilon_b': Epsilon_b,
                'A_star': -1,
                'R0': 1e9,

                'z': 1,
                'Ebv': 0,
                'theta_v': 0,
                'Num_phi': 1,
                'Num_threads': Num_threads,
                'Num_gam_e': 101,
                'Num_R': 100,
                'weno5': False,
                'reverse': False,
                'index_Y': 2,   # 1 full numerical  2 Nakar  3 Fan

                'plot_LC': False,
               }

    redchi = fit(**parameters)

    if np.isnan(redchi):
        redchi = np.inf

    return -0.5 * redchi


def lnprob(theta):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta)

    
ndim, nwalkers, nsteps = 8, 50, 4000
filename = "chain_fe_z=1.h5"
backend = emcee.backends.HDFBackend(filename)
try:
    print("Initial size: {0}".format(backend.iteration))
except:
    print("new backend created.")

with MPIPool() as pool:

    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, backend=backend,
                                    moves=[(emcee.moves.DEMove(), 0.7), 
                                    (emcee.moves.DESnookerMove(), 0.2),
                                    (emcee.moves.StretchMove(), 0.1), ],
                                    pool=pool)

    print("Running MCMC...")

    initial = [(np.random.uniform(50., 55.5),  #Eiso
                np.random.uniform(2., 2.9),  #p
                np.random.uniform(1.5, 3),  #Gamma
                np.random.uniform(-4, 3),     #n0
                np.random.uniform(-3., 0.),  #theta_j
                np.random.uniform(-3, -0.1),  #epe
                np.random.uniform(-7., -0.1),  #epB
                np.random.uniform(-3, 0),    #fe
                ) for i in range(nwalkers)]  
                 
    sampler.run_mcmc(None, 0*(nsteps-backend.iteration), progress=True)


print("Done.")


parameters = ["$log_{10}E_{\\rm k,iso}$", 
              "$p$", 
              "$log_{10}\\Gamma$", 
              "$log_{10}n_0$",
              "$log_{10}\\theta_{\\rm j}$", 
              "$log_{10}\\epsilon_{e}$", 
              "$log_{10}\\epsilon_{B}$",
              "$log_{10}f_{e}$"
              ]
              
names = ["log_{10}E_{\\rm k,iso}", 
              "p", 
              "log_{10}\\Gamma", 
              "log_{10}n_0",
              "log_{10}\\theta_{\\rm j}", 
              "log_{10}\\epsilon_{e}", 
              "log_{10}\\epsilon_{B}",
              "log_{10}f_{e}"
              ]

burnin = np.int32(0.7 * (nwalkers * nsteps))
samples = sampler.flatchain[burnin:, ]

E, p, Gamma, n0, thetaj, ee, eb, fe = map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))

results = [E, p, Gamma, n0, thetaj, ee, eb, fe]

np.savetxt("readm.txt", results, fmt='%s')

'''
samples = sampler.get_chain()
flat_samples = sampler.get_chain(discard=200, thin=100, flat=True)
from getdist import plots, MCSamples
import getdist

samples1=MCSamples(samples=flat_samples[:,], names=names, labels=names)
g=plots.get_subplot_plotter(subplot_size=3)
g.settings.figure_legend_frame = False
g.settings.alpha_filled_add = 0.2
g.settings.title_limit_fontsize = 28
g.settings.title_limit_labels = True
g.settings.axes_labelsize = 30
g.settings.axes_fontsize = 22
g.settings.plot_meanlikes = False
g.settings.tight_layout  = True

g.triangle_plot(roots=samples1,
                labels=parameters,
                filled=True,
    #            contour_colors='g',
                title_limit=1,
                shaded=False)
                
plt.savefig("Ryde-triangle.pdf", dpi=300)

'''
    # Plot 1 sigma error onto the data.
fig = corner.corner(samples,
                        labels=parameters,
                        truths=[E[0], p[0], Gamma[0], n0[0], thetaj[0], ee[0], eb[0], fe[0]],
                        quantiles=[0.16, 0.5, 0.84], 
                        plot_datapoints=False, 
                        show_titles=True,
                        title_fmt=".2f", 
                        smooth=True, 
                        smooth1d=True,
                        truth_color='darkorange', 
                        bins=30,
                        range=[0.99]*8,
                        title_kwargs={"fontsize":24},
                        label_kwargs={"fontsize":18}
                        )  # 1sigma error

fig.savefig("Ryde-triangle.pdf", dpi=300, bbox_inches='tight') # to show full title of the most top subfig.
#show()

    # --
print("Mean acceptance fraction:", np.mean(sampler.acceptance_fraction))

    # If you have installed acor (http://github.com/dfm/acor), you can estimate
    # the autocorrelation time for the chain. The autocorrelation time is also
    # a vector with 10 entries (one for each dimension of parameter space).
try:
    print("Autocorrelation time:", sampler.acor)
except ImportError:
    print("You can install acor: http://github.com/dfm/acor")
