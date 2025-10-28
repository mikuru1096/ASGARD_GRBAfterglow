import matplotlib.pyplot as plt
import numpy as np
import sys
from src import *
from extinc import opt_extinction
from astropy.cosmology import FlatLambdaCDM
from astropy import units
from scipy.interpolate import interp1d
from cal_chi2_spectrum import cal_chi2_spectrum as chispec
from cal_chi2_light_curve import cal_chi2_light_curve as chilc

import gc
import os
    
def plot_data(axes, Tobs, flux, color, label, scaling_factor=1):
    """Utility function to plot data on log-log scale."""
    line = axes.loglog(Tobs, flux * scaling_factor, color=color, linewidth=2.5, label=label, alpha=1)
    return line

def plot_syn_lc(Tobs, Flux_optr, Flux_optz, Flux_opti, Flux_optg, Flux_9GHz, Flux_5GHz, Flux_3GHz, Flux_XRT, Flux_GeV, Flux_TeV):

    fig, axes = plt.subplots(figsize=(10, 8))

    # Plot light curves
    plot_data(axes, Tobs, Flux_XRT, 'k', "$0.5-10$ keV", scaling_factor=1e0)
    
    plot_data(axes, Tobs, Flux_optr, '#FFD700', "r $\\times$ 0.3e-16", scaling_factor=0.3 * 1e-16)
    
    plot_data(axes, Tobs, Flux_optz, '#FF8C00', "z $\\times$ 0.1e-16", scaling_factor=0.1 * 1e-16)
    
    plot_data(axes, Tobs, Flux_opti, '#FF4500', "i $\\times$ 0.03e-16", scaling_factor=0.03 * 1e-16)
    
    plot_data(axes, Tobs, Flux_optg, '#FF0000', "g $\\times$ 0.01e-16", scaling_factor=0.01 * 1e-16)
    
    plot_data(axes, Tobs, Flux_9GHz, '#FF1493', "9GHz $\\times$ 1e-19", scaling_factor=1 * 1e-19)
    plot_data(axes, Tobs, Flux_5GHz, '#8A2BE2', "5.5GHz $\\times$ 0.3e-19", scaling_factor=0.3 * 1e-19)
    plot_data(axes, Tobs, Flux_3GHz, '#4B0082', "3GHz $\\times$   0.1e-19", scaling_factor=0.1 * 1e-19)
    
    plot_data(axes, Tobs, Flux_GeV, '#9400D3', "1e24 $\\times$   1e4", scaling_factor=1 * 1e4)
    
    plot_data(axes, Tobs, Flux_TeV, '#FF00FF', "1e27 $\\times$   1e5", scaling_factor=1 * 1e5)
    
    # Set limits
    axes.set_ylim([1.0e-22, 1.0e-4])
    
    axes.yaxis.set_tick_params(labelsize=14)
    axes.xaxis.set_tick_params(labelsize=14)
    
    axes.set_xlim([1.0e2, 1.0e8])

    # Set labels and title
    axes.set_xlabel('Time [s]', fontsize=14)
    axes.set_ylabel('Flux erg/cm$^2$/s or $\\mu$Jy', fontsize=14)
    axes.set_title('LCs', fontsize=18)

    axes.tick_params(
    which='both',
    direction='in',
    bottom=True,
    top=True,
    left=True,
    right=True
    )

    #plt.xticks(size=16)
    #plt.yticks(size=16)

    # Creating the legend
    handles, labels = axes.get_legend_handles_labels()
    axes.legend(handles, labels, bbox_to_anchor=(1, 1), loc="upper right", ncol=1, 
                borderaxespad=0.5, fancybox=True, fontsize=10)

    plt.xscale('log')
    plt.yscale('log')

    plt.show()

    fig.savefig("Radiation_Lightcurves.pdf", dpi=300)
    plt.clf()
    plt.close()
    gc.collect()

    return 0


def fit(*args,**kwargs):

    Eta_0 = kwargs.get('Eta_0')  # Initial Lorentz factor
    Epsilon_e = kwargs.get('Epsilon_e')  # Electron energy fraction of jet
    Epsilon_b = kwargs.get('Epsilon_b')  # Magnetic energy fraction of jet
    p = kwargs.get('p')  # Spectral index of electrons
    theta_v = kwargs.get('theta_v', 0)  # Viewing angle of observer
    OpeningAngle_jet = kwargs.get('OpeningAngle_jet')  # Openning angle of jet
    f_e = kwargs.get('f_e',1.0)  # Non-thermal fraction of electrons
    dNe = kwargs.get('dNe')  # Uniform ambient medium density
    A_star = kwargs.get('A_star', -1)  # Steller wind parameter, < 0 for Uniform ambient medium. If enabled dNe & A_star > 0 simultaneously, it enters the Wind-to-ISM density jump scenario. Set dNe be a very small value 1e-15 for Wind-dominated scenario.
    E_iso = kwargs.get('E_iso')  # Isotropic kinetic energy of jet
    z = kwargs.get('z')  # Redshift
    
    # Openmp threads number, used in Fortran subroutine. 
    # Maximum value same as your CPU cores.
    Num_threads = kwargs.get('Num_threads', 8)
    
    # Number of grid points of \theta dirction. 
    # Warning!! Minimum bins 250, lower value lead to morphological deviation of spec and lc.
    Num_theta = kwargs.get('Num_theta', 300)
    
    # Number of grid points of \phi dirction. 
    # Warning!! When theta_v=0 set it to one, otherwize set it to few of tens.
    Num_phi = kwargs.get('Num_phi', 1)
    
    # Number of grid points of electron energy.
    Num_gam_e = kwargs.get('Num_gam_e', 101)  
    
    # Number of grid points of INTRINSIC photon energy.
    Num_nu = kwargs.get('Num_nu', 200)  
    
    # Number of grid points of jet dynamic time (radius).
    Num_R = kwargs.get('Num_R', 300)  
    
    # Flag for the method of calculating the Compton cooling factor. 
    # 1 full numerical  2 Nakar 2007  3 Fan 2003
    index_Y = kwargs.get('index_Y', 2) 
    
    # Flag for the numerical method of calculating the Synchrotron Radiation. 
    # 1 trapezoidal O(h^2) 2 composite Simpson O(h^4)
    index_syn_intger = kwargs.get('index_syn_intger', 2) 
    
    # Flag for the numerical method of calculating the jet dynamics. 
    # 1 Huang  2 Pe'er  3 Zhang/Nava
    index_dyn = kwargs.get('index_dyn', 3) 
    
    # Flag for the numerical method of electron spectrum calculation. 
    # WENO5 (Weighted Essentially Non-Oscillatory scheme, 5th order) is a explicit 5th-order scheme. (not used)
    weno5 = kwargs.get('weno5', False)  
    
    plot_LC = kwargs.get('plot_LC', False)  # Flag for light curve plot
    
    R0 = kwargs.get('R0', 1e9)  # If you need an Uniform-to-Wind ambient medium, set this parameter as the transition radius.
    
    '''
    Start Calculation!
    '''
    
    # The observation frequencies
    # Name of the data file for fitting. It must correspond one-to-one with the files in the data_light_curve directory. 
    # Since optical data requires special processing, the filename must begin with "opt".
    band = ['xrt', 'optr', 'optz', 'opti', 'optg', '9GHz', '5.5GHz', '3GHz']
    
    # Set Calculation Frequencies
    Num_XRT = 8
    opt = [4.63e14, 3.39e14, 4.01e14, 6.42e14]
    radio = [9e9, 5.5e9, 3e9]
    XRT_low = 0.5 * constants.para_kev2hz
    XRT_high = 10 * constants.para_kev2hz
    XRT = np.logspace(np.log10(XRT_low), np.log10(XRT_high), Num_XRT)
    GeV_Frequency = [1e24]
    TeV_Frequency = [1e27]
    
    Frequencies = np.concatenate((XRT, opt, radio, GeV_Frequency, TeV_Frequency), axis=0)
    Num_nu_obs = len(Frequencies)

    # parameters of afterglow
    
    T_log10_duration = 8.0  # Time calculated to 1e8 seconds after triggering.
    m_0 = 1.0e12  # Initial mass of jet shell
    u_0 = 1.0e13  # Initial internal energy of jet shell
    r_0 = 1.0e14  # Initial external-forward shock radius of jet shell

    cosmo = FlatLambdaCDM(H0=67.8, Om0=0.308)
    dL = cosmo.luminosity_distance(z=z).to(units.cm).value  # Luminosity distance

    # DO NOT CHANGE!!
    Boundary = np.array([Eta_0, m_0, u_0, r_0, Epsilon_e, Epsilon_b, p, z, OpeningAngle_jet, theta_v, dNe, A_star, dL, E_iso, T_log10_duration, f_e, R0])
    
    # Time serial for Observer
    Num_Tobs = 200
    Tobs = np.logspace(2.0, 8.0, Num_Tobs)
    Fluxes_final = np.zeros([len(band), Num_Tobs])

    # DO NOT CHANGE unless you know what happens.
    V_seed_min = np.log10(1.0e-8 * constants.para_ev2hz)
    V_seed_max = np.log10(1.0e3 * constants.para_tev2hz)
    V_seed = np.logspace(V_seed_min, V_seed_max, Num_nu)

    # Through certain settings, it can be used to generate energy spectra.
#    Num_V_obs = 80
#    V_obs_min = np.log10(1.0e-6 * constants.para_ev2hz)
#    V_obs_max = np.log10(1.0e-3 * constants.para_tev2hz)
#    V_obs = np.logspace(V_obs_min, V_obs_max, Num_V_obs)

    # Calculate the dynamics, generating the sequence of redshifted time, Lorentz factor, and jet radius in the Rest Frame of burst.
    R_Tobs, R_Gamma, R = Dynamics.dynamics_forward(Boundary, Num_R, index_dyn)
    
    # Calculate the electron spectrum, generating the electron energy, electron number spectrum, 
    # intrinsic synchrotron radiation luminosity, and synchrotron photon number density.
    # Uses a fully implicit scheme by default, with first-order accuracy.
    gam_e, dN_gam_e_fs, L_syn_spec_f, seed_syn_f = Electron.fs_electron_fullhide(Boundary, R_Tobs, R_Gamma, R, V_seed, Num_gam_e, index_Y, index_syn_intger, Num_threads)
    
    # Call the WENO5 scheme for comparison 
#    gam_e1, dN_gam_e_fs1, L_syn_spec_f1, seed_syn_f1 = Electron.fs_electron_weno5(Boundary, R_Tobs, R_Gamma, R, V_seed, Num_gam_e, index_Y, Num_threads)
#    plt.loglog(gam_e,np.abs(dN_gam_e_fs[:,270]-dN_gam_e_fs1[:,270])/dN_gam_e_fs1[:,270])
#    plt.show()
    
    # Calculate the synchrotron self-Compton radiation, and generate the intrinsic luminosity and photon number density.
    L_SSC_spec_f,seed_ssc_f = Radiation.ssc_spec(R, gam_e, dN_gam_e_fs, V_seed, seed_syn_f, Num_threads)
    
    L_tot = L_syn_spec_f + L_SSC_spec_f
    
    
    # Calculate the absorption of Photon-Photon Annihilation, generating the extinction fraction.
    absorption = Radiation.annihilation(R_Gamma, R, V_seed, seed_syn_f, seed_ssc_f, Num_threads)

    F_tot_abs = L_tot*absorption/(4.0*np.pi*dL**2)*(1.0+z)

    # Calculate the equal arrival time surface effect and the Doppler effect, and output the observational results.
    # The top-hat jet uses this module, while the structured jet uses a separate module.
    Fluxes = Interpolation.sed_interpolation(Boundary, R_Tobs, R_Gamma, R, F_tot_abs, V_seed, Frequencies, Tobs, Num_theta, Num_phi, Num_threads)
    
    # If needs EBL absorption, set it by yourself.
#    ebl_absorption_z = Radiation.cal_ebl(z, V_obs, model="Saldana21.out")


#    f = interp1d(V_obs, F_tot_obs.T, kind='cubic') # 

#    Fluxes=f(Frequencies).T

    # units in erg/cm^2/s
    Flux_XRT = np.trapz(Fluxes[:Num_XRT].T,Frequencies[:Num_XRT])
    
    # units in \muJy
    opt_fluxes = Fluxes[Num_XRT:Num_XRT+4, :] * 1e29
    radio_fluxes = Fluxes[Num_XRT+4:Num_XRT+7, :] * 1e29
    
    # units in erg/cm^2/s
    GeV_fluxes = Fluxes[Num_XRT+7, :] * GeV_Frequency
    TeV_fluxes = Fluxes[Num_XRT+8, :] * TeV_Frequency
    
    plot_syn_lc(Tobs, *opt_fluxes, *radio_fluxes, Flux_XRT, GeV_fluxes, TeV_fluxes)

    gc.collect()
    
    return 0
