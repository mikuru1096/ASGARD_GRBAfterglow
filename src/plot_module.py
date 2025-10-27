import numpy as np
import matplotlib.pyplot as plt
from Dynamics import dynamics
from Seed_dyn import seed
from Constants import constants
from SSC_spec_dyn import ssc_spec
from Annihilation_dyn import annihilation
from SED_interpolartion import sed_interpolation
from Cal_ebl import cal_ebl
from astropy.cosmology import FlatLambdaCDM
from astropy import units

E_iso = 1e51
Eta_0 = 55
m_0 = 1e12
u_0 = 1e13
r_0 = 1e14
Epsilon_e = 1e-1
Epsilon_b = 1e-3
p = 2.2
z = 0.0785
OpeningAngle_jet = 0.1
theta_v = 0.0
dNe = 1.0
A_star = -1

cosmo = FlatLambdaCDM(H0=67.8, Om0=0.308)
dL = cosmo.luminosity_distance(z=z).to(units.cm).value

T_log10_duration = 7.0

Boundary = np.array([Eta_0, m_0, u_0, r_0, Epsilon_e, Epsilon_b, p, z, OpeningAngle_jet, theta_v, dNe, A_star, dL, E_iso, T_log10_duration])
Num_gam_e = 81

R_Tobs, R_Gamma, R, gam_e, dN_gam_e_fs = dynamics(Boundary, Num_gam_e)

V_obs_min = np.log10(1.0e-6 * constants.para_ev2hz)
V_obs_max = np.log10(1.0e3 * constants.para_tev2hz)
V_obs = np.logspace(V_obs_min, V_obs_max, 200)
Num_nu_obs = V_obs.size

V_seed_min = np.log10(1.0e-9 * constants.para_ev2hz)
V_seed_max = np.log10(1.0e6 * constants.para_tev2hz)
V_seed = np.logspace(V_seed_min, V_seed_max, 300)
Num_nu = V_seed.size

Num_R = 900
F_syn_spec,seed_syn = seed(Boundary, R_Gamma, R, gam_e, dN_gam_e_fs, V_seed, Num_nu, Num_R, Num_gam_e)

F_SSC_spec,seed_ssc = ssc_spec(R, gam_e, dN_gam_e_fs, V_seed, seed_syn, Num_nu, Num_R, Num_gam_e)

absorption = annihilation(R_Gamma, R, V_seed, seed_syn, seed_ssc, Num_nu, Num_R)

seed_tot = seed_syn + seed_ssc
F_tot = F_syn_spec + F_SSC_spec
F_tot = F_tot/(4.0*np.pi*dL**2)*(1.0+z)

Num_Tobs = 50
Tobs = np.logspace(-1,7,Num_Tobs)
Num_Theta = 300

F_tot_obs = sed_interpolation(Boundary, R_Tobs, R_Gamma, R, F_tot, V_seed, V_obs, Tobs, Num_nu, Num_nu_obs, Num_Tobs, Num_Theta, Num_R)

fig = plt.figure()
axes = fig.add_axes([0.15, 0.1, 0.7, 0.8])
line = axes.plot(V_seed,V_seed*F_tot[800,:]*absorption[800,:])
line2 = axes.plot(V_seed,V_seed*F_tot[800,:])

plt.xscale('log')
plt.yscale('log')

fig1 = plt.figure()
axes1 = fig1.add_axes([0.15, 0.1, 0.7, 0.8])
line1 = axes1.plot(R_Tobs,V_seed[100]*F_tot[:,100]*absorption[:,100])

plt.xscale('log')
plt.yscale('log')

ebl_absorption_z = cal_ebl(z, V_obs, model="Dominguez11.txt")

fig2 = plt.figure()
axes2 = fig2.add_axes([0.15, 0.1, 0.7, 0.8])
#line3 = axes2.plot(V_obs/ constants.para_ev2hz, V_obs*F_tot_obs[40,:]*ebl_absorption_z)
#line4 = axes2.plot(V_obs/ constants.para_ev2hz, V_obs*F_tot_obs[40,:])
line5 = axes2.plot(Tobs,V_obs[100]*F_tot_obs[:,100])
plt.xscale('log')
plt.yscale('log')
plt.show()



















