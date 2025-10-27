from mergered import fit
import numpy as np

E_iso = 10 ** 53
Eta_0 = 10 ** 2
n_0 = 10 ** -2
A_star = -1
theta_j = 10 ** -1
theta_v = 0 * theta_j
p = 2.2
Epsilon_e = 10 ** -1
Epsilon_b = 10 ** -3
f_e = 10 ** -1
z = 0.4

params={
        'Num_threads': 8,
        
        'index_dyn': 3,          #   1: Huang                2: Pe'er                      3: Zhang/Nava
        'index_Y': 2,            #   1: full numerical       2: Nakar                      3: Fan
        'index_syn_intger': 2,   #   1: trapezoidal O(h^2)   2: composite Simpson O(h^4)
        
        'Num_gam_e': 101,
        'Num_nu': 101,
        'Num_R': 300,
        'Num_theta': 300,
        'Num_phi': 1,
        
        'z': z,
        'Eta_0': Eta_0,
        'Epsilon_e': Epsilon_e,
        'Epsilon_b': Epsilon_b,
        'p': p,
        'OpeningAngle_jet': theta_j,
        'f_e': f_e,
        'E_iso': E_iso,
        'dNe': n_0,
        'A_star': A_star,
        'R0': 1e9,
        
        'Ebv': 2.93,
        'Lyman_Ar': 0,
        'f_sys': -1,
        'theta_v': theta_v,
        
        'weno5': False,
        'reverse': False,

        'plot_LC': True,
       }

fit(**params)







