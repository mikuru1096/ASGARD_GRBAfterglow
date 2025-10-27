import numpy as np
#from HE_pipline.My_constants import constants
from os import path
import matplotlib.pyplot as plt


def plotdata(data="XRT.txt"):
    d = path.abspath(__file__)
    parent_path = path.dirname(d)
    a = path.join(parent_path, "data/", data)
    table = np.loadtxt(a)
    time = table[:, 0]

    err_t_up = table[:, 1]
    err_t_down = -table[:, 2]
    flux = table[:, 3]
    flux_up = table[:, 4]
    flux_down = -table[:, 5]
    
    fig = plt.figure()
    axes = fig.add_axes([0.15, 0.1, 0.7, 0.8])
    line = axes.errorbar(time, flux, xerr = [err_t_down, err_t_up], yerr = [flux_down, flux_up], fmt='r.')
    
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
    
