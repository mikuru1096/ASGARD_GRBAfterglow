from extinction import fitzpatrick99 as f99
import numpy as np

def opt_extinction(mag_data,mag_err,frequency,Rv,Ebv,zeropointflux):
    
    wave = np.array([2.997e10/frequency*1e8])
    Av = Rv * Ebv
    
    mag_data_deredden = mag_data - f99(wave, Av, Rv)
    flux_data_deredden = 10**(0.4*(zeropointflux-mag_data_deredden))
    flux_data_err = 0.4*np.log(10.0)*flux_data_deredden*mag_err
    
    return flux_data_deredden, flux_data_err
