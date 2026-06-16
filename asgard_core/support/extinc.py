from extinction import fitzpatrick99 as f99
import numpy as np
from .extinction_cur import get_abs

def _frequency_to_wavelength_micron(frequency):
    return np.array([2.997e10 / frequency * 1e4], dtype=float)


def _flux_from_mag(mag_data, mag_err, zeropointflux):
    flux_data_deredden = 10 ** (0.4 * (zeropointflux - mag_data))
    flux_data_err = 0.4 * np.log(10.0) * flux_data_deredden * mag_err
    return flux_data_deredden, flux_data_err


def opt_extinction(
    mag_data,
    mag_err,
    frequency,
    Rv,
    Ebv,
    zeropointflux,
    redshift=None,
    lyman_ar=0.0,
    law="fitzpatrick99",
):
    Av = Rv * Ebv

    if law == "fitzpatrick99":
        wave = np.array([2.997e10 / frequency * 1e8], dtype=float)
        mag_data_deredden = mag_data - f99(wave, Av, Rv)
        if redshift is not None and lyman_ar != 0.0:
            wave_in_mu_m = _frequency_to_wavelength_micron(frequency)
            if np.any((wave_in_mu_m > 0.6) & (wave_in_mu_m < 0.68)):
                mag_data_deredden = mag_data_deredden - lyman_ar
        return _flux_from_mag(mag_data_deredden, mag_err, zeropointflux)

    if law == "zhou_smc":
        return opt_extinction_zhou(mag_data, mag_err, frequency, Rv, Ebv, zeropointflux, redshift, lyman_ar)

    if law == "pei92":
        return opt_extinction_pei92(mag_data, mag_err, frequency, "SMC", Rv, Ebv, zeropointflux, redshift)

    raise ValueError(f"Unsupported extinction law: {law}")
    
    
def opt_extinction_zhou(mag_data,mag_err,frequency,Rv,Ebv,zeropointflux,redshift,Lyman_Ar):
    
    wave_in_mu_m = _frequency_to_wavelength_micron(frequency)
    Av = Rv * Ebv
    
    mag_data_deredden = mag_data - get_abs(wave_in_mu_m, Av, redshift, Lyman_Ar)
    return _flux_from_mag(mag_data_deredden, mag_err, zeropointflux)
    
def opt_extinction_pei92(mag_data,mag_err,frequency,model,Rv,Ebv,zeropointflux,redshift):
    
    wave_in_mu_m = _frequency_to_wavelength_micron(frequency)
    wave_in_mu_m_redshift = wave_in_mu_m / (1.0 + redshift)
    
    mag_data_deredden = mag_data - pei92(wave_in_mu_m_redshift, Rv, Ebv, model) - pei92(wave_in_mu_m, 3.08, 0.29, 'MW')
    return _flux_from_mag(mag_data_deredden, mag_err, zeropointflux)
    
def pei92(wave_in_mu_m, Rv, Ebv, model='SMC') -> float:
    """
    ported from XSPEC originally by
    Martin.Still@gsfc.nasa.gov

    """

    if model=='MW':
        a=np.array([165.0, 14.0, 0.045, 0.002, 0.002, 0.012])
        lamb=np.array([0.047, 0.08, 0.22, 9.7, 18.0, 25.0])
        b=np.array([90.0, 4.0, -1.95, -1.95, -1.8, 0.0])
        n=np.array([2.0, 6.5, 2.0, 2.0, 2.0, 2.0])

    if model=='LMC':
        a=np.array([175.0, 19.0, 0.023, 0.005, 0.062, 0.02])
        lamb=np.array([0.046, 0.08, 0.22, 9.7, 18.0, 25.0])
        b=np.array([90.0, 5.5, -1.95, -1.95, -1.8, 0.0])
        n=np.array([2.0, 4.5, 2.0, 2.0, 2.0, 2.0])


    if model=='SMC':
        a=np.array([185.0, 27.0, 0.005, 0.01, 0.012, 0.03])
        lamb=np.array([0.042, 0.08, 0.22, 9.7, 18.0, 25.0])
        b=np.array([90.0, 5.5, -1.95, -1.95, -1.8, 0.0])
        n=np.array([2.0, 4.0, 2.0, 2.0, 2.0, 2.0])

    a_b = Ebv * (1.0 + Rv)

    # compute terms of sum

    ratio = wave_in_mu_m / lamb

    term = np.power(ratio, n)

    inv_term = 1.0 / term

    bottom = term + inv_term + b

    xi = np.sum(a / bottom)

    return a_b * xi
    
