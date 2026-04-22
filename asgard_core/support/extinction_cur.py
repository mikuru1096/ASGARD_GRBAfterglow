import numpy as np

def extinction_curve(um, av, ext_type='SMC'):
    pars = {
        'MW': (14.3, 6.49, 2.02, 0.0514),
        'SMC': (39.4, 3.89, 6.31, 0),
        'LMC': ...,
    }
    
    outArr = np.zeros_like(um)
    umSel = um >= 0.1216
    goodUm = um[umSel]
    outArr[~umSel] = np.nan

    def _general_func(par_tuple):
        c1, c2, c3, c4 = par_tuple
        p1 = c1 / ((goodUm / 0.08)**c2 + (0.08 / goodUm)**c2 + c3)
        p2 = 233 * (1 - c1 / (6.88**c2 + 0.145**c2 + c3) - c4/4.60) \
            / ((goodUm / 0.046)**2 + (0.046 / goodUm)**2 + 90)
        p3 = c4 / ((goodUm / 0.2175)**2 + (0.2175 / goodUm)**2 - 1.95)

        return av * (p1 + p2 + p3)

    outArr[umSel] = _general_func(pars[ext_type])

    return outArr

def get_abs(um, av, z, ar):

    umRest = um/(1.0+z)
    hostAbs = np.ones_like(um)
    umRestSel = umRest > 0.1216
    hostAbs[umRestSel] = extinction_curve(umRest[umRestSel], av, 'SMC')

    # 这里是手动指定 r 波段的吸收
    rBandSel = (um > 0.6) & (um < 0.68)
    hostAbs[rBandSel] = ar

    return hostAbs
