import numpy as np
from scipy.interpolate import interp1d, interp2d

para_tev2hz = 2.418e26
def cal_ebl(z, v_obs, model="Dominguez11.txt"):
    file_path = Path(__file__).parent / "EBL" / model
    table = np.loadtxt(file_path)
    redshifts = table[0, 1:]
    energies = table[1:, 0] * para_tev2hz
    tau_values = table[1:, 1:]

    tau_interp = interp2d(redshifts, energies, tau_values, kind='linear')
    
    index = np.argmax(redshift > z)
    tau_z = tau_interp(z, energies).flatten()
    
    # 创建tau关于能量的插值函数
    tau_z_interp = interp1d(energies, tau_z, kind='linear', 
                           bounds_error=False, 
                           fill_value=(tau_z[0], tau_z[-1]))
    
    # 计算吸收系数
    absorption = np.exp(-tau_z_interp(v_obs))
    
    # 处理超出范围的值
    absorption[v_obs < energies[0]] = 1.0
    absorption[v_obs > energies[-1]] = 1.0e-30
    
    return absorption





