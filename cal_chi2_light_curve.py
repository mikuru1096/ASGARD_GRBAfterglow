import numpy as np
from pathlib import Path
from scipy import interpolate
import warnings

from extinc import opt_extinction

def cal_chi2_light_curve(bands_fit, model_curves, model_serial, frequency, 
                        Rv, Ebv, zeropointflux, f_sys):

    data_dir = Path(__file__).parent / "data_light_curve"
    
    if not data_dir.exists():
        raise FileNotFoundError(f"数据目录不存在: {data_dir}")
    
    chi2_total = 0.0
    
    # 预创建插值器（避免在循环中重复创建）
    model_interpolators = []
    for i in range(len(bands_fit)):
        interp = interpolate.interp1d(model_serial, model_curves[i, :], 
                                    kind='linear', bounds_error=False, 
                                    fill_value=np.nan)
        model_interpolators.append(interp)
    
    
    for data_file in data_dir.glob("*"):
        if not data_file.is_file():
            continue
            
        name = data_file.stem
        if name not in bands_fit:
            continue
            
        try:
            table = np.loadtxt(data_file)
            if table.ndim == 1:
                table = table.reshape(1, -1)
            
            range_data, flux_data, flux_err = _parse_observation_data(table, name)
            
            band_idx = bands_fit.index(name)
            if name.startswith('opt'):
                flux_data, flux_err = opt_extinction(
                    flux_data, flux_err, frequency[band_idx], Rv, Ebv, 
                    zeropointflux[band_idx])
            
            range_data = _convert_time_units(range_data, name)
            
            _validate_model_range(range_data, model_serial)
            
            fit_flux = model_interpolators[band_idx](range_data)
            if np.any(np.isnan(fit_flux)):
                raise ValueError("部分数据点超出模型范围")
            
            sigma = _get_uncertainties(flux_data, flux_err, fit_flux, table.shape[1])
            variance = _calculate_variance(flux_data, sigma, f_sys)
            
            temp_chi2 = np.sum((fit_flux - flux_data)**2 / variance)
            chi2_total += temp_chi2
            
        except Exception as e:
            warnings.warn(f"处理文件 {data_file.name} 时出错: {str(e)}")
            continue
    
    return chi2_total

def _parse_observation_data(table, name):
    """解析观测数据格式"""
    n_cols = table.shape[1] if table.ndim > 1 else table.size
    
    if n_cols == 6:
        range_data, flux_data = table[:, 0], table[:, 3]
        flux_err_up, flux_err_down = table[:, 4], -table[:, 5]
        flux_err = np.where(flux_data > table[:, 3], flux_err_up, flux_err_down)
    elif n_cols == 4:
        range_data = 0.5*(table[:, 0]+table[:, 1]) if table[0,0]<table[0,1] else table[:, 0]
        flux_data, flux_err = table[:, 2], table[:, 3]
    elif n_cols == 3:
        range_data, flux_data, flux_err = table[:, 0], table[:, 1], table[:, 2]
    elif n_cols == 2:
        range_data, flux_data = table[:, 0], table[:, 1]
        flux_err = flux_data * 0.1
    else:
        raise ValueError(f'观测数据应为2-6列，当前为{n_cols}列')
    
    return range_data, flux_data, flux_err

def _convert_time_units(range_data, name):
    """转换时间单位"""
    return range_data if name.endswith('keV') else range_data * 86400

def _validate_model_range(range_data, model_serial):
    """验证模型范围是否覆盖数据"""
    if np.min(range_data) < model_serial[0] or np.max(range_data) > model_serial[-1]:
        raise ValueError('模型曲线不能完全覆盖数据范围')

def _get_uncertainties(flux_data, flux_err, fit_flux, n_cols):
    """获取不确定性估计"""
    return np.where(fit_flux > flux_data, flux_err, -flux_err) if n_cols == 6 else flux_err

def _calculate_variance(flux_data, sigma, f_sys):
    """计算方差"""
    return (flux_data * 0.1)**2 if f_sys <= 0 else (flux_data * f_sys)**2 + sigma**2
    
    
