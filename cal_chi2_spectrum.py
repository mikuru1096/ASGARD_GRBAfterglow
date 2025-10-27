import numpy as np
from pathlib import Path
from scipy import interpolate
import warnings
import gc

def cal_chi2_spectrum(name_fit, model_curves, model_serial):

    data_dir = Path(__file__).parent / "data_spectrum"
    
    if not data_dir.exists():
        raise FileNotFoundError(f"数据目录不存在: {data_dir}")
    
    chi2_total = 0.0
    
    # 预创建插值器
    model_interpolators = []
    for i in range(len(name_fit)):
        interp = interpolate.interp1d(model_serial, model_curves[i, :], 
                                    kind='linear', bounds_error=False, 
                                    fill_value=np.nan)
        model_interpolators.append(interp)

    for data_file in data_dir.glob("*"):
        if not data_file.is_file():
            continue
            
        name = data_file.stem
        if name not in name_fit:
            continue

        try:
            table = np.loadtxt(data_file)
            if table.ndim == 1:
                table = table.reshape(1, -1)
            
            range_data, flux_data, flux_err = _parse_spectrum_data(table, name)
            
            _validate_model_range(range_data, model_serial)
            
            band_idx = name_fit.index(name)
            fit_flux = model_interpolators[band_idx](range_data)

            if np.any(np.isnan(fit_flux)):
                raise ValueError("部分数据点超出模型范围")
            
            sigma = _get_spectrum_uncertainties(flux_data, flux_err, fit_flux, table.shape[1])
            
            num_data = len(range_data)
            temp_chi2 = np.sum(((fit_flux - flux_data) / sigma) ** 2) / num_data
            chi2_total += temp_chi2
            
        except Exception as e:
            warnings.warn(f"处理文件 {data_file.name} 时出错: {str(e)}")
            continue
    
    gc.collect()
    return chi2_total

def _parse_spectrum_data(table, name):
    """解析能谱观测数据格式"""
    n_cols = table.shape[1] if table.ndim > 1 else table.size
    
    if n_cols == 6:
        range_data, flux_data = table[:, 0], table[:, 3]
        flux_err_up, flux_err_down = table[:, 4], -table[:, 5]
        flux_err = np.where(flux_data > table[:, 3], flux_err_up, flux_err_down)
    elif n_cols == 4:
        range_data = table[:, 0]
        flux_data, flux_err = table[:, 2], table[:, 3]
    elif n_cols == 3:
        range_data, flux_data, flux_err = table[:, 0], table[:, 1], table[:, 2]
    elif n_cols == 2:
        range_data, flux_data = table[:, 0], table[:, 1]
        flux_err = flux_data * 0.1
    else:
        raise ValueError(f'观测数据应为2-6列，当前为{n_cols}列')
    
    return range_data, flux_data, flux_err

def _validate_model_range(range_data, model_serial):
    """验证模型范围是否覆盖数据"""
    if np.min(range_data) < model_serial[0] or np.max(range_data) > model_serial[-1]:
        raise ValueError('模型曲线不能完全覆盖数据范围')

def _get_spectrum_uncertainties(flux_data, flux_err, fit_flux, n_cols):
    """获取能谱不确定性估计"""
    if n_cols == 6:
        # 对于有上下误差的数据，根据拟合值与观测值的关系选择误差方向
        flux_err_up, flux_err_down = flux_err, -flux_err
        return np.where(fit_flux > flux_data, flux_err_up, flux_err_down)
    else:
        # 对于其他格式，直接使用给定的误差
        return flux_err

# 测试代码
if __name__ == '__main__':
    # 示例数据
    time = ['spectrum1']
    energy_curves = np.array([[0.01, 10]])
    energy_serial = np.array([0.01, 20])
    
    result = cal_chi2_spectrum(time, energy_curves, energy_serial)
    print(f"卡方值: {result}")
