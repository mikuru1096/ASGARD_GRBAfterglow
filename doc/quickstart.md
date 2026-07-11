# 快速开始

先按 [安装指南](installation.md) 构建 `electron_forward_fullhide_1d`。下面的脚本使用当前公开签名创建 top-hat/ISM 模型。

## 建立模型

```python
import numpy as np
from asgard_core import (
    Hadronic, Model, Numerics, Observer, ObserverGrid, Radiation,
    ReverseShock, SolverOptions, UniformMedium, top_hat_jet,
)

rad = Radiation(
    epsilon_e=0.1, epsilon_B=1e-3, p=2.3,
    proton_energy_fraction=0.0,
    magnetic_decay_alpha_t=0.0, magnetic_decay_t0_s=1.0,
    accelerated_electron_fraction=0.1, thermal_electrons=False,
    include_ssc=True, proton_synch=False, include_pgamma=False,
    bethe_heitler=False, hadronic_inverse_compton=False,
    pp=False, neutrino=False, acceleration_efficiency=1.0,
    reverse_proton_energy_fraction=0.0, pgamma_scheme="disabled",
    pair_production=False,
)

model = Model(
    jet=top_hat_jet(
        energy_iso_erg=1e52, initial_lorentz_factor=300.0,
        opening_angle_rad=0.1, shell_duration_s=None,
        magnetar=None, spreading=False,
    ),
    medium=UniformMedium(number_density_cm3=1.0),
    observer=Observer(
        z=0.1, viewing_angle_rad=0.0, viewing_azimuth_rad=0.0,
        luminosity_distance_cm=None,
    ),
    fwd_rad=rad,
    numerics=Numerics(
        num_radius=120,
        structured_num_theta=1, structured_num_phi=1,
        eats_num_theta=120, eats_num_phi=1,
        downstream_num_chi=None,
        num_observer_time=120, num_electron_gamma=121,
        num_photon_frequency=121, num_threads=8,
        electron_adaptive_substeps=False, electron_substep_rtol=0.02,
        electron_substep_min=100, electron_substep_max=1000,
        initial_radius_cm=1e14,
    ),
    observer_grid=ObserverGrid(time_min_s=1e2, time_max_s=1e7),
    solver_options=SolverOptions(
        electron_solver="fullhide_1d",
        dynamics_solver="forward_legacy",
        geometry_projection="sed_legacy",
        electron_photon_coupling="separated",
        ssc_cooling_mode="nakar_y_thomson",
        synchrotron_integration="fixed_grid",
        cooling_kernel="legacy",
        structured_backend="fortran_1d",
        patch_sampling="uniform",
        patch_sampling_pilot_theta=0, patch_sampling_num_times=12,
        patch_sampling_beaming_factor=3.0,
        patch_sampling_beaming_resolution=8.0,
        structured_parallel_mode="outer",
        structured_outer_threads=None, structured_inner_threads=None,
        fullhide2d_transport_model="legacy",
        fullhide2d_stochastic_accel_norm=0.0,
        fullhide2d_escape_mode="closed",
    ),
    reverse_shock=ReverseShock(
        enabled=False, shell_duration_s=10.0, upstream_sigma=0.0,
        include_cross_zone_ic=False, include_ssc=False,
    ),
    hadronic=Hadronic(
        enabled=False, solver="legacy_1d", num_proton_gamma=161,
        num_neutrino_frequency=121, pgamma_scheme="disabled",
        pair_cascade_iterations=1,
    ),
)
```

`structured_num_*` 控制 jet patch；`eats_num_*` 控制 EATS 角积分；`downstream_num_chi` 只服务有限 q-shell 路径。

## 查询观测量

```python
times = np.logspace(2, 7, 80)
nu = np.array([1e9, 1e14, 1e18])

grid = model.flux_density_grid(times, nu)
print(grid.total.shape)       # (frequency, time)
print(grid.fwd.sync.shape)

points = model.flux_density(times[:3], nu)
sed = model.spectrum(1e4, np.logspace(8, 25, 200))
band = model.flux(time_s=1e4, nu_min_hz=1e14, nu_max_hz=1e18)
track = model.details()
```

`flux_density` 按元素解释配对的 `(time, frequency)`；`flux_density_grid` 计算笛卡尔积。固定时刻宽频谱用 `spectrum`，频段积分用 `flux`。

## 单位与检查

公开 API 使用 cgs：时间 s、频率 Hz、距离 cm、角度 rad、能量 erg。未显式给 luminosity distance 时按红移计算。

首次结果至少检查：

- `total`、分量和内部粒子分布有限且非负；
- 光变和 `details()` 中的动力学量随时间连续；
- 改变网格后主要峰时、峰值和谱形收敛；
- 启用 RS 时，`upstream_sigma -> 0` 回到非磁化极限。

出现无物理来源的跳变时，应检查动力学、输运、源项和投影，不做 smoothing。更多 selector 与输出字段见 [公开 API](public_api.md)，可直接运行的查询片段见 [示例](examples.md)。
