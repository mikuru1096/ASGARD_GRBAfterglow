# 拟合工作流

本文从零说明如何把观测数据放进 ASGARD、定义参数、检查 likelihood、运行采样，并判断结果是否可信。拟合不会改变物理模型；拟合器只做一件事：把一组参数写入同一个 `Model`，计算模型通量，再与数据比较。

如果只想查每个 API 选项能选什么，先看 `doc/public_api.md`。如果已经会建模，只想比较 `emcee` 和 `pymultinest`，看 `doc/mcmc_fitting.md`。

## 1. 拟合在做什么

ASGARD 当前拟合对象是多波段通量密度、固定时刻谱和频段积分通量。最常见的是通量密度点：

\[
\{t_i,\nu_i,F_{\nu,i},\sigma_i\}.
\]

给定参数 \(\boldsymbol{\theta}\)，`Fitter` 计算

\[
F_{\nu,{\rm model}}(t_i,\nu_i|\boldsymbol{\theta})
\]

并使用独立高斯误差：

\[
\ln\mathcal{L}(\boldsymbol{\theta})
=
-\frac{1}{2}
\sum_i
\left[
\frac{
F_{\nu,{\rm model}}(t_i,\nu_i|\boldsymbol{\theta})
-
F_{\nu,{\rm obs},i}
}{
\sigma_i
}
\right]^2 .
\]

这里没有额外的误差重标定、系统误差膨胀或 outlier rejection。若数据误差、单位或筛选有问题，应该在数据读入边界解决，不在 likelihood 内加补丁。

## 2. 准备数据

通量密度数据需要四个同长度数组：

```python
times_s = np.array([...], dtype=float)
nu_hz = np.array([...], dtype=float)
flux = np.array([...], dtype=float)
flux_err = np.array([...], dtype=float)
```

单位固定为：

| 数组 | 单位 | 常见来源 | 进入 ASGARD 前要做的事 |
| --- | --- | --- | --- |
| `times_s` | s | MJD、天、秒 | 转成爆发后秒。 |
| `nu_hz` | Hz | radio/optical/X-ray band | 每个点给一个有效频率。 |
| `flux` | cgs flux density | mJy、uJy、mag | 转成 \({\rm erg\,s^{-1}\,cm^{-2}\,Hz^{-1}}\)。 |
| `flux_err` | 同 `flux` | 测光误差 | 不能为 0；上限数据不要直接当高斯点。 |

常用转换：

\[
1\,{\rm Jy}=10^{-23}\,
{\rm erg\,s^{-1}\,cm^{-2}\,Hz^{-1}},
\qquad
1\,{\rm mJy}=10^{-26}.
\]

如果观测给的是星等，应先选定零点系统并转换为 flux density；如果给的是能段 flux，应使用 `add_flux`，不要假装它是某个单频点。

## 3. 建立最小模型

初学者先用 top-hat jet、均匀介质、正向激波同步辐射和 SSC。这个模型最容易判断单位、数据绑定和参数边界是否正确。

```python
import numpy as np

from asgard_core import (
    Fitter,
    Hadronic,
    Model,
    Numerics,
    Observer,
    ObserverGrid,
    Param,
    Radiation,
    ReverseShock,
    Scale,
    SolverOptions,
    UniformMedium,
    top_hat_jet,
)

model = Model(
    jet=top_hat_jet(
        energy_iso_erg=1.0e52,
        initial_lorentz_factor=300.0,
        opening_angle_rad=0.1,
        shell_duration_s=None,
        magnetar=None,
        spreading=False,
    ),
    medium=UniformMedium(number_density_cm3=1.0),
    observer=Observer(
        z=0.1,
        viewing_angle_rad=0.0,
        viewing_azimuth_rad=0.0,
        luminosity_distance_cm=None,
    ),
    fwd_rad=Radiation(
        epsilon_e=0.1,
        epsilon_B=1.0e-3,
        p=2.3,
        proton_energy_fraction=0.0,
        epsilon_b_floor=None,
        magnetic_decay_alpha_t=0.0,
        magnetic_decay_t0_s=1.0,
        accelerated_electron_fraction=0.1,
        thermal_electrons=False,
        include_ssc=True,
        include_kn_correction=False,
        proton_synch=True,
        include_pgamma=False,
        bethe_heitler=False,
        hadronic_inverse_compton=False,
        pp=False,
        neutrino=False,
        acceleration_efficiency=1.0,
        reverse_proton_energy_fraction=0.0,
        pgamma_scheme="disabled",
        pair_production=False,
    ),
    numerics=Numerics(
        num_radius=120,
        num_theta=120,
        num_phi=1,
        num_observer_time=120,
        num_electron_gamma=121,
        num_photon_frequency=121,
        num_chi=None,
        num_threads=8,
        electron_adaptive_substeps=False,
        electron_substep_rtol=0.02,
        electron_substep_min=100,
        electron_substep_max=1000,
        initial_radius_cm=1.0e14,
    ),
    observer_grid=ObserverGrid(time_min_s=1.0e2, time_max_s=1.0e7),
    solver_options=SolverOptions(
        electron_solver="fullhide_1d",
        dynamics_solver="forward_legacy",
        geometry_projection="sed_legacy",
        electron_photon_coupling="separated",
        ssc_cooling_mode="nakar_y_thomson",
        synchrotron_integration="fixed_grid",
        cooling_kernel="legacy",
        radiation_kernel="legacy",
        structured_backend="fortran_1d",
        patch_sampling="uniform",
        patch_projection="auto",
        patch_sampling_pilot_theta=0,
        patch_sampling_num_times=12,
        patch_sampling_beaming_factor=3.0,
        patch_sampling_beaming_resolution=8.0,
        structured_parallel_mode="outer",
        structured_outer_threads=None,
        structured_inner_threads=None,
        fullhide2d_transport_model="legacy",
        fullhide2d_stochastic_accel_norm=0.0,
        fullhide2d_escape_mode="closed",
    ),
    reverse_shock=ReverseShock(
        enabled=False,
        shell_duration_s=10.0,
        upstream_sigma=0.0,
        include_cross_zone_ic=False,
        include_ssc=False,
    ),
    hadronic=Hadronic(
        enabled=False,
        solver="legacy_1d",
        num_proton_gamma=161,
        num_neutrino_frequency=121,
        pgamma_scheme="disabled",
        pair_cascade_iterations=1,
    ),
)
```

上面的网格（120 个半径点、120 个角向点）用于演示完整模型构造逻辑；实际第一次拟合时，先使用第 6 节脚本里的 quick grid（72 个半径点、1 个角向点）检查单位、边界和参数退化。quick grid 和 formal grid 的关系见 `doc/mcmc_fitting.md` 第 5 节。

先不要打开 wind、结构化喷流、反向激波、强子或 \(\chi\) 分辨投影。只有 baseline 不能解释残差，且残差对应明确物理问题时，再加复杂模块。

## 4. 建立 Fitter

```python
fitter = Fitter(model=model)
fitter.add_flux_density(times_s, nu_hz, flux, flux_err)
```

`add_flux_density` 适合逐点光变数据。还有两个入口：

| 入口 | 数据形式 | 适合什么 |
| --- | --- | --- |
| `add_flux_density(times_s, frequencies_hz, flux, flux_err)` | 每个点一个 \(t,\nu,F_\nu\) | 多波段光变、最常见拟合。 |
| `add_spectrum(time_s, frequencies_hz, flux, flux_err)` | 一个时刻的多频谱 | SED 拟合或检查。 |
| `add_flux(nu_min_hz, nu_max_hz, time_s, flux, flux_err, num_points=64)` | 频段积分通量 | X-ray/高能频段 flux。 |

同一个 `Fitter` 可以加入多个数据块。likelihood 会把所有块的 \(\chi^2\) 相加。

`Fitter` 也可以在构造时直接接收参数列表：

```python
fitter = Fitter(model=model, params=[
    Param("logE", "jet.energy_iso_erg", 50.0, 54.0, Scale.LOG10),
])
```

教学文档常把 `fitter.params = [...]` 单独写出来，是为了让“建模、加数据、选参数”三个动作分开看清楚；两种写法进入同一 `fit()` 路径。

## 5. 选择拟合参数

参数用 `Param` 绑定到 `Model` 的字段路径：

```python
fitter.params = [
    Param("logE", "jet.energy_iso_erg", 50.0, 54.0, Scale.LOG10),
    Param("logn", "medium.number_density_cm3", -4.0, 2.0, Scale.LOG10),
    Param("logepsB", "fwd_rad.epsilon_B", -6.0, -1.0, Scale.LOG10),
    Param("p", "fwd_rad.p", 2.01, 2.8, Scale.LINEAR),
]
```

每个 `Param` 有五个含义：

| 位置 | 例子 | 意思 |
| --- | --- | --- |
| `name` | `"logE"` | 采样器看到的参数名。 |
| `path` | `"jet.energy_iso_erg"` | 写入 `Model` 的字段路径。 |
| `lower` | `50.0` | 采样变量下界。 |
| `upper` | `54.0` | 采样变量上界。 |
| `scale` | `Scale.LOG10` | 写入模型前是否做 \(10^x\) 变换。 |

`Scale.LOG10` 表示采样器中的 \(x\) 会写成 \(10^x\)。上例 `logE=52` 会把 `jet.energy_iso_erg` 设为 \(10^{52}\) erg。`Scale.LINEAR` 直接写入数值。`Scale.FIXED` 或 `lower == upper` 表示固定参数。

推荐的第一轮参数：

| 参数 | 路径 | 建议尺度 | 为什么先拟合 |
| --- | --- | --- | --- |
| \(E_{\rm iso}\) | `jet.energy_iso_erg` | `LOG10` | 控制整体能量和峰值。 |
| \(n\) | `medium.number_density_cm3` | `LOG10` | 控制减速、峰时和冷却。 |
| \(\epsilon_B\) | `fwd_rad.epsilon_B` | `LOG10` | 控制磁场、\(\nu_c\)、射电/光学比例。 |
| \(p\) | `fwd_rad.p` | `LINEAR` | 控制谱斜率和衰减斜率。 |

不要第一轮同时释放十几个参数。参数太多时，采样器会把单位错误、模型缺项和退化方向混在一起。

## 6. 一段完整的 emcee 拟合脚本

下面是一段可以从头阅读的最小脚本。它先用一个已知模型生成合成数据，再用 `emcee` 把能量、密度和 \(\epsilon_B\) 拟合回来。真实数据拟合时，只需要把 `times_s`、`nu_hz`、`flux`、`flux_err` 换成观测数据。

运行前需要当前环境能导入 `emcee`。若项目环境没有安装，可以用：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run --with emcee --with matplotlib python your_fit_script.py'
```

```python
import numpy as np
import matplotlib.pyplot as plt

from asgard_core import (
    Fitter,
    Hadronic,
    Model,
    Numerics,
    Observer,
    ObserverGrid,
    Param,
    Radiation,
    ReverseShock,
    Scale,
    SolverOptions,
    UniformMedium,
    top_hat_jet,
)


def build_quick_model() -> Model:
    return Model(
        jet=top_hat_jet(
            energy_iso_erg=1.0e52,
            initial_lorentz_factor=300.0,
            opening_angle_rad=0.1,
            shell_duration_s=None,
            magnetar=None,
            spreading=False,
        ),
        medium=UniformMedium(number_density_cm3=1.0),
        observer=Observer(
            z=0.1,
            viewing_angle_rad=0.0,
            viewing_azimuth_rad=0.0,
            luminosity_distance_cm=None,
        ),
        fwd_rad=Radiation(
            epsilon_e=0.1,
            epsilon_B=1.0e-3,
            p=2.3,
            proton_energy_fraction=0.0,
            epsilon_b_floor=None,
            magnetic_decay_alpha_t=0.0,
            magnetic_decay_t0_s=1.0,
            accelerated_electron_fraction=0.1,
            thermal_electrons=False,
            include_ssc=True,
            include_kn_correction=False,
            proton_synch=True,
            include_pgamma=False,
            bethe_heitler=False,
            hadronic_inverse_compton=False,
            pp=False,
            neutrino=False,
            acceleration_efficiency=1.0,
            reverse_proton_energy_fraction=0.0,
            pgamma_scheme="disabled",
            pair_production=False,
        ),
        rvs_rad=None,  # 反向激波独立辐射参数；普通公开路径保持 None。
        numerics=Numerics(
            num_radius=72,
            num_theta=1,
            num_phi=1,
            num_observer_time=72,
            num_electron_gamma=81,
            num_photon_frequency=81,
            num_chi=None,
            num_threads=1,
            electron_adaptive_substeps=False,
            electron_substep_rtol=0.02,
            electron_substep_min=100,
            electron_substep_max=1000,
            initial_radius_cm=1.0e14,
        ),
        observer_grid=ObserverGrid(time_min_s=1.0e2, time_max_s=1.0e7),
        solver_options=SolverOptions(
            electron_solver="fullhide_1d",
            dynamics_solver="forward_legacy",
            geometry_projection="sed_legacy",
            electron_photon_coupling="separated",
            ssc_cooling_mode="nakar_y_thomson",
            synchrotron_integration="fixed_grid",
            cooling_kernel="legacy",
            radiation_kernel="legacy",
            structured_backend="fortran_1d",
            patch_sampling="uniform",
            patch_projection="auto",
            patch_sampling_pilot_theta=0,
            patch_sampling_num_times=12,
            patch_sampling_beaming_factor=3.0,
            patch_sampling_beaming_resolution=8.0,
            structured_parallel_mode="outer",
            structured_outer_threads=None,
            structured_inner_threads=None,
            fullhide2d_transport_model="legacy",
            fullhide2d_stochastic_accel_norm=0.0,
            fullhide2d_escape_mode="closed",
        ),
        reverse_shock=ReverseShock(
            enabled=False,
            shell_duration_s=10.0,
            upstream_sigma=0.0,
            include_cross_zone_ic=False,
            include_ssc=False,
        ),
        hadronic=Hadronic(
            enabled=False,
            solver="legacy_1d",
            num_proton_gamma=161,
            num_neutrino_frequency=121,
            pgamma_scheme="disabled",
            pair_cascade_iterations=1,
        ),
    )


rng = np.random.default_rng(1234)
truth_model = build_quick_model()

times_s = np.logspace(3.0, 6.0, 12)
nu_hz = np.full(times_s.shape, 1.0e14)
truth = truth_model.flux_density(times_s, nu_hz).total
flux_err = 0.08 * truth
flux = truth + rng.normal(0.0, flux_err)

# 上面只用一个频率演示。真实多波段数据只需让 nu_hz 成为每个数据点
# 对应的频率数组；add_flux_density 支持每个点不同的 (t_i, nu_i)。

fit_model = build_quick_model()
fitter = Fitter(model=fit_model)
fitter.add_flux_density(times_s, nu_hz, flux, flux_err)
fitter.params = [
    Param("logE", "jet.energy_iso_erg", 50.0, 54.0, Scale.LOG10),
    Param("logn", "medium.number_density_cm3", -3.0, 2.0, Scale.LOG10),
    Param("logepsB", "fwd_rad.epsilon_B", -6.0, -1.0, Scale.LOG10),
]

result = fitter.fit(
    sampler="emcee",
    total_steps=256,
    burn_frac=0.5,
    thin=2,
    nwalkers=16,
)

print(result.best_params)
print(result.best_loglike)
print(result.labels)

t_plot = np.logspace(np.log10(times_s.min()), np.log10(times_s.max()), 80)
best_curve = fitter.flux_density_grid(result.best_params, t_plot, np.array([1.0e14])).total[0]

plt.errorbar(times_s, flux, yerr=flux_err, fmt="o", label="synthetic data")
plt.plot(t_plot, best_curve, label="best fit")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("observer time [s]")
plt.ylabel(r"$F_\nu$ [erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$]")
plt.legend()
plt.tight_layout()
plt.savefig("asgard_emcee_fit.png", dpi=220)
```

如果安装了 `corner`，可追加：

```python
import corner

fig = corner.corner(result.samples, labels=result.labels)
fig.savefig("asgard_emcee_corner.png", dpi=220)
```

合成数据脚本的意义不是证明模型正确，而是验证拟合管线：数据单位、`Param` 路径、likelihood、采样器和输出读取都能闭合。真实数据拟合前必须先完成下一节的 likelihood 手动检查。

## 7. 先手动检查 likelihood

运行采样前，先检查一个参数点：

```python
trial = {"logE": 52.0, "logn": 0.0, "logepsB": -3.0}
ll = fitter.loglike(trial)
print(ll)
```

如果你的 `fitter.params` 中包含 `Param("p", "fwd_rad.p", ...)`，则 `trial` 也必须加入 `"p": 2.3`。`trial` 的键要和 active 参数名一致。

再画同一参数点的模型曲线：

```python
times_plot = np.logspace(np.log10(times_s.min()), np.log10(times_s.max()), 120)
bands = np.unique(nu_hz)
# 如果 nu_hz 包含多个频率，bands 会是多元素数组，
# pred 的形状会变成 (len(bands), len(times_plot))。
pred = fitter.flux_density_grid(trial, times_plot, bands)
```

这一步要回答四个问题：

| 问题 | 失败时优先查什么 |
| --- | --- |
| `loglike` 是否是有限数？ | 数据数组、误差是否为 0、Fortran 扩展是否已构建。 |
| 模型量级是否接近数据？ | flux 单位、Jy/mJy/uJy 转换、红移和距离。 |
| 光变是否连续？ | 动力学、电子输运、密度剖面和投影路径。 |
| 每个频段是否映射到正确频率？ | `nu_hz` 是否按 Hz 输入，X-ray 是否误当成频率点。 |

## 8. 用 emcee 做第一轮采样

`emcee` 是默认选择，适合先看 posterior 形状、参数退化和模型是否明显不匹配。

```python
result = fitter.fit(
    sampler="emcee",
    total_steps=512,
    burn_frac=0.5,
    thin=2,
    nwalkers=24,
)
```

这些参数的含义：

| 参数 | 作用 | 选择建议 |
| --- | --- | --- |
| `total_steps` | 每个 walker 走多少步。 | quick check 用 128 到 512；正式运行按收敛诊断增加。 |
| `burn_frac` | 丢弃前多少比例的 flattened chain。 | 初学者用 0.5；正式结果检查 trace 后再定。 |
| `thin` | 每隔多少个样本保留一次。 | 只影响输出样本密度，不解决未收敛。 |
| `nwalkers` | walker 数量。 | 至少 \(2N_{\rm dim}+2\)，常用 4 到 8 倍维度。 |

当前实现会在每个 active 参数边界的 20% 到 80% 区间内初始化 walker，并在采样时拒绝越界点。

## 9. 用 PyMultiNest 做 evidence 和多峰检查

PyMultiNest 适合需要 Bayesian evidence、强多峰 posterior 或复杂退化结构的情况：

```python
result = fitter.fit(
    sampler="pymultinest",
    outputfiles_basename="chains/asgard-",
    verbose=True,
)
```

当前 `Fitter` 会调用 `pymultinest.solve.solve`。使用前必须确认 Python 包 `pymultinest` 和底层 MultiNest 动态库都能在当前环境中加载。`pymultinest` 路径使用每个 `Param` 的上下界作为均匀先验：

\[
x_i = x_{i,\min}+u_i(x_{i,\max}-x_{i,\min}),
\qquad
u_i\in[0,1].
\]

注意事项：

- 当前 `sampler="pymultinest"` 要求所有 `Param` 都是 active 参数；固定参数请直接写在 `model` 里，不要放进 `params`。
- `total_steps`、`burn_frac`、`thin`、`nwalkers` 是 emcee 参数，对 PyMultiNest 不起作用。
- `outputfiles_basename` 会生成 MultiNest chain 文件；输出目录应可写，正式结果要记录路径。
- evidence 只能比较同一数据、同一 likelihood、不同模型假设；不能用来比较单位不同或数据筛选不同的运行。

## 10. 读取结果

`fit` 返回 `FitResult`：

```python
print(result.best_params)
print(result.best_loglike)
print(result.labels)
```

字段含义：

| 字段 | emcee | PyMultiNest | 含义 |
| --- | --- | --- | --- |
| `best_params` | 有 | 有 | 当前样本中最大 likelihood 参数。 |
| `best_loglike` | 有 | 有 | 最大 log likelihood。 |
| `samples` | 有 | 有 | 采样变量数组，不是物理量数组；log 参数仍是 log 值。 |
| `log_prob` | 有 | 有 | 每个样本的 log likelihood/log probability。 |
| `labels` | active 参数名 | 参数名 | 样本列名。 |
| `logz`, `logzerr` | 无 | 有 | Bayesian evidence 及误差。 |

如果安装了 `corner`：

```python
import corner

fig = corner.corner(result.samples, labels=result.labels)
fig.savefig("posterior_corner.png", dpi=220)
```

图像是输出 artifact。只有脚本、数据、命令、git 版本和构建状态都可追溯时，才应提交图像。

## 11. 拟合后的物理检查

不要只看 `best_loglike`。最小检查包括：

- best-fit 光变和 SED 随时间/频率连续平滑。
- \(\Gamma(R)\)、\(B'(R)\)、\(\nu_m(R)\)、\(\nu_c(R)\)、\(\nu_a(R)\) 无孤立尖峰。
- posterior 不贴着无物理意义的边界取胜。
- 改变一个参数时，峰时和峰值连续变化。
- 如果启用反向激波，`upstream_sigma -> 0` 回到非磁化基线。
- 如果启用强子，列清 \(p\gamma\)、BH、pp、\(\gamma\gamma\)、neutrino 哪些过程打开，并检查 energy budget。

拟合中出现非平滑结构时，优先查动力学、电子输运、辐射源项和投影；不要改 likelihood 或对输出光变做 smoothing。

## 12. 什么时候增加复杂模型

复杂模块必须回答明确问题：

| 模块 | 什么时候加 | 加完必须检查 |
| --- | --- | --- |
| `WindMedium` | 早晚期斜率支持 wind-like 介质。 | \(k=2\) 边界、密度 floor/cap 是否影响结果。 |
| 结构化喷流 | off-axis 峰时/峰宽无法由 top-hat 解释。 | 角向 patch 收敛、viewing angle 连续性。 |
| `ReverseShock` | 早期光学/射电 excess 有 crossing 物理动机。 | `rev.sync` 是否来自区域 3，\(\sigma\to0\) 是否回归。 |
| `Hadronic` | 高能谱、neutrino 或二级反馈是科学目标。 | proton/secondary/photon sink/source 能量预算。 |
| `chi_eats_2d` | 研究有限厚壳层投影或 \(\chi\) 收敛。 | 只解释 FS synch+SSA；强子和 SSC 仍是壳层级。 |

不要为了“看起来完整”穷举无决策价值的 ablation。每次增加模块前先写清它要验证的物理假设。
