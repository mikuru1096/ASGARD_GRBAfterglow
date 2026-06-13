# 拟合工作流

本文说明如何用 `ASGARD.Fitter` 从多波段数据构建 likelihood、运行采样并检查物理结果。拟合不会改变物理模型；它只把参数映射到同一个 `Model` 求解链。

## 1. 数据准备

光变数据应整理为同长度数组：

```python
times_s = np.array([...], dtype=float)
nu_hz = np.array([...], dtype=float)
flux = np.array([...], dtype=float)
flux_err = np.array([...], dtype=float)
```

所有量必须使用 ASGARD 公开单位：秒、Hz、cgs flux density。若观测给出星等、\(\mu{\rm Jy}\) 或能段 flux，应在进入 `Fitter` 前转换。

## 2. 建立最小模型

拟合应从最小物理模型开始：

```python
from ASGARD import Fitter, ISM, Model, Observer, Param, Radiation, Scale, Setups, TophatJet

model = Model(
    jet=TophatJet(E_iso=1.0e52, Gamma0=300.0, theta_j=0.1),
    medium=ISM(n_ism=1.0),
    observer=Observer(z=0.1, theta_obs=0.0),
    fwd_rad=Radiation(eps_e=0.1, eps_B=1.0e-3, p=2.3, xi_N=0.1, ssc=True),
    setups=Setups(
        electron_solver="fullhide_1d",
        num_r=120,
        num_tobs=120,
        num_gam_e=121,
        num_nu=121,
    ),
)
```

先拟合正向激波 baseline。只有 residual 或物理目标要求时，才加入 wind、结构化喷流、反向激波或强子路径。

## 3. 绑定参数

```python
fitter = Fitter(model=model)
fitter.add_flux_density(times_s, nu_hz, flux, flux_err)

fitter.params = [
    Param("logE", "jet.E_iso", 50.0, 54.0, Scale.LOG10),
    Param("logn", "medium.n_ism", -4.0, 2.0, Scale.LOG10),
    Param("logepsB", "fwd_rad.eps_B", -6.0, -1.0, Scale.LOG10),
    Param("p", "fwd_rad.p", 2.0, 2.8, Scale.LINEAR),
]
```

`Param` 的第二个参数是 `Model` 内部字段路径。`Scale.LOG10` 表示采样变量 \(x\) 会被转换为 \(10^x\) 后写入模型。

## 4. 计算 likelihood

```python
ll = fitter.loglike({"logE": 52.0, "logn": 0.0, "logepsB": -3.0, "p": 2.3})
print(ll)
```

当前 likelihood 是独立高斯误差的 \(\chi^2\) 形式：

\[
\ln\mathcal{L}
=-\frac{1}{2}
\sum_i
\left[
\frac{F_{\nu,{\rm model}}(t_i,\nu_i)-F_{\nu,{\rm obs},i}}
{\sigma_i}
\right]^2.
\]

教程脚本生成的合成数据检查图：

![ASGARD 合成数据拟合检查](assets/tutorials/synthetic_fit.png)

若误差来自外部观测文件，应只在数据读入边界做单位和缺失值检查；内部模型路径不加兜底。

## 5. 运行采样

小规模检查可用内置 `fit`：

```python
result = fitter.fit(
    sampler="emcee",
    total_steps=256,
    burn_frac=0.5,
    thin=2,
    nwalkers=24,
)
print(result.best_params)
print(result.best_loglike)
```

正式采样时应把 quick grid 和 formal grid 分开：quick grid 用于检查参数范围、退化方向和明显 bug；formal grid 用于最终 posterior 产品。不要把 cached query timing 当成 cold solve timing。

## 6. 读取 posterior

`FitResult` 的主要字段：

- `best_params`：最大 likelihood 参数。
- `best_loglike`：对应 likelihood。
- `samples`：采样后的二维数组。
- `log_prob`：每个样本的 log probability。
- `labels`：样本列名。
- `top_k_params`：当前记录的优选参数集合。

如果安装了 `corner`，可以画 posterior：

```python
import corner

fig = corner.corner(result.samples, labels=result.labels)
fig.savefig("posterior_corner.png")
```

图像是 artifact。只有能由已提交脚本复现并记录命令时，才纳入版本库。

## 7. 拟合后的物理检查

拟合结果必须通过物理验收：

- best-fit 光变和 SED 连续平滑。
- \(\Gamma(R)\)、\(B(R)\)、\(\nu_m\)、\(\nu_c\)、\(\nu_a\) 没有孤立尖峰。
- 参数 posterior 不应靠无物理意义的边界取胜。
- 若需要反向激波，`reverse_sigma -> 0` 极限应回到非磁化结果。
- 若启用强子路径，应明确列出启用的 \(p\gamma\)、BH、pp、\(\gamma\gamma\) 或 neutrino 过程。

拟合中出现非平滑结构时，优先检查模型和数值路径，而不是修改 likelihood 或后处理。
