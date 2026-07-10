# 未修复缺陷

## γγ/电磁级联的壳层几何与传输所有权未闭合

### 第一性原理

当前 shell-level `annihilation` 与 pair cascade 固定采用共动路径
`R/(12 Gamma)`。该式不是一般 Blandford--McKee 壳宽，只能由 `k=0`
homogeneous slab 的质量守恒推出。若瞬时激波后 proper density 取
`n2'=4 Gamma n1(R)`，则

```text
Msw = 4 pi R^2 Delta' mp n2',
Delta' = Msw/[16 pi R^2 Gamma n1(R) mp].
```

对从原点延伸到当前半径的幂律外介质 `n1=A R^-k`，

```text
Msw = 4 pi mp n1(R) R^3/(3-k),
Delta' = R/[4(3-k) Gamma].
```

所以仅 `k=0` 给出 `Delta'=R/(12 Gamma)`；纯风 `k=2` 给出
`Delta'=R/(4 Gamma)`。对 density jump、tabulated profile 和
transrelativistic 流，`Msw+n1(R)` 至多定义把全部旧物质瞬时压到当前
激波后密度的 one-zone 等效宽度，不能替代演化后的多区体积与空间柱密度。

光深与 transfer 也必须分层。局域 γγ opacity 由目标光子柱密度
`Nnu=nnu Delta'` 决定；共空间、均匀发射与吸收使用
`psi(tau)=(1-exp(-tau))/tau`，前景屏使用 `exp(-tau)`，多区问题需要按
空间顺序积分。不能用一个固定路径或额外常数同时代表这三种几何。

### 受影响的真实调用

- `src/Radiation/pair_absorption.f90::annihilation`：从 `R,Gamma` 猜测
  `R/(12 Gamma)`，并在同一过程内同时组装 opacity 与有限单元 transfer。
- `asgard_core/asgard_state.py`：forward、reverse 和 hadronic target 可被合并后
  共用一个 absorption；pair cascade 开启时局域 target 置零，实际结果完全由
  `tau_extra` 中的 cascade 几何控制。
- `src/Structured/structured_jet_1d.f90::apply_absorption`：已有 `R_mass` 与局域
  density 状态，但仍只向辐射核传 `R,Gamma`；reverse/hadronic 输出复用同一个
  transfer。
- `prompt/radiation.py::compute_branch_radiation`：internal-shock branch 已有
  `comoving_volume_cm3`，仍错误套用外激波的 `R/(12 Gamma)` 路径。
- `src/Hadronic/hadronic_cascade.f90::cascade_seq`：当前存在三套没有共同所有者的
  时钟。pair electron 状态以 `Delta tobs/nsub` 推进，`electron_loss` 使用
  `R/(Gamma c)`，而 `tesc=R/(12 Gamma c)` 控制 `tau_pair`、光子源归一化和
  光子逃逸。

`src/Electron/electron_forward_transport_2d.f90` 通过显式 `dx_comov_hist`
调用 `pair_tau`，不属于这个固定路径缺陷；它说明路径应由几何所有者提供。

### 修复时序

该缺陷必须在强子阶段最后做一次原子修复。若只修改 shell annihilation，
pair-cascade 开启的公开路径仍由旧 `tau_extra` 控制；若只修改 cascade，普通
leptonic、structured 和 prompt 路径仍不一致。禁止在中间阶段加入 `k` switch、
经验系数、clamp、平滑或用 `R,Gamma` 猜测历史相关壳宽。

### 验收条件

- 建立纯函数式 opacity 核：输入目标光子柱密度或显式 density/path，输出
  `tau_gamma_gamma`；核内不读取 density profile，不决定 shell 几何，也不应用
  transfer。
- 几何所有者分别构造 FS、RS、prompt、structured 与 cascade 的光子柱密度；
  支持的 jump/多区路径必须使用演化后的分区体积或 `dx_comov`，没有物理状态的
  组合在公开边界明确拒绝，不能退回 `R/(12 Gamma)`。
- transfer 所有者明确区分共空间 `psi(tau)`、前景 `exp(-tau)` 与有序多区积分；
  FS+RS target 不再无条件合并后套一个共同路径。
- 保留现有 `Radiation.annihilation` 的 public f2py 参数顺序与数组 shape；通过薄包装
  明确保留其 `k=0 homogeneous slab` 语义，真实生产调用使用新的显式几何核心。
- `hadronic_cascade` 的 pair opacity、光子源与逃逸、pair electron 演化和
  `electron_loss` 时钟统一由显式几何/时钟所有者定义，不再各自使用互不相干的
  `Delta tobs`、`R/(Gamma c)` 或固定 `R/(12 Gamma c)`。
- 通过真实构建入口 `pair_absorption`、`structured_jet_1d`、
  `hadronic_forward_1d` 及受影响 source closure 的严格
  `-Wline-truncation -Werror=line-truncation` 检查。
- 直接运行并比较至少 `k=0`、wind、prompt、FS+RS、density jump、pair-cascade
  开/关路径；验证 `k=0` homogeneous 极限、wind 解析比例、prompt 体积路径和
  cascade 光深闭合。光深、频谱与时序曲线不得出现人为跳变，性能比较使用至少
  3 次 median。

## `total_only` 在分量求和前做对数频率插值，拟合总通量与公开分量总和不相等

### 第一性原理

对正的相邻谱点，投影核在对数频率上使用幂律插值。若两个物理分量分别为
`F1(nu)` 与 `F2(nu)`，一般有

```text
Ilog[F1 + F2] != Ilog[F1] + Ilog[F2].
```

每个辐射分量在一个频率单元内可分别近似为幂律；逐分量对数插值后求和能精确
保持这两个局部幂律，而先求和再把总谱强制近似成单一幂律会在分量交叉区产生
系统误差。`paper/main.tex` 的 observer contract 也明确规定返回总通量是 active
projected components 的代数和，而不是 nodewise total spectrum 的单次投影。

### 受影响的真实调用

- `asgard_core/asgard_state.py::observe_setup`：`mode="total_only"` 直接投影
  `components.total`；`mode="full_components"` 则逐分量投影后调用 `_sumobs`。
- `asgard_core/asgard_state.py::_observechi`：同样先合并 non-chi shell 分量，再做
  一次 shell-level 对数频率插值。
- `asgard_core/api_fit.py::_eval_cfg` 与
  `asgard_core/api_model.py::Model._total_matrix`：生产拟合和总通量快速路径均显式
  请求 `mode="total_only"`，因此该差异会进入 likelihood 与用户返回值。

当前真实 `standard` ISM、`fullhide_1d`、`num_nu=201`、正式 `band_freqs()` 入口的
直接比较得到：最大相对差 `5.617584054906704e-3`，active 点的 95% 分位差
`2.3393221690766247e-3`；最大差点的逐分量总和为
`8.86950613238617e-33`，`total_only` 为 `8.919612807067092e-33`。同一状态下
projection median 分别为 `5.013004e-3 s` 与 `2.321128e-3 s`，说明直接删除快速
路径虽能修正结果，却会把投影时间提高约 2.16 倍，不能作为最终算法。

### 验收条件

- `total_only` 与 `full_components["total"]` 必须采用同一个“逐分量频率插值后
  求和”离散定义；按相同线程数比较时只允许浮点归约顺序误差。
- 非强子 active 分量应由一个批量 Fortran EATS 核共享角向几何、arrival-time
  括号、Doppler/redshift 与目标频率定位，同时保留每个分量独立的正值对数插值；
  禁止继续投影预先求和的节点谱。
- 保留现有 `sed_interpolation` public f2py ABI、参数顺序与输出 shape；生产 Python
  路径可调用新的批量核心，单分量入口仍由同一核心退化得到。
- 受影响构建至少包括 `SED_interpolation`、`SED_interpolation_structured` 与引用
  投影源文件的 `structured_jet_1d`，并通过干净 source closure 的
  `-Wline-truncation -Werror=line-truncation` 检查。
- 直接运行 `Model.flux_density_grid`、拟合 `_eval_cfg`、prompt EATS、structured
  leptonic 入口；比较升序/乱序频率和观测时刻、单分量、FS synchrotron+SSC、
  FS+RS/cross-IC。总通量误差达到浮点舍入量级，频谱与光变保持非负、连续、平滑。
- 性能用相同已求解状态至少 3 次交错 median 比较；批量路径不得退化为逐分量重复
  调用，目标是在修正精度的同时把默认两分量投影保持在当前 `total_only` 同一量级。

## 自适应 theta midpoint 会永久丢失窄 EATS 支撑

### 第一性原理

`sed_adaptive_theta::integrate_theta_cell` 用父单元中点形成 coarse rule，再用两个子单元
中点形成 fine rule。这两个 midpoint rule 不嵌套：父中点不属于任一子规则。当前接受条件
还把 `flux_norm == 0` 直接视为收敛。因此，只要观测时刻的非零 EATS 角向支撑覆盖父中点、
却没有覆盖两个四分点，就会出现

```text
coarse_obs > 0, fine_obs == 0, err_norm > 0,
```

但代码仍接受 `fine_obs == 0`。即使只删除 `flux_norm == 0` 条件并继续细分，后续 dyadic
midpoint 也不再访问父中点；足够窄且位于父中点附近的支撑仍可在所有允许深度内被完全
漏掉。因此该缺陷不能用一个条件补丁修复，必须让嵌套求积保留已有非零样本，或先解析
确定当前观测时刻与径向 arrival-time 段在 theta 上的支撑交区，再在交区内积分。

### 受影响的真实调用与直接复现

- `asgard_core/asgard_postprocess.py::observe_flux` 在
  `geometry_kernel="sed_adaptive_theta"` 时调用该核。
- `prompt/eats.py::project_branch_flux` 在 `adaptive_max_depth > 0` 时调用同一核。
- `src/Interpolation/SED_interpolation.f90::sed_adaptive_theta` 的每个 theta/phi 基础单元均
  受影响；离轴窄 beaming、短径向发射窗、谱边界或 EATS 刚进入/离开一个角单元时风险最高。

直接 Fortran public 入口构造了一个物理正半径的窄 arrival window：单个 theta 单元
`[0,0.1]`，`R≈1e16 cm`，相邻径向到达时宽 `10 s`，查询时刻位于
`theta=0.05` 的父中点支撑内。结果为

```text
fixed NumTheta=1       10.201872969693007
fixed NumTheta=1000     0.061236149447951777
adaptive depth=1        0
adaptive depth=2        0
adaptive depth=4        0
adaptive depth=8        0
```

高分辨率固定积分证明真实角向积分非零；粗网格值虽不收敛，但它提供了必须保留的非零
支撑证据。自适应结果在所有深度精确为零，证明问题来自节点丢失而不是容差。

### 验收条件

- 采用严格嵌套的 theta 求积，或先解析裁剪 arrival-time 支撑后求积；禁止只修改
  `flux_norm == 0`、加入 flux floor、clamp、smoothing 或事后补点。
- 若采用 adaptive Simpson，应在 `mu=cos(theta)` 测度上定义被积函数，父单元端点与中点
  在两个子单元中复用；fine rule 必须包含 coarse rule 的全部采样信息，误差判据在
  `fine == 0, coarse != 0` 时不得接受。
- 若采用支撑交区法，应从每个径向 arrival-time 段的半开参数区间严格求出允许的 `mu`
  区间，保持现有 `0 <= Ratio < 1` 端点所有权；非单调径向段的多重交区不得合并或漏计。
- 上述窄窗复现必须收敛到高分辨率固定 theta 积分的非零结果；再比较 on-axis/off-axis、
  宽/窄径向窗、谱边界、乱序/重复观测时刻以及 prompt branch。
- 输出必须有限、非负，并随观测时刻和 viewing angle 连续平滑；对 theta 分辨率与最大深度
  做至少三层收敛比较，不允许用平滑后处理隐藏单元进入/离开造成的断崖。
- 构建 `SED_interpolation`，运行干净 source closure 的
  `-Wline-truncation -Werror=line-truncation`；性能以相同误差目标至少 3 次交错 median 比较，
  记录每个基础单元的实际谱投影评价数，而不是只比较最大深度。
