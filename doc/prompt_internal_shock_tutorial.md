# Prompt 内部激波快照教程

`prompt/` 用于两壳碰撞的快照、FS/RS jump、辐射与 EATS 诊断。它不是稳定 `asgard_core` 公共 API，也不进入正式拟合链。

## 1. 当前入口

- `prompt.internal_shock`：壳层运动学、碰撞与 jump；
- `prompt.radiation`：电子、同步、SSC 和 gamma-gamma；
- `prompt.eats`：脉冲到达时投影。

先查看函数签名和模块内示例；这些入口允许随诊断需求演化。

## 2. 两壳初态

设慢壳和快壳的实验室系宽度、质量、洛伦兹因子与发射间隔为

\[
(\Delta_s,M_s,\Gamma_s),\qquad
(\Delta_f,M_f,\Gamma_f),\qquad \Delta t_{\rm ej},
\]

且 \(\Gamma_f>\Gamma_s\)。速度由

\[
\beta_i=\sqrt{1-\Gamma_i^{-2}}
\]

精确计算。对很大的 \(\Gamma\)，直接相减两个接近 1 的速度会损失精度；实现使用稳定的相对运动表达式。

若壳宽定义为实验室系宽度，未扩展近似下追赶半径由两壳前后沿轨迹相交确定。常见尺度

\[
R_{\rm col}\sim
\frac{2\Gamma_s^2c\Delta t_{\rm ej}}
{1-\Gamma_s^2/\Gamma_f^2}
\]

只用于估算，实际例程保留有限宽度和边界位置。

碰撞前必须满足壳层次序、正宽度和正质量。没有追赶事件时应暴露无激波状态，而不是生成弱伪激波。

## 3. 壳层密度与磁化

冷壳共动静质量密度由其质量与共动体积确定：

\[
\rho_i'=\frac{M_i}{\Omega R^2\Delta_i'},
\qquad \Delta_i'=\Gamma_i\Delta_i
\]

其中 \(\Omega\) 是壳占据立体角。若输入采用各向同性等效能量，质量与 \(\Omega\) 的换算必须只做一次。

上游磁化定义为

\[
\sigma_i=\frac{B_i'^2}{4\pi\rho_i'c^2}.
\]

若总能量包含物质与 Poynting 分量，物质质量不能再按总能量重复归一。代码中的能量拆分必须与该 \(\sigma\) 定义相配。

## 4. 相对洛伦兹因子

两区域相对洛伦兹因子采用

\[
\Gamma_{ab}=\Gamma_a\Gamma_b(1-\beta_a\beta_b).
\]

超相对论近似

\[
\Gamma_{ab}\simeq\frac12
\left(\frac{\Gamma_a}{\Gamma_b}+
\frac{\Gamma_b}{\Gamma_a}\right)
\]

可用于检查，但不能覆盖接近等速壳时的精确分支。

接触面洛伦兹因子 \(\Gamma_c\) 位于 \(\Gamma_s\) 与 \(\Gamma_f\) 之间。其两侧分别形成进入慢壳的 FS 和进入快壳的 RS。

## 5. 压力平衡

接触面两侧下游总压力满足

\[
p_{2,\rm tot}'(\Gamma_{cs},\rho_s',\sigma_s)
=p_{3,\rm tot}'(\Gamma_{fc},\rho_f',\sigma_f).
\]

总压力包含热压与由横向场压缩产生的磁压。求根变量为 \(\Gamma_c\)，物理解必须：

- 位于两壳洛伦兹因子之间；
- 两侧均存在压缩激波；
- 在 \(\sigma\to0\) 时连续回到流体 jump；
- 压力残差随求根收敛。

若一侧相对运动不足以形成激波，应返回无效 shock branch，而不是用任意下限制造耗散。

## 6. FS/RS jump

每侧 jump 由冷上游质量、相对洛伦兹因子和磁化决定。输出至少包含：

\[
\rho_d',\quad e_{\rm th,d}',\quad B_d',\quad
u_{u,s},\quad u_{d,s}.
\]

下标 \(s\) 表示激波静止系。质量通量连续给出

\[
\rho_u'u_{u,s}=\rho_d'u_{d,s}.
\]

磁化 jump 使用冷上游、横向场的一维理想 MHD 条件。\(\sigma=0\) 精确落到流体分支；有限 \(\sigma\) 改变可供非热粒子使用的热能。

微物理参数 \(\epsilon_e,\epsilon_B,\xi_e,p\) 应作用于 jump 后可用热能。上游有序场与微湍动 `epsilon_B` 场不能重复计入。

## 7. Crossing 与 branch history

FS 穿过慢壳、RS 穿过快壳的时间由激波系速度与有限壳宽决定。每条 branch 保存其活动区间、扫掠质量和下游状态。

在 crossing 前，新的上游物质持续注入电子；crossing 后，注入终止，已有粒子只经历辐射和绝热冷却。

branch history 必须保持：

- 时间严格递增；
- 扫掠质量单调且不超过壳质量；
- crossing 前后状态连续；
- 无效激波分支不产生辐射。

快照工作流没有把分支历史升级为长期 afterglow 动力学。

## 8. 电子与磁场

每个有效激波区按其热化能注入

\[
Q_e(\gamma)=Q_0\gamma^{-p},
\qquad \gamma_m\le\gamma\le\gamma_M.
\]

\(Q_0\) 同时满足电子数和能量预算。最小电子洛伦兹因子依赖每个加速电子可得能量；最大值由加速和损失时间竞争确定。

辐射磁场由压缩有序场与允许的微湍动分量共同定义。两者的能量来源必须区分，避免同一磁能重复用于 jump 和辐射。

## 9. 同步与 SSA

局域同步频率为

\[
\nu_c'=\frac{3eB'}{4\pi m_ec}\gamma^2\sin\alpha_p.
\]

同步发射率与 SSA 系数由同一电子谱积分。均匀区传输采用

\[
T_{\rm ssa}=\frac{1-e^{-\tau_{\rm ssa}}}{\tau_{\rm ssa}}.
\]

FS 和 RS 使用各自的体积、路径长度、电子谱和磁场，不共享一个经验自吸频率。

## 10. SSC 与 gamma-gamma

SSC 目标场来自同一区域的同步光子。若计算 cross-zone IC，必须显式变换另一侧光子场的能量、方向与几何稀释。

gamma-gamma 光深由局域目标光子谱与路径长度计算。吸收后的 pair 只有在相应二级链实现时才进入后续辐射；不能把被吸收能量自动等同于可见级联。

## 11. EATS 脉冲

对每个发射事件，观测时间为

\[
t_{\rm obs}=(1+z)\left(t-\frac{R\cos\alpha}{c}\right),
\]

Doppler 因子为

\[
\delta_D=[\Gamma(1-\beta\cos\alpha)]^{-1}.
\]

EATS 对壳面角度和 branch history 积分。壳宽、激波 crossing 与角延迟共同决定脉冲宽度；文档不把它们统一替换成 \(R/(12\Gamma)\) 或单一角时间尺度。

## 12. 最小诊断流程

```python
from prompt.internal_shock import simulate_internal_shock
from prompt.radiation import compute_prompt_observed_flux

shock = simulate_internal_shock(...)
flux = compute_prompt_observed_flux(shock, ...)
```

函数名与参数以当前源码签名为准；该示例表达调用顺序，不承诺 public API 稳定性。

## 13. 运行检查

```bash
uv run python -m prompt.run_snapshot
```

检查结果应包括：

1. \(R_{\rm col}\)、\(\Gamma_c\) 有限且处于物理解区间；
2. FS/RS 压力残差收敛；
3. 扫掠质量不超过各自壳质量；
4. crossing 前后谱和光变连续；
5. 磁化趋零连续回到流体结果；
6. 无效 shock branch 的耗散与辐射为零；
7. 频谱阈值、SSA turnover 和 gamma-gamma cutoff 符合物理次序；
8. EATS 通量有限、非负，时间轴严格递增。

## 14. 当前边界

- 这是两壳、局域一维 jump 的诊断模型。
- 不包含多壳长期碰撞树或全局 MHD 演化。
- 不自动与 afterglow 拟合参数共享能量预算。
- 不把快照 pair 产生解释为完整含时级联。
- 不通过 clamp、平滑或经验 pulse 模板修补不连续。

总物理约定见[物理模型](physics_model.md#20-prompt)。
