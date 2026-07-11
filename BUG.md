# 未修复缺陷

这里只记录已经由代码路径或数值结果确认、尚未修复的问题。修复并验收后删除对应条目。

## 5. Wind fullhide 2D 未回到 1D 极限

在相同 wind 物理参数下，`fullhide_2d` 的薄层/1D 极限相对 `fullhide_1d` 的有效时频单元误差为 median 0.934 dex、p95 2.542 dex。同一 runtime 配置中，1D electron solve 约 0.99 s，2D electron solve 约 0.08 s；这不是可接受的 2D 加速，而是两条路径没有解同一个极限问题。

修复前不得用 wind 1D/2D 的时间比宣称 scaling 收益。验收必须在相同半径、电子网格、冷却模式和投影几何下，先比较每个半径步的电子分布、同步辐射和 SSA，再比 observer flux。

## 6. Skymap 1D/2D 强度不闭合到同一物理极限

Skymap 像素积分与其自身 sample 总通量逐时刻精确闭合，视场也无漏光；但同一输入下 1D/2D 图像通量相差约 5--28 倍。该差异与 `fullhide_2d` 的 1D 极限失败同源，修复第 5 条之前不得把当前 1D/2D skymap 并列作为物理验证。

## 1. Detailed pp gamma 模型存在分段跳变

`pp_gamma_model` 的 Geant4、SIBYLL、QGSJET 和 Pythia8 路径忠实实现 AM3/Kafexhiu 分段参数化，但模型切换点本身不连续。已测得 1、4、20 GeV 附近约 7.44%、19.88%、15.51% 的功率跳变；100 GeV 处 Geant4 约 8.53%，SIBYLL/QGSJET 更大。

- 默认保持 `delta`；详细模型仅显式 opt-in。
- `qnu`、`qpair`、`ploss` 不受 gamma selector 影响。
- 禁止用 smoothing 或事后归一化消除跳变。
- 验收：采用有文献依据的连续模型，或从原始产生模型重新推导匹配条件，同时保持能量预算。

当前证据来自真实 `hadronic_forward_1d` 入口，而非孤立公式绘图。重新评估时必须同时比较微分谱、积分 gamma 功率、阈值行为和对应 pp loss；只让曲线视觉连续不能关闭本条。

当前选择是保留模型原式并清楚标注限制。

## 2. γγ/电磁级联的壳层几何与传输所有权未闭合

当前 formal hadronic 路径把 shell-level 光度、光子数密度、吸收深度和 observer flux 分散在多个层次。对于有限壳厚、RS/FS 跨区种子场和 pair cascade，路径长度、体积、共动/观测者坐标及转移次数没有统一所有者。

受影响入口包括 forward/reverse hadronic、cross-zone IC、pair injection 与 structured projection。继续扩展前必须：

1. 明确每个数组是 luminosity、emissivity、density 还是 observer flux。
2. 只在一个边界完成单位与坐标变换。
3. 由几何层唯一提供 path length、volume 与 SSA/γγ transfer。
4. 用无吸收极限、薄/厚极限和能量闭合验证真实入口。

不得通过重复乘壳厚、体积或 Doppler 因子补救输出。

建议修复顺序：先冻结单壳层坐标和单位契约，再闭合 FS 内部 transfer，然后加入 RS 跨区种子场，最后才允许 structured/χ 扩展。每一步都需要比较 absorption disabled 的旧路径，避免几何重构同时改变发射核。

关闭条件包括：

- 每个 transfer 数组在接口处标明 shape、单位、坐标系和所有者；
- SSA 与 γγ optical depth 都只由一次 path integral 产生；
- observer 层只做投影，不重新解释局部 density；
- pair cascade 每轮的吸收功率与次级注入功率闭合。

## 3. Joint secondary pair 坐标单位与状态/辐射重复所有权

joint cooling 中 secondary pair 既作为电子状态参与下一步冷却，又从辐射层的源项重新组装。部分路径以 `dN/dγ` 表示状态，部分中间量接近 `Q(γ)` 或单位体积率，生命周期和归一化边界不够明确，存在重复注入或漏乘时间步的风险。

修复要求：

- 输运层唯一拥有 pair state，辐射层只返回有明确单位的 source。
- 在一个位置完成 source × proper-time 到 state increment 的积分。
- separated 与 joint 对同一冻结背景的单步结果可逐数组比较。
- 总 pair energy 不超过 BH/pγ/γγ 注入预算，随时间连续。

还需覆盖 pair source 为零、只启用单个源项和联合源项三种真实配置。零源项必须自然给出零增量，不能依赖阈值 clamp；联合结果应等于相同冻结背景下各线性源项之和，直到反馈真正改变背景。

关闭本条前要记录电子状态在 radius step 前后分别代表什么时刻，并确认 observer radiation 读取的是更新后的唯一状态快照。

## 4. Reverse-shock BH electron grid 未覆盖完整运动学支撑

反向激波 Bethe–Heitler 次级电子目前复用有限 electron gamma grid；高能质子与 RS photon field 的组合可把次级电子推到网格上界之外，截断部分能量而没有完整的守恒诊断。

修复要求：由质子与光子运动学确定 secondary grid，上界收敛时 BH luminosity 和 pair energy 同时收敛；不得以末格堆积、外推或经验补偿恢复能量。

验收至少比较三组逐步扩大的上界，记录总 pair 数、总 pair energy 与辐射功率的收敛率；同时检查最低能格，避免扩大上界时隐式降低低能分辨率。FS 与 RS 使用同一运动学公式，但各自 photon field 和体积不得共享状态。
