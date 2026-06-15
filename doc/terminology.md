# 术语表

本文固定 ASGARD 中文文档的术语。代码标识、API 名、文件名和常用缩写保持英文；解释性文本使用中文。

## 激波与几何

| 英文/缩写 | 中文写法 | 说明 |
| --- | --- | --- |
| forward shock, FS | 正向激波 | 不再使用“正激波”或“前向激波”。 |
| reverse shock, RS | 反向激波 | 不再使用“反激波”“后向激波”或“逆向激波”。 |
| ejecta reverse shock | 抛射物反向激波 | 指原始 ejecta 被反向激波加热形成区域 3 的分支；避免写成“原始反向激波”。 |
| secondary reverse shock | 次级反向激波 | 首次出现可写“密度增强触发的次级反向激波”；不写“secondary reverse shock”混排。 |
| equal-arrival-time surface, EATS | 等到达时间面 | 用于观测者投影。 |
| Doppler factor | Doppler 因子 | 物理量写 \(\delta=[\Gamma(1-\beta\mu)]^{-1}\)；插值代码变量 `doppler` 表示 \(D=\delta^{-1}\)。 |
| structured jet | 结构化喷流 | 使用 `top_hat_jet`、`gaussian_jet`、`power_law_jet` 等函数式构造器。 |
| shell | 壳层 | 指半径网格上的动力学/输运单元。 |
| shell-level | 壳层级 | 指每个半径壳层只有一个局域物理状态；不要写成 shell level 或 shell contract。 |
| patch | 角向面元 | 结构化喷流或偏振投影中的角向单元。 |
| \(\chi\)-resolved | \(\chi\) 分辨 | 正文/公式用 \(\chi\)，代码名保留 `chi` 和 `chi_eats_2d`。 |
| finite-thickness shell | 有限厚壳层 | 指 `chi_eats_2d` 的正向激波同步辐射+SSA 投影几何。 |
| density bump | 密度增强 | 首次可注明 smooth density bump；避免中英文混用。 |

## 辐射与粒子过程

| 英文/缩写 | 中文写法 | 说明 |
| --- | --- | --- |
| synchrotron | 同步辐射 | 代码中 `synch` 可保留。 |
| synchrotron self-absorption, SSA | 同步自吸收 | 缩写 SSA 可保留。 |
| synchrotron self-Compton, SSC | 同步自康普顿 | 缩写 SSC 可保留。 |
| inverse Compton, IC | 逆康普顿 | cross-zone IC 可写作跨区逆康普顿。 |
| gamma-gamma pair production | \(\gamma\gamma\) 对产生 | 中文标题可写“对产生”。 |
| pair cascade | 对级联 | 当前主链是 \(\gamma\gamma\) 对产生/同步辐射级联。 |
| hadronic | 强子 | hadronic path 写作强子路径。 |
| Bethe-Heitler, BH | Bethe-Heitler 过程 | 缩写 BH 可保留。 |
| photopion, p-gamma | 光介子过程，\(p\gamma\) | 公式中写 \(p\gamma\)。 |
| polarization | 偏振 | Stokes \(I,Q,U\) 使用 LaTeX。 |

## 数值和接口

| 英文/缩写 | 中文写法 | 说明 |
| --- | --- | --- |
| public API | 公开 API | `Model`、`Fitter`、`Param` 等保持英文。 |
| backend | 后端 | 指 Fortran/Python 内部实现路径。 |
| kernel | 数值核 | 高代价 Fortran 算子。 |
| solver | 求解器 | `fullhide_1d` 等名称保持英文。 |
| benchmark | 基准测试 | 图像/CSV 必须可由脚本复现。 |
| smoke test | 冒烟测试 | 用于快速验证接口和物理边界。 |
| fallback | 回退/兜底 | 项目约束中禁止无物理依据的 fallback；若是适用域外的明确物理分支，应写明边界。 |

## 公式书写

公式统一使用 LaTeX：

- 行内公式使用 `\(...\)`，例如 \(\Gamma_0=300\)。
- 独立公式使用 `\[...\]`，例如

\[
\frac{\mathrm{d}R}{\mathrm{d}t_{\rm lab}}=\beta c.
\]

不得用纯文本伪公式替代数学表达式。代码片段中的变量名不视为公式。
