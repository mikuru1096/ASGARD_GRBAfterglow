# 术语与记号

## 写法

- 首次出现写中文名、英文名与缩写，之后使用缩写。
- FS/RS 分别指 forward/reverse shock；不要用“前/后激波”混淆传播方向。
- observer time、lab time、comoving time 必须显式区分。
- 上标 prime 表示共动系；`obs` 表示观测者系。

## 激波与几何

| 术语 | 含义 |
| --- | --- |
| FS / RS | 正向/反向激波 |
| shocked region | 已激波区，不等同整个 ejecta shell |
| EATS | 等到达时间面 |
| patch | 结构化喷流角网格单元 |
| q-shell / χ | 有限下游壳层及其自相似坐标 |
| LOS | 视线方向 |

## 粒子与辐射

| 记号 | 含义 |
| --- | --- |
| `dN/dgamma` | 每洛伦兹因子的粒子数分布 |
| synchrotron / SSC | 同步辐射 / 同步自康普顿 |
| SSA | 同步自吸收 |
| γγ | 双光子对产生 |
| pγ / pp | 光强子 / 质子—质子相互作用 |
| BH | Bethe–Heitler 对产生 |
| KN | Klein–Nishina 修正 |

`luminosity`、`emissivity`、`number density` 和 `observer flux` 不是同义词；文档和代码注释必须写明数组单位与坐标系。

## 接口与数值

- `Model`：不可观测物理配置与查询入口的组合。
- `Numerics`：网格规模与步进控制，不承载物理 selector。
- `SolverOptions`：算法、投影与并行 selector。
- `FluxResult`：分量化观测 flux density。
- `CharTrack` / `TrackBundle`：内部特征线和粒子状态诊断。
- separated/joint：光子—粒子源项分离或联合反馈。
- structured backend：结构化喷流 patch 的执行后端，不是 jet profile。

变量名中的 `_cm`、`_s`、`_hz`、`_erg` 是公开 API 的单位契约。代码内部旧短名只用于 ABI 或已冻结状态，不应扩散到用户文档。
