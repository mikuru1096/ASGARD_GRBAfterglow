# ASGARD 待优化项

## 原则

- 不新建文件，不引入工具模块
- 不定义单行包装函数
- 不用 compose / pipe / 注册表 / Option 等抽象层
- 物理学内核代码完全不动（`*_kernel.f90`、`*_common.f90` 中的数值计算）
- 每次 Fortran 改动后做 compile check（`-Wline-truncation`）

## 问题 1: 2D 电子 solver 单行包装（-33 行）

**文件**: `src/Electron/FS_electron_charint_2d.f90`（19 行）、`FS_electron_fullhide_2d.f90`（24 行）

两个文件各是一个子程序，内容只有一行 `call fs_electron_transport_2d_core(..., .true., 'charint_2d')` 和 `call fs_electron_transport_2d_core(..., .false., 'fullhide_2d')`。唯一的差异是布尔标志和字符串标签。

**做法**: 删除包装，Python 直接调用 `fs_electron_transport_2d_core`，传入标志和标签参数。更新 `build_extensions.py` 的 `F2PY_ENTRYPOINTS` 和 `DIRECT_ORDERED_BUILD_MODULES`。

## 问题 2: `electron_transport_2d_core` 的 Boundary 解包重复

**文件**: `src/Electron/electron_transport_2d_kernel.f90` 和 `electron_seed_history_kernel.f90`

这两个 2D 专用 kernel 可能也有独立的 Boundary 解包（与 1D solver 相同逻辑）。需要验证后合并到 `electron_common.f90`。

## 问题 3: `_COOLING_KERNEL_ALIASES` — 死抽象（-8 行）

**文件**: `asgard_core/asgard_runtime.py:75-77, 504-508`

```python
_COOLING_KERNEL_ALIASES = {"legacy": "legacy"}

def _resolve_cooling_kernel(config: FitConfig) -> str:
    cooling_kernel = _COOLING_KERNEL_ALIASES.get(config.cooling_kernel.lower())
    if cooling_kernel is None:
        raise ValueError(...)
    return cooling_kernel  # 返回值从未被使用！
```

一个 dict 只有一条映射，函数返回值被丢弃。它唯一的作用是验证 `cooling_kernel == "legacy"`。

**做法**: 替换为 `if config.cooling_kernel.lower() != "legacy": raise ValueError(...)`，删除 dict 和函数。

## 问题 4: `_HADRONIC_SOLVER_ALIASES` — 微缩别名表（-6 行）

**文件**: `asgard_core/asgard_runtime.py:99-104, 511-515`

```python
_HADRONIC_SOLVER_ALIASES = {
    "legacy": "legacy_1d", "legacy_1d": "legacy_1d",
    "am3": "am3_1d", "am3_1d": "am3_1d",
}
```

4 行 dict + 5 行函数只为映射 2 个 canonical 值。`_resolve_hadronic_solver` 仅被调用 1 次。

**做法**: 内联到 `solve_hadronic` 中，三行 if/elif。

## 问题 5: `_PGAMMA_SCHEME_LEGACY_ALIASES` 两阶段查找（-3 行）

**文件**: `asgard_core/asgard_runtime.py:110-125, 518-525`

两阶段别名查找：先查 canonical（3 entries），miss 则查 legacy（8 entries）。可合并为一个 dict（11 entries），单次 `.get()`。

**做法**: 合并两个 dict，`_resolve_pgamma_scheme` 改为单次查找。

## 问题 6: `ELECTRON_1D_SOURCES` — 无用别名（-1 行）

**文件**: `build_extensions.py:86`

```python
ELECTRON_1D_SOURCES = ELECTRON_COMMON_SOURCES
```

该变量定义后从未被引用——1D solver 直接使用 `ELECTRON_COMMON_SOURCES`。

**做法**: 删除该行。

## 不做

- `DIRECT_ORDERED_BUILD_MODULES` 与 `F2PY_ENTRYPOINTS` 重复 — f2py 结构性约束，改它风险高于收益
- `asgard_runtime.py` `solve_electron` 7 路 dispatch — 各 solver 参数签名不同，强行统一会引入脆弱抽象
- `*_kernel.f90` 和 `*_common.f90` 中的所有数值计算子程序
- `implicit real(8)(A-H,O-Z)` 风格 — 丑陋但无害

## 预期效果

| 问题 | 文件 | 净删除 |
|------|------|--------|
| 2D solver 单行包装 | `FS_electron_*_2d.f90` | -33 |
| 2D Boundary 解包 | `electron_transport_2d_kernel.f90` 等 | ~-10 |
| `_COOLING_KERNEL_ALIASES` 死抽象 | `asgard_runtime.py` | -8 |
| `_HADRONIC_SOLVER_ALIASES` 微缩表 | `asgard_runtime.py` | -6 |
| pgamma 两阶段查找 | `asgard_runtime.py` | -3 |
| `ELECTRON_1D_SOURCES` 无用别名 | `build_extensions.py` | -1 |
| **合计** | | **~-60** |

## 验证

```bash
rtk python -m pytest tests/ -x --timeout=120 -q
rtk python tests/readme_smoke_bench.py
```

---

# OOP 消除方案

## 原则

- 类只应用于**必须封装可变状态**或**必须多态分发**的场景
- 能用 `(tuple, dict, namedtuple, function)` 就不用 `class`
- 消灭 type alias（`X = Y` 重命名）、无信息量 property、重复类型定义

## 问题 1: 5 个无用 type alias（-5 行，更新 ~30 处引用）

**文件**: `ASGARD/api_model.py:666-670`

```python
BranchView = FluxPair       # 仅用作 FluxPair 的别名
DetailView = CharTrack      # 仅用作 CharTrack 的别名
ModelDetails = TrackBundle  # 仅用作 TrackBundle 的别名
ModelFluxResult = FluxResult # 仅用作 FluxResult 的别名
ObsData = ObsSet            # 仅用作 ObsSet 的别名
```

每个别名被 5-10 个文件 import。它们不增加语义——`BranchView` 和 `FluxPair` 是同一个东西，有两个名字只会让人困惑该用哪个。

**做法**: 删除所有 5 个别名，所有 import 和调用处改用原名。受影响的文件：`api_observe.py`, `api_adaptive.py`, `api_fit.py`, `asgard_fit.py`, `ASGARD/__init__.py`, 若干 `tests/*.py`。

## 问题 2: `Scale` Enum 的重复成员（-5 行）

**文件**: `ASGARD/api_model.py:14-22`

```python
class Scale(str, Enum):
    LINEAR = "linear"
    linear = "linear"    # 重复！Enum 不允许同值成员
    LOG = "log"
    log = "log"
    LOG10 = "log10"
    log10 = "log10"
    FIXED = "fixed"
    fixed = "fixed"
```

10 个成员实际只有 5 个不同值。因为 `Scale(str, Enum)` 继承 str，`Scale.LINEAR == "linear"` 已经为 True，不需要小写版本。Enum 规范明确反对重复值。

**做法**: 删除 5 个小写成员。下游代码（`api_fit.py`, `asgard_fit.py`, tests）使用大写形式。

## 问题 3: `FluxPair.synch` — 冗余 property（-3 行）

**文件**: `ASGARD/api_model.py:445-447`

```python
@property
def synch(self) -> np.ndarray:
    return self.sync
```

字段名是 `sync`（缩写），property 是 `synch`（全拼）。二选一即可。

**做法**: 删除 property。如果偏好 `synch` 拼写，直接把字段名改为 `synch`。下游的 `result.fwd.synch` 调用改为 `result.fwd.sync`。

## 问题 4: 三个 `rvs` property — 无信息量别名（-6 行）

**文件**: `ASGARD/api_model.py`

`TrackBundle.rvs` → `return self.rev`（line 502）
`FluxResult.rvs` → `return self.rev`（line 514）
`CharTrack` 无关，`FluxPair` 无关。

两个 property 只是把 `rev` 改名叫 `rvs`。无信息增益。

**做法**: 删除两个 `rvs` property。下游 `.rvs` 调用改为 `.rev`。

## 问题 5: `ObsSet` — 可变列表包装（-50 行）

**文件**: `ASGARD/api_model.py:567-662`

`ObsSet` 是一个可变 dataclass，包含 3 个 list 和 3 个 `add_*` 方法。每个 `add_*` 方法内容相同：构造一个 dict append 到对应 list。

实际上整个使用模式是：创建空的 `ObsSet()`，逐个 `add_flux_density(...)` 追加数据，然后传给 `compile_problem`。这是典型的 builder 模式，用函数式写法更简洁：

```python
# 之前（OOP，可变）
data = ObsData()
data.add_flux_density(t=times, nu=freqs, f_nu=flux, err=errs)

# 之后（FP，不可变）
data = dict(flux_density=[
    dict(times_s=times, frequencies_hz=freqs, flux=flux, flux_err=errs)
])
```

**做法**: 删除 `ObsSet` 类及其 `add_*` 方法。`_compile_obs` 直接接受一个 dict。可选保留一个纯函数 `make_flux_density_entry(times_s, freqs, flux, flux_err)` 替代三个 `add_*` 方法。`api_fit.py` 中的 `Fitter.data` 改为存 dict 而非 `ObsData()`。

## 问题 6: 两套重复的 `FluxPair` / `FluxResult` / `CharTrack`（不做，记录）

**文件**: `asgard_core/asgard_types.py` vs `ASGARD/api_model.py`

两套独立的 dataclass 描述同样的物理量，但字段名不同（`time_s` vs `t_obs`, `gamma` vs `Gamma`）。内部模块用 `asgard_types.py` 版本，公共 API 用 `api_model.py` 版本，中间有转换。

**做法**: 这次不做——两套类型的字段语义不同（内部用 SI，API 用用户习惯单位/命名），转换有物理意义。但长期应该统一到一套类型。

## 预期效果

| 问题 | 删除行 | 新增行 | 净变化 |
|------|--------|--------|--------|
| 5 个 type alias | -5 | 0 | -5 |
| `Scale` 重复成员 | -5 | 0 | -5 |
| `FluxPair.synch` property | -3 | 0 | -3 |
| `rvs` property ×2 | -6 | 0 | -6 |
| `ObsSet` → 纯 dict | -50 | +10 | -40 |
| **合计** | **-69** | **+10** | **~-60** |

## 验证

所有改动是机械替换（重命名、删除 property、去 class 包装），不改变运行时行为：

```bash
rtk python -c "from ASGARD import Model, ISM, Wind, TophatJet, Observer; print('imports ok')"
rtk python -m pytest tests/ -x --timeout=120 -q
```
