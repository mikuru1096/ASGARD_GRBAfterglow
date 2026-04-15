# ASGARD 下一阶段主线计划

## 1. 当前基线

- `slc1_mmg2` 的 SSC 主链已经切到 `aux-grid -> ssc_spec(...)`。
- `mmg2` 激波前沿第一轮保峰修正已经落地，`tests/mmg2_front_sharpness_check.py` 可运行。
- `FS_electron_slc1` 的构建链已经有可复现入口：
  - 优先走原来的 `numpy.f2py -c`
  - 若 Windows + Python 3.12 + meson 再次触发 Fortran module 排序问题，则自动退到“按依赖顺序手工编译对象文件 + 最小 `pyf` + `f2py` 链接”
- `tests/order_convergence_check.py` 已升级为正式阶数诊断，输出：
  - 电子谱链阶数
  - 辐射谱链阶数
  - 独立动力学阶数

## 2. 当前真实结论

- 当前主线还不能宣称“电子谱和辐射谱都达到二阶以上”。
- 第一轮阶数诊断表明：
  - 电子谱链整体未达标。
  - 辐射谱链部分项目已达标，但并不完整。
  - 降阶不只存在于 `slc1_mmg2`，`fullhide/slc1` 的 observer-side 与部分 SSC/SSA 项也低于二阶。

## 3. 优先任务

### 3.1 电子谱链降阶定位

优先检查以下问题：

- `slc1_mmg2-electron-support-low`
  - 这是当前最差项，先查低能支撑边界口径。
- `slc1_mmg2-electron-shape-aligned`
  - 若峰位对齐后仍降阶，说明问题在谱形重建本身。
- `fullhide/slc1` 的电子谱低阶项
  - 先排除诊断口径误差，再决定是否需要改数值层。

重点排查层次：

1. `dN_x -> dN/dgamma` 的恢复口径
2. `work-grid -> public-grid` 投影
3. 低能边界与注入/冷却切换处的离散误差

### 3.2 辐射谱链降阶定位

优先检查以下问题：

- `observer-side total`
- `bands_flux`
- `nu_a`
- `fullhide` 的 `SSC`

重点排查层次：

1. synchrotron 共动谱积分
2. SSC 谱积分与投影
3. observer-side 插值与 EATS 映射
4. `nu_a` 搜索与插值口径

## 4. 执行顺序

1. 先只修一个明确降阶热点，不同时改多个层次。
2. 每做完一轮数值修改，先跑最小相关检查：
   - `python build_extensions.py --module FS_electron_slc1 --force`
   - `python tests/order_convergence_check.py`
   - `python tests/mmg2_front_sharpness_check.py`
3. 若阶数改善但物理曲线不连续或不平滑，则视为错误实现，继续回查。
4. 只有当电子谱链和辐射谱链都满足 `> 2`，这项任务才算完成。

## 5. 当前下一步

直接从 `slc1_mmg2` 的低能支撑边界降阶开始查：

- 对照 `work-grid` 与 `public-grid`
- 对照 `dN_x` 与 `dN/dgamma`
- 查低能边界重建是否在当前口径下引入系统降阶
