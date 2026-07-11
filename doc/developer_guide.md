# 开发与验证指南

本文补充根目录 `AGENTS.md`，后者始终优先。

## 1. 固定环境

唯一工作树是：

```text
/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow
```

默认使用 WSL Ubuntu、项目 `uv` 和 `/usr/bin/gfortran`。所有命令以 `rtk`
开头；非交互 WSL 先 `source ~/.wsl_env`。开始前记录 branch、HEAD 和工作树：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && git status --short --branch && git log -1 --oneline'
```

## 2. 变更原则

- 从物理或接口问题出发，不创建没有明确效用的路径；
- 保留用户已有改动，不使用 destructive Git 命令；
- Python 负责编排，Fortran 负责数值核；
- 只在系统边界校验，内部状态遵守已建立 contract；
- 不添加 fallback、smoothing、经验归一化或事后修补；
- 新 bug 先登记 `BUG.md`，修复并验收后同轮删除；
- 临时 driver、profile 和 benchmark 放 `/tmp`，不新增测试文件。

## 3. 修改前

1. 写清收益假设、物理假设和可观察验收量；
2. 从公开入口沿调用链定位到数值核；
3. 搜索 Python caller、f2py export 和 build source closure；
4. 记录旧实现输出、运行时间和峰值临时内存；
5. 确认未支持边界，避免无意扩大 public API。

删除代码必须同时证明没有运行 caller、构建节点和外部契约。删除文档必须修复导航
与链接，不保留“已删除页面”的占位文件。

## 4. Fortran 修改

优先精简循环、临时数组和薄 wrapper，不引入 manager/context。保持 ABI 参数顺序、
shape 和 dtype，除非任务明确修改公开接口。

声明块保持紧凑，同类型、同语义变量合并；新局部变量名最多一个下划线。数值表达式
的求值顺序若改变，必须用物理预算和精度比较验收。

受影响扩展强制构建：

```bash
TMPDIR=/tmp uv run python build_extensions.py --module MODULE_NAME --force
```

随后从干净 `/tmp` 按 `build_extensions.py` 的真实有序 source closure 执行
`-Wline-truncation`。不要在仓库根目录拾取旧 `.mod`，也不要用邻近文件集合代替
真实闭包。

## 5. Python 修改

- 字符串 selector 只在公开边界映射一次；
- 保留承载单位、坐标、线程 callback 或 f2py ABI 的函数边界；
- 只在首个运算已产生独立数组后使用原位累加；
- 不回写 public 输入数组；
- cache key 必须覆盖所有改变结果的参数；
- dataclass 字段有物理意义时不为缩行删除。

公开参数或返回值变化时同步 `public_api.md` 和 backend limits。

## 6. 数值验收

结果比较至少覆盖：

- shape、dtype、内存布局和输入不变性；
- callback 次数、参数和 report 顺序；
- 禁用新过程时回到 baseline；
- 低阶/薄厚/无吸收等解析极限；
- 粒子数、能量和 photon survival 预算；
- 时间、半径、能量及观测曲线的连续平滑性。

若优化只引入末位浮点差异，可以接受，但必须定位到等价求值顺序或更少临时量。
性能至少运行三次，使用 median；大型 observer/hadronic 汇总同时记录峰值内存。
无可测收益的复杂改动撤回。

## 7. 最小验证矩阵

| 变更 | 最小真实入口 |
|---|---|
| Dynamics | forward/reverse 强制构建与 reverse smoke |
| Electron 1D | 对应 solver 构建与公开 `Model` 查询 |
| Electron 2D | 2D 构建与 χ state/projection 查询 |
| Radiation | electron/radiation extension 与 SED/lightcurve |
| Hadronic | forward、必要时 reverse/structured formal 入口 |
| Observer | lightcurve、SED、total/components |
| 文档 | encoding、strict MkDocs、diff check |

跨阶段修改运行所有受影响行，不为“完整”执行没有决策价值的穷举。

## 8. 文档验证

```bash
uv run python tools/check_text_encoding.py
uv run --with "mkdocs<2" --with "mkdocs-material>=9.5" mkdocs build --strict --site-dir /tmp/asgard_mkdocs_site
git diff --check
```

文档描述当前实现，不把计划写成事实。源码结构只在 `code_overview.md` 维护，算法
总式只在 `project_algorithm_design.md` 维护，未完成项只在 `BUG.md`/`TODO.md` 维护。

## 9. Review

修改后先查 bug，再从第一性原理审视：

1. 是否直接解决原始问题；
2. 是否还有可删除的中间层、复制数组或重复公式；
3. 单位、坐标和数组 owner 是否唯一；
4. 并行 scratch、cache 和 callback 是否仍线程安全；
5. 是否意外扩大 ABI 或 public selector；
6. diff 是否只含本任务文件；
7. 真实运行和物理诊断是否支持结论。

## 10. 提交

只 stage 本任务文件。提交前查看 `git diff --stat`、`git diff --check` 和 status；
完成状态仅记录 Git 提交号。推送前再次确认 remote、branch、HEAD 和工作树。

## 11. 常用真实入口

```bash
TMPDIR=/tmp uv run python build_extensions.py --module Dynamics_forward --module Dynamics_reverse --force
TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_fullhide_1d --force
TMPDIR=/tmp uv run python build_extensions.py --module hadronic_forward_1d --force
```

命令是入口示例，不替代 source closure 审计。模块名、编译顺序和公开 wrapper 以
当前 `build_extensions.py` 为准。

文档-only 变更不编译 Fortran，但仍需 strict MkDocs；代码-only 变更也必须检查受
影响文档是否引用被删除的字段、文件或 selector。

验证失败时修复根因，不降低检查等级。
