# 故障排查

本文列出 ASGARD 文档和常用运行路径的高频问题。

## 1. 找不到 Fortran 扩展

现象：导入或运行 `Model.flux_density_grid` 时提示缺少 `.so`。

处理：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_fullhide_1d --force'
```

若启用了反向激波或强子路径，按 `doc/validation_and_benchmarks.md` 构建对应模块。

## 2. MkDocs 构建失败

先运行：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run --with "mkdocs<2" --with "mkdocs-material>=9.5" mkdocs build --strict'
```

常见原因：

- `mkdocs.yml` 的 `nav` 引用了不存在或未提交的文件。
- Markdown 表格或代码块没有闭合。
- LaTeX 公式分隔符不成对。

## 3. LaTeX 公式检查失败

公式必须使用 `\(...\)` 或 `\[...\]`。不要把纯文本伪公式写成代码块。若公式需要多行对齐，使用

\[
\begin{aligned}
a &= b+c,\\
d &= e+f.
\end{aligned}
\]

如果公式在网页中能显示但 LaTeX 编译失败，优先检查下标、反斜杠和 `aligned` 环境。

## 4. 光变不平滑

真实物理量随时间或半径演化应连续平滑。例外必须来自明确物理事件，例如 density jump、energy injection 或 shock crossing。

排查顺序：

1. 检查动力学轨道 \(\Gamma(R)\)、\(R(t_{\rm obs})\)。
2. 检查 \(B(R)\)、\(\nu_m\)、\(\nu_c\)、\(\nu_a\)。
3. 检查电子或质子谱是否相邻壳层突变。
4. 检查 observer projection 的时间/频率网格是否覆盖观测点。
5. 检查是否误把 unsupported backend 当成已支持能力。

不要用 smoothing、经验时间因子或后处理裁剪掩盖问题。

## 5. 拟合不收敛

先确认：

- 数据单位是否为秒、Hz 和 cgs flux density。
- `flux_err` 是否非零且代表真实观测误差。
- quick grid 下 best-fit 光变是否已能大致穿过数据。
- 参数范围是否过宽或把物理上不相关的开关同时打开。

拟合应从最小正向激波模型开始。只有 residual 支持时，才加入 wind、结构化喷流、反向激波或强子过程。

## 6. 私有网页不可见或公开

Pages 访问权限由 GitHub 仓库设置控制，不在源码中保存。检查：

1. `mikuru1096/ASGARD_private -> Settings -> Pages` 的 Source 是否为 `GitHub Actions`。
2. Pages visibility 是否为 `Private`。
3. 合作者是否在仓库 collaborators/teams 中。

若 Pages 设置页没有 `Private` 选项，说明当前仓库或组织没有可用的私有 Pages 能力。
