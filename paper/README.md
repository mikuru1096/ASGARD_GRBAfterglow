# ASGARD ApJ manuscript

本目录保存 AASTeX 稿件、图、source data 与可复现构建脚本。它不是项目 TODO；论文任务在 Git 中追踪。

## 固定论点

ASGARD 的核心论点是：同一套守恒、粒子输运、辐射和等到达时间计算链可以在明确误差与后端边界下，连接 GRB 余辉物理、可观测量和参数推断。正文只陈述已由代码、源数据或文献直接支持的结果。

## 证据规则

- 每幅定量图都必须由已提交的 source benchmark 入口生成 `paper/source_data/`，再由正式绘图入口只读重建。
- benchmark 必须记录版本、网格、硬件、线程和统计量，不能只给单次 wall time。
- 物理 claim 对应方程、真实代码入口和验证结果；未闭合功能明确写成边界。
- 不用平滑、经验归一化或未说明的筛点改善图形。
- 正文中的数值必须能在表格或 source data 中逐项找到。

## 冻结 baseline

发布候选稿必须记录 Git commit、依赖锁、构建命令及 figure/source-data checksum。当前代码基线以仓库 HEAD 与 `uv.lock` 为准；论文不得引用工作树未提交结果。已知缺陷见根目录 `BUG.md`。

## 构建

在项目根目录执行：

```bash
rtk wsl.exe -d Ubuntu -- bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/generate_cross_code_benchmarks.py'
rtk wsl.exe -d Ubuntu -- bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/generate_paper_figures.py'
rtk wsl.exe -d Ubuntu -- bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && latexmk -cd -pdf -interaction=nonstopmode -halt-on-error -outdir=build paper/main.tex'
```

`tests/generate_cross_code_benchmarks.py`、`tests/vegas_afterglow_comparison.py`
及各专门 benchmark 脚本是 source-data 入口；它们记录物理假设、版本、网格和
运行统计量。`tests/generate_paper_figures.py` 是唯一正式只读绘图入口，读取已跟踪的
`paper/source_data/` 并生成正文 17 张图，不在绘图阶段制造或改写数值证据。正式图同时导出
SVG、PDF 和 600-dpi TIFF；禁止以手工图片、探索性 `output/` 产物或另一个
TeX 主文件替代该链。最终检查包括引用、单位、交叉引用、图例、字体嵌入、
source-data 对应关系和干净重建。

Fig. 1 的浅色 GRB 场景是唯一生成式视觉资产，只承担非定量背景。选定原图保存在
`paper/figures/assets/fig1_afterglow_source.png`，完整 prompt、模式和用途边界保存在
`paper/source_data/fig1_visual_overview_imagegen.json`；所有文字、箭头和壳层次序仍由
正式 Python 入口叠加。其余图件只使用已跟踪数值 source data。

表格介质的 wind-termination、Gaussian jump 与 early PION CSM 响应由同一公共
`TabulatedMedium + reverse-shock` 入口产生。Fig. 10--11 是受控密度结构，
Fig. 12 使用 1869.597 AD 的 PION 1024-cell 剖面及解析自由风内延。正式 artifact
保存在 `paper/source_data/benchmarks/density_structure/`。Observer 光变只区分
total、FS 与 primary RS；secondary RS 仅由 event、branch internal energy、耗散和
局部同步谱光度诊断，不声称存在独立 observer secondary-RS flux。

## 目录约定

- `main.tex`：唯一主稿入口。
- `figures/`：由脚本生成的投稿图。
- `source_data/`：由正式 benchmark 入口生成、与图和表逐项对应的机器可读数据。
- `tests/generate_paper_figures.py`：唯一正式只读图件生成入口。
- `scripts/`：只包含生成 source data 或重建其他稿件 artifact 的脚本。

临时草稿、profile 和探索图放 `/tmp`，不提交。


## Circumstellar-medium source data

`source_data/csm/` contains the 1869.597 AD, 1024-cell PION Eta Car profile, the analytic free-wind extension from $10^{13}$ to $1.32\times10^{16}$ cm, and the deterministic 96-point ASGARD interface. `PION_PROVENANCE_1870.txt` records the run and reduction; `SHA256SUMS` covers the tracked inputs and derived tables. No 1024/2048 convergence claim is made.
