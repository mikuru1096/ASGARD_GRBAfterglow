# ASGARD ApJ manuscript

本目录保存 AASTeX 稿件、图、source data 与可复现构建脚本。它不是项目 TODO；论文任务在 Git 中追踪。

## 固定论点

ASGARD 的核心论点是：同一套守恒、粒子输运、辐射和等到达时间计算链可以在明确误差与后端边界下，连接 GRB 余辉物理、可观测量和参数推断。正文只陈述已由代码、源数据或文献直接支持的结果。

## 证据规则

- 每幅定量图都必须由已提交脚本和 `paper/source_data/` 重建。
- benchmark 必须记录版本、网格、硬件、线程和统计量，不能只给单次 wall time。
- 物理 claim 对应方程、真实代码入口和验证结果；未闭合功能明确写成边界。
- 不用平滑、经验归一化或未说明的筛点改善图形。
- 正文中的数值必须能在表格或 source data 中逐项找到。

## 冻结 baseline

发布候选稿必须记录 Git commit、依赖锁、构建命令及 figure/source-data checksum。当前代码基线以仓库 HEAD 与 `uv.lock` 为准；论文不得引用工作树未提交结果。已知缺陷见根目录 `BUG.md`。

## 构建

在项目根目录执行：

```bash
cd paper
latexmk -pdf -interaction=nonstopmode main.tex
```

`main.tex`、投稿图和 source data 已跟踪；当前仓库没有统一的 figure-generation 脚本，因此图件尚不能从一个正式入口全量重建。发布前必须补齐该链，不能用手工图片或另一个 TeX 主文件代替。最终检查包括引用、单位、交叉引用、图例、字体嵌入、source-data 对应关系和干净重建。

## 目录约定

- `main.tex`：唯一主稿入口。
- `figures/`：由脚本生成的投稿图。
- `source_data/`：图和表对应的机器可读数据。
- `scripts/`：只包含重建稿件 artifact 的脚本。

临时草稿、profile 和探索图放 `/tmp`，不提交。
