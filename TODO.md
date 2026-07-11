# ASGARD TODO

这里只记录尚未开始或正在推进的工作；已确认缺陷放在 `BUG.md`，完成项只以 Git 提交号表示。

## 当前工作

- 沿 Fortran 计算链继续审计：Dynamics → Electron → Radiation/EATS → structured/2D → Hadronic。
- 从 Fan & Piran 原文重推 KN 修正的 Y 参数，闭合自洽冷却后再开放 Python selector。
- 收敛相互绑定的 Python/Fortran selector，公开层只保留物理上可独立选择的组合。
- 为 detailed pp gamma 模型建立连续性、能量预算和性能基准。

## 长期边界

- 完成 chi-local photon/hadron density、secondary feedback 和 observer projection 后，才扩展 2D hadronic。
- 闭合 γγ/电磁级联的几何、单位和状态所有权。
- 扩展 reverse-shock BH electron grid 的运动学覆盖。
- 公开自定义介质、非 `k=2` wind 和 spreading 前，先完成真实 kernel 契约。

## 不做

- 不以 fallback、clamp、smoothing 或经验归一化掩盖物理或数值错误。
- 不为没有决策价值的组合穷举测试。
