# ASGARD 当前冻结主线

## 1. 当前唯一基线

- 当前工作树就是主线基线，不再继续沿用 Claude 的 Phase 1/2/3/4 计划口径。
- 若 `AGENTS.md` 与其他说明文件冲突，以 `AGENTS.md` 的“2026-04-17 当前冻结基线”为准。

## 2. 当前只做的事

1. 保持 `run_fit(config)`、`Fitter.loglike()` 和 compiled fast-path 可运行。
2. 保持当前构建链可复现：
   - `build_extensions.py`
   - Windows ordered-object fallback
   - `electron_get_y` 构建链
3. 保持 benchmark / plot / demo 入口可运行。
4. 若后续发现物理结果不连续或不平滑，优先回查：
   - `Dynamics_reverse.f90`
   - `calling_modules.f90`
   - `electron_get_y.f90`
   - observer-side 插值与 frame 口径

## 3. 当前不做的事

- 不继续推进单独的“Phase 3 API cleanup”计划。
- 不继续扩张新的抽象层，只接受能够替代现有重复代码的整理。
- 不新增新的 profile/check/doc 辅助脚本。
- 不把阶段性重构说明继续当作项目主计划。

## 4. 维护顺序

1. 先更新 `AGENTS.md`。
2. 再更新 `plan.md`。
3. 再做必要代码与构建修正。
4. 最后只保留与当前主线一致的说明文件。
