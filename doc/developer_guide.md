# 开发指南

本文档记录 ASGARD 当前开发工作流。它补充 `AGENTS.md`，不替代其中的硬性约束。

## 工作区基线

默认工作区：

```text
/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow
```

每次开始工作先执行：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && git status --short --branch && git log -1 --oneline'
```

不要使用旧路径 `~/projects/ASGARD_GRBAfterglow`。

## Git 原则

- 不回滚用户已有改动，除非用户明确要求。
- 不用 `git reset --hard` 或 `git checkout --` 清理未确认改动。
- 不盲目 `git add .`。
- 生成 artifact 前后记录 HEAD、status、diff stat、命令和验收口径。
- 完成任务后提交应只包含本任务文件。

## 代码分层

Python 负责：

- public API
- config / dataclass contracts
- runtime dispatch
- wrapper / orchestration
- benchmark / plotting
- fitting

Fortran 负责：

- dynamics kernels
- electron transport
- radiation integrals
- hadronic microphysics
- observer interpolation kernels

规则：

- 重要数值物理应在 Fortran。
- Python 只做 orchestration、wrapping、benchmark、API glue。
- 不创建大概率无效的代码。
- 不添加 fallback、heuristic post-processing 或非物理 smoothing。

## 修改 Fortran 的流程

1. 明确物理动机和受影响 kernel。
2. 阅读调用链和已有 helper。
3. 做最小改动。
4. 检查声明块长度和 line truncation。
5. 强制编译受影响 extension。
6. 跑最小 smoke/benchmark。
7. 做 review 查 bug，再从第一性原理确认实现是否仍是最简单稳健路径。

声明块规则：

- 同类型声明合并到一行。
- 语义相近的量分组。
- B 类子程序声明块不超过 15 行。
- 禁止每行只声明一个变量。

## 修改 Python 接口的流程

1. 确认 public API 名称和返回类型。
2. 不改变现有兼容入口，除非任务要求。
3. API 边界只验证用户输入、外部 API 或文件输入。
4. 内部状态不做防御性编程。
5. 更新 `doc/public_api.md` 和相关 smoke。

## 物理结果自检

必须检查：

- 时间演化是否连续。
- 频谱是否存在非物理断崖。
- density jump 或 injection event 是否是唯一允许的突变来源。
- `sigma -> 0`、禁用 process、低阶极限是否回到 baseline。
- 图像 artifact 是否由脚本生成。

不能做：

- 通过 smoothing 隐藏错位。
- 通过经验时间因子修正峰时。
- 通过数值 floor 掩盖负值或振荡，除非 floor 是明确物理边界。
- 把 comparison backend 当作物理真值。

## 文档更新规则

改动以下内容时需要同步文档：

- Public API 参数或返回类型。
- 新 solver 或删除 solver。
- 新 physical switch。
- 新 benchmark artifact。
- 新 unsupported boundary。
- 构建命令变化。
- 运行路径或环境变化。

文档入口：

- `doc/index.md`
- `doc/public_api.md`
- `doc/physics_model.md`
- `doc/numerical_methods.md`
- `doc/validation_and_benchmarks.md`

## 最小提交前检查

文档-only：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && git diff --check'
```

Fortran：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module MODULE_NAME --force'
```

然后跑 line-truncation 和最小 smoke。

## Review 检查表

提交前逐项确认：

- 变更是否直接服务原始需求。
- 是否删掉了调试代码、临时文件和失败 artifact。
- 是否无 fallback、无经验补丁、无不必要防御性代码。
- 物理路径是否仍由已有优秀算法或正式 kernel 承担。
- 新文档是否描述当前实现，而非计划。
- 受影响 benchmark 是否可复现。
- `git status --short --branch` 是否只显示预期文件。

## 发布与推送

推送前：

1. 确认目标 remote。
2. 确认 branch。
3. 确认 HEAD 和远端 ref。
4. 确认没有未提交改动。

示例：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && git status --short --branch && git remote -v && git log -1 --oneline'
```
