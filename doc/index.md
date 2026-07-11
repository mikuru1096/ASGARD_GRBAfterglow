# ASGARD 文档入口

ASGARD 用 Fortran 数值核计算 GRB 余辉动力学、粒子输运和辐射，用 Python 组织物理配置、观测投影与参数推断。

## 第一次使用

1. [安装与构建](installation.md)
2. [快速开始](quickstart.md)
3. [示例](examples.md)
4. [公开 API](public_api.md)
5. [拟合工作流](fitting_workflow.md)

## 理解物理与算法

- [物理模型](physics_model.md)
- [数值方法](numerical_methods.md)与[算法设计](project_algorithm_design.md)
- [电子求解器](electron_solver_algorithms.md)
- [公开后端边界](public_backend_limits.md)

专题文档覆盖 2D q-shell、磁化反向激波、joint secondary feedback、pγ 与 prompt internal shock。专题入口只代表对应研究路径，不自动表示 public API 的所有组合都已验收。

## 开发者

- [代码总览](code_overview.md)
- [开发指南](developer_guide.md)
- [术语表](terminology.md)

## 事实来源

公开签名以 `asgard_core/__init__.py` 和 `asgard_core/api_model.py` 为准；构建模块以 `build_extensions.py` 为准；未完成工作与已知缺陷分别只维护在根目录 `TODO.md` 和 `BUG.md`。文档与代码冲突时先修文档或实现，不增加兼容叙述掩盖差异。

所有公开量默认使用 cgs；时间 s、频率 Hz、距离 cm、角度 rad、能量 erg。物理演化应连续平滑，出现无物理来源的跳变先按 bug 检查。
