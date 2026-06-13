# 基准刷新协议

本文档固定 ASGARD benchmark 重新生成的最低基线。目标是让每张图、每个 CSV 和每个 comparison report 都能回溯到代码状态、命令、build 状态和物理验收口径。

## 必须记录的元数据

每次刷新 benchmark 前后都记录：

- Git HEAD：`git rev-parse --short HEAD`
- 工作树状态：`git status --short --branch`
- Tracked diff：`git diff --stat`
- 命令原文，包括 `rtk bash -lc 'source ~/.wsl_env && cd "...ASGARD_GRBAfterglow" && ...'`
- 是否重编译 Fortran extension；若 benchmark 依赖 Fortran 改动，记录具体 `build_extensions.py --module ... --force` 命令。
- 输出路径；图像和 CSV 必须由脚本复现，不手工编辑。
- 验收口径：物理平滑性、过程开关、comparison backend 的角色，以及不作为失败条件的模型差异。

## 标准命令

Vegas baseline full comparison：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run --extra compare python tests/vegas_afterglow_comparison.py --scenario baseline'
```

该命令运行 `tests/vegas_afterglow_comparison.py` 的全部 builder。当前 baseline 会刷新 `output/asgard_doc/vegas_afterglow_compare/` 下由脚本生成的 comparison figures；`magnetic_decay_2d` 会额外输出 broadband spectrum 图。

RS / Vegas comparison：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run --extra compare python tests/vegas_afterglow_comparison.py --scenario baseline --only reverse_shock_lc reverse_shock_thermal'
```

Magnetized RS sigma scan：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python scripts/benchmarks/magnetized_rs_sigma_benchmark.py --medium ism --mode quick'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python scripts/benchmarks/magnetized_rs_sigma_benchmark.py --medium wind --mode quick'
```

当前 magnetized RS scan 使用同一组爆发参数以让射电和光学 RS 都显著：`E_iso=1e54 erg`, `Gamma0=50`, `duration=500 s`, `rvs_eps_e=0.3`, `rvs_eps_B=0.1`。纯 ISM 使用 `n_ism=0.1 cm^-3`；纯 wind 使用 `Astar=0.05` 且 `n_ism=1e-15 cm^-3` 作为 ASGARD/Vegas 中的 ISM floor，不加入 wind-to-ISM 转变。

Magnetized RS ASGARD/Vegas light-curve overlay：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python scripts/benchmarks/magnetized_rs_vegas_lc_compare.py --medium ism'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python scripts/benchmarks/magnetized_rs_vegas_lc_compare.py --medium wind'
```

这两条命令输出：

- `output/asgard_doc/magnetized_rs_sigma_benchmark/magnetized_rs_sigma_ism_vegas_lc_compare.png`
- `output/asgard_doc/magnetized_rs_sigma_benchmark_wind/magnetized_rs_sigma_wind_vegas_lc_compare.png`

Triple density-jump RS+FS top-hat light curves：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python scripts/benchmarks/triple_density_jump_rs_fs_tophat.py --mode formal'
```

该命令复用 ASGARD/Vegas baseline 的典型 top-hat 余辉参数和网格，固定 `E_iso=1e52 erg`、`Gamma0=100`、RS upstream `sigma=0.1`、`theta_c=0.1`、`n_ism=1 cm^-3`、`z=0.1`、`d_L=1e26 cm`；开启 RS+FS synchrotron，关闭 SSC/IC cooling。密度跳变半径为 `1e14, 1e15, 1e16 cm`，半径空间 normal Gaussian 相对标准差 `sigma_R/R_jump=0.1`，增强因子 `100`；初始半径前移到 `1e13 cm`，保证第一个宽跳变在计算域内完整展开。输出为 `output/asgard_doc/reverse_density_jump_tests/triple_density_jump_rs_fs_tophat.png/.svg/.pdf/.tiff`，图中包含密度剖面以及 1 GHz、`1e14 Hz`、`1e18 Hz` 光变，对比 no-jump total 与 triple-jump total/FS/RS。

VegasAfterglow 2.0.5 的 stock sdist 在纯 ISM、`sigma >= 1e-5` 时会给出全 0 reverse-shock synchrotron。刷新 overlay 前先对 2.0.5 sdist 应用本仓库记录的 patch，并把 patched wheel 安装进当前 uv 环境：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/_external/vegasafterglow_205_src/vegasafterglow-2.0.5" && git apply --unidiff-zero "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow/doc/patches/vegasafterglow-2.0.5-rvs-ism-compression.patch"'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv pip install --force-reinstall --no-deps "/mnt/c/Users/jia/Documents/New project/_external/vegasafterglow_205_src/vegasafterglow-2.0.5"'
```

如果本地 Vegas source tree 已经应用过该 patch，第一条命令会因为补丁已存在而失败；此时应检查 `src/dynamics/reverse-shock.tpp` 中 `compression_advance <= 0` 边界和 `src/config/simulation-defaults.h` 中 `max_ode_steps = 1000000` 是否存在，而不是重复打补丁。

Lan 2023 polarization overlay：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/polarization_literature_overlay.py'
```

2D chi-resolved EATS benchmark：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python scripts/benchmarks/chi_eats_2d_benchmark.py --mode quick --solver fullhide_2d --medium both --only all'
```

该命令固定 `num_chi=24`，全图使用 `num_r=300`、`num_theta=300`、偏轴 `num_phi=50`；对轴 `theta_v=0` 使用 `num_phi=1`。观测时间覆盖 `1e2-1e9 s`，SED 频段覆盖 `1e6-1e28 Hz`。典型余辉参数为 `E_iso=1e53 erg`、`Gamma0=300`、`epsilon_e=0.1`、`epsilon_B=1e-3`、`epsilon_B_floor=epsilon_B`、`magnetic_decay_alpha_t=0`、`p=2.3`、`theta_j=0.1`、`z=0.1`；该设置要求 2D 壳层内部磁场在每个 shell 内为常数，不使用 χ 方向的特殊磁场剖面。对 `fullhide_2d` 在 ISM (`n=1 cm^-3`) 和纯 wind (`A*=0.1`, `n_ism=1e-15 cm^-3`) 下生成 `sed_legacy` 和 `chi_eats_2d` 的总通量对比图；其中总通量包含现有 shell-level forward SSC，`chi_eats_2d` 一期只替换 FS synchrotron+SSA observer projection。输出位于 `output/asgard_doc/chi_eats_2d_benchmark/`，包括 `chi_eats_2d` vs `sed_legacy` light curve/SED、χ 几何/辐射诊断、2D/1D SED、以及 top-hat `theta_v/theta_j = 0, 0.5, 1, 1.5, 3, 5` light-curve 对比 PNG/PDF；数值摘要写入 `chi_eats_2d_summary.csv` 和 `chi_eats_2d_metadata.json`。SSC、hadronic、pair cascade 仍按 shell-level contract，不能解读为 chi-resolved 闭环。偏轴 EATS 禁止使用 `num_phi=1`，因为 `num_phi=1` 是仅适用于 `theta_v=0` 的轴对称 φ 折叠。

ISM χ grid convergence scan：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python scripts/benchmarks/chi_eats_2d_benchmark.py --mode quick --solver fullhide_2d --medium ism --only chi-grid-scan'
```

该命令固定 ISM 与 `fullhide_2d + chi_eats_2d`，使用 `theta_v/theta_j=0.5`，除 χ 以外的网格相对 quick benchmark 减半：`num_gam_e=16`、`num_nu=21`、`num_r=150`、`num_theta=150`、`num_phi=25`、`num_tobs=64`；扫描 `num_chi=32,64,128,256,512`，图中以 `num_chi=512` 为参考显示三个频段的 light-curve 收敛和 wall time。runtime 计时前先执行一次不入图的 `num_chi=32` warm-up，以排除首次调用的冷启动开销。输出为 `fullhide_2d_ism_chi_grid_scan.png/pdf`、`chi_grid_scan_ism_summary.csv`、`chi_grid_scan_ism_lightcurves.csv` 和 `chi_grid_scan_ism_metadata.json`。

## 构建门槛

文档或 plotting-only 脚本改动只需：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && git diff --check'
```

Fortran 或物理路径改动必须按受影响范围执行：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module Dynamics_forward --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module Dynamics_reverse --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_charint_1d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_charint_2d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_fullhide_1d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_slc1_1d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module hadronic_forward_1d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_reverse_kernel --force'
```

Fortran line-truncation 检查使用 WSL gfortran，并在 touched source group 上启用 `-Wline-truncation`。如果相关 extension 已改动但未重建，不应把 benchmark 视为已刷新。

## 物理验收

- RS 图必须说明 ASGARD 使用局部 `gamma34` 注入和显式 `U3/V3` thermal state。VegasAfterglow 是 comparison backend，不是目标答案。
- 如果 ASGARD 与 Vegas 的差异来自 Vegas 使用 averaged `Gamma_th` 而 ASGARD 使用局部 shock-front injection，应记录为模型差异。
- Magnetized RS sigma scan 必须检查 `B3`, `gamma34`, `U3/V3`, `nu_m`, `nu_c`, `nu_a` 的平滑性和 `sigma -> 0` 极限。
- Polarization overlay 必须分开报告 peak time 和 peak amplitude。不要用经验 time shift 或 smoothing 修正 timing mismatch。
- Hadronic report 必须写明启用过程，以及路径是 formal 1D shell transport 还是已记录边界。
- 2D chi-resolved EATS 图必须使用 `num_chi=24`、`num_r=300`、`num_theta=300`、偏轴 `num_phi=50`、对轴 `num_phi=1`、`t_obs=1e2-1e9 s`、SED `1e6-1e28 Hz`，分别检查 ISM 和纯 wind 环境下的 `chi_eats_2d / sed_legacy`、`2d / 1d` SED、以及不同 `theta_v/theta_j` 下 top-hat 光变 ratio 连续平滑；χ 诊断应显示 projection 网格随当前 shell 的正半径 BM 壳层域自适应，晚期不得出现负半径、全部 `chi_dvolume_weight` 同时归零或源项截断导致的断崖。射电 SSA 敏感测试应检查 transport-to-projection χ remap 的 `sum(P*Delta chi)` 与 `sum(tau)` 守恒，以及 emitting cell 使用 τ 坐标平均 escape probability 后不再由 photosphere 落在哪个 χ cell 边界决定。图中 SSC 只是 shell-level forward SSC 对总通量的贡献，不属于 chi-resolved EATS 验收项。
- ISM χ grid convergence scan 必须检查 `num_chi=32..512` 的 `F_chi/F_512` 曲线随时间连续，wall time 随 `num_chi` 递增趋势合理；CSV 保留完整 ratio min/median/max，不用显示范围裁剪替代原始数据。

## 产物策略

- 提交用于复现 benchmark 的脚本和文档。
- 只有当 benchmark figures 是预期文档资产并可由记录命令复现时才提交。
- 不提交 `.buildcache/`、临时 debug 脚本、partial CSV、失败占位图或一次性 notebook export。
- Untracked benchmark output directories 可以保留在本地，直到有意提升为文档资产。
