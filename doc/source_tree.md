# ASGARD Source Tree

本文档只保留当前真正有用的源码入口，不再维护整仓库逐文件树。

代码内容总览见 `doc/code_overview.md`。本文档只作为精简目录索引。

## Python Runtime

- `ASGARD/`
  - 外部用户 API
  - 重点文件：
    - `api_model.py`
    - `api_observe.py`
- `asgard_core/`
  - 内部编排与物理耦合
  - 重点文件：
    - `asgard_runtime.py`
    - `asgard_state.py`
    - `asgard_setup.py`
    - `asgard_ssc.py`
    - `asgard_types.py`
    - hadronic Python wrappers / transitional references:
      - `hadronic_pgamma.py`
      - `hadronic_hummer.py`
      - `hadronic_am3_solver.py`
      - `hadronic_bethe_heitler.py`
      - `hadronic_hadronic_ic.py`
      - `hadronic_pp.py`
      - `hadronic_pair_production.py`
      - `hadronic_species_transport.py`
      - `hadronic_secondary_radiation.py`
      - `hadronic_acceleration.py`
  - 说明：
    - 这些 Python hadronic 模块当前承担编排、benchmark、过渡性微观核参考
    - 按当前项目约束，凡是直接从 AM3/Hümmer/KA2008 转写的最终微观物理实现，都必须迁到 `src/Hadronic/`

## Fortran Runtime

- `src/Constants.f90`
- `src/Dynamics/`
  - `dynamics_common.f90`
  - `Dynamics_forward.f90`
  - `Dynamics_reverse.f90`
- `src/Electron/`
  - 1D/2D 轻子输运核
  - 重点文件：
    - `electron_common.f90`
    - `electron_radiation_kernel.f90`
    - `electron_cooling_kernel.f90`
    - `electron_seed_history_kernel.f90`
    - `electron_transport_2d_kernel.f90`
    - `FS_electron_fullhide_1d.f90`
    - `FS_electron_fullhide_2d.f90`
    - `FS_electron_charint_1d.f90`
    - `FS_electron_charint_2d.f90`
- `src/Radiation/`
  - `radiation_common.f90`
  - `SSC_spec.f90`
  - `Annihilation.f90`
  - `Seed_reverse.f90`
- `src/Hadronic/`
  - 当前正式 hadronic Fortran 入口
  - 重点文件：
    - `hadronic_common.f90`
    - `hadronic_transport_kernel.f90`
    - `hadronic_radiation_kernel.f90`
    - `hadronic_interaction_kernel.f90`
    - `hadronic_decay_kernel.f90`
    - `FS_hadronic_1d.f90`
  - 说明：
    - 后续所有 AM3-derived 最终微观核都应落在这里
    - Python 不再是最终 AM3 转写物理代码的归宿

## Build Entrypoint

- `build_extensions.py`
  - 当前登记的 active hadronic Fortran extension：
    - `FS_hadronic_1d`

## Tests / Benchmarks

- `tests/hadronic_1d_smoke.py`
- `tests/hadronic_proton_synch_1d_diag.py`
- `tests/hadronic_pg_neutrino_1d_diag.py`
- `tests/hadronic_bethe_heitler_smoke.py`
- `tests/hadronic_hadronic_ic_smoke.py`
- `tests/hadronic_pp_smoke.py`
- `tests/hadronic_pair_production_smoke.py`
- `tests/hadronic_pgamma_benchmark_report.py`
- `tests/vegas_afterglow_comparison.py`
