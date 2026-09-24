# PM 保存产物独立审计

审计检查：876/876 通过；完整性：全六例。

这里的通过表示保存产物、独立重建和失败报告一致。科学精度门限未通过不被改写为通过，也不因此自动视为审计错误。

| 案例 | 终点 a | Δ200 事件 a | 科学门限失败 | 完整终态 |
| --- | ---: | ---: | --- | --- |
| ref_force128 | 0.5326680358 | 0.5326680358 | nonhomology | 无 |
| ref_main | 0.5324377348 | 0.5324377348 | nonhomology, PM_vs_ODE_event, PM_vs_ODE_radius | 有 |
| ref_halfstep | 0.53242186 | 0.53242186 | nonhomology | 无 |
| ref_particles128 | 0.5147345794 | 0.5147345794 | nonhomology | 无 |
| ref_shift | 0.5114809715 | 0.5114809715 | nonhomology, axis_deviation | 无 |
| clock_main | 0.5325213869 | 0.5325213869 | nonhomology | 有 |

独立重建包括初始粒子标签、共同动量、核心半径分位数、奇异值轴比、实际包围计数、固定窗口径向剖面及完整深度二维 CIC 投影。时间顺序、保存步长的宇宙时间积分、聚合极值、阈值事件与所有失败列表逐项核对。

- No particle trajectory was reintegrated; historical drifts/forces use saved diagnostics only.
- Only ref_main and clock_main archive full final states; other cases permit aggregate checks only.
- Final fullFFT uses the pinned legacy CIC operator and an independent complex FFT path; it is not an independent continuum gravitational solver.
- Radial profile window is 0<=r<=.45 box, initial-label shells 0<=q<.4; projection sums the entire periodic z depth.
- Saved projection_axis coordinates (j+.5)/128 are display pixel centres; physical CIC nodes are j/128. No physical phase inference uses display axes.
- Equal second-moment axes do not exclude higher-order cubic anisotropy; regular cubic symmetry can preserve axis ratio one despite nonhomologous radial collapse.
- Implementation/statistic audit success does not imply resolved halos, converged collapse, virialization, stars, or a resolved clock response.

全部数值、源哈希和历史审计失败保存在 `results/pm_artifact_validation.json`。

复现：`python code/validate_pm_artifacts.py`。运行中的可用子集可用 `--available`；`--skip-fullfft` 会明确记录未做终态力重算。
