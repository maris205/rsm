# 匹配尺度 PM 产物独立审计

审计检查 440/440 通过；完整：True。

审计通过表示保存结果及失败报告一致，不表示科学门限通过。

| 案例 | 事件 a | 科学门限失败 |
| --- | ---: | --- |
| ref_matched64 | 0.5109537192205208 | PM_vs_ODE_event, PM_vs_ODE_radius |
| ref_matched128 | 0.5015828888584287 | PM_vs_ODE_radius |

独立旧 fullFFT 求值次数：初态 0，终态 2；复用已存初态力向量 2 份，不重跑粒子轨迹。

汇总核验：`{"scientific_case_failures": 3, "scientific_global_failures": ["matched128_event_vs_matched64"], "high_to_low_event_fractional_difference": -0.018339880912086604}`

- Scientific gate failures remain scientific results; artifact audit PASS does not override them.
- Historical drift/force diagnostics are checked as saved records; intermediate particle states are not archived.
- Initial force vectors are reused from hashed continuum-diagnostic artifacts; no new initial force call. Their force moments/Q4 are independently reconstructed, while raw deposition diagnostics are cross-file comparisons.
- One final fullFFT per case is an audit evaluation separate from the initial diagnostic and production trajectory.
- Legacy fullFFT shares the fixed CIC operator and checks discrete execution, not continuum gravitational accuracy.
- Q4 supplements second-moment axes but does not exhaust angular anisotropy or establish three-dimensional convergence.
- Two matched scales and epsilon=0 do not establish phase/time convergence, resolved halo profiles, galaxies, or clock response.

完整数值和哈希见 `independent_artifact_audit.json`。复现：`python code/audit_matched_artifacts.py`。
