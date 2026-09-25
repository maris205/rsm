# 带电质量阈值与慢时钟：最小量子反馈实验

**已完成8条预定背景轨迹及独立复算。** 指定的χ⁻²质量项可以传递为低能α的主阶响应；强量子反馈会使时钟反转，而较重带电质量产生的未保护势能尺度远大于慢变宇宙背景。这是一份有明确适用条件的理论试验，没有新的观测拟合，也没有从黎曼零点推导出带电粒子。

- [中文结果与推导](reports/results_cn.md)：作用量、全部结果、形状区别、辐射尺度和后续理论任务。
- [图1：时钟与电磁响应](figures/clock_and_threshold.png)／[矢量PDF](figures/clock_and_threshold.pdf)。形式低能响应；示例粒子质量不满足光学重阈值条件。
- [图2：反馈力与势能预算](figures/feedback_and_budget.png)／[矢量PDF](figures/feedback_and_budget.pdf)。预算是指定方案的诊断，不是实验排除图。
- [生产前协议](protocol.md)、[输入](results/inputs.json)、[主结果及67项检查](results/summary.json)、[完整轨迹CSV](results/trajectories.csv)。
- [独立宇宙时间/动量复核](reports/independent_audit.md)：8/8案例通过，最大归一差异5.121×10⁻¹⁰；[原始检查](results/independent_checks.json)。
- [作用量与一圈系数审计](../../reports/charged_threshold_action_audit_20260925.md)：10项代数检查，带电阈值、有效势、动能与源项符号。
- [来源和范围登记](sources.json)。

比较的是两种不同动能完成（常χ与常R），各有树级probe和三种反馈强度。probe关闭量子反馈，是形式对照；有限反馈示例为约0.79—7.89 meV的假设带电阈值，尚未证明满足粒子实验与热史约束。两条强反馈轨迹的反转作为结果完整保留。

在仓库根目录复现，需要Python、NumPy、SciPy、Matplotlib（解析检查另需SymPy）：

```bash
timeout 180 python experiments/charged_threshold_v01/code/run_threshold.py
timeout 180 python experiments/charged_threshold_v01/code/independent_check.py
python experiments/charged_threshold_v01/code/plot_results.py
python reports/charged_threshold_action_audit_20260925_checks.py
```

主求解器只复用同仓库旧背景的初始参考，不依赖相邻项目，不启动PM生产任务。已归档的日志位于`results/`；重新运行会重写相应结果，比较历史运行前请保留副本。协议与代码哈希随每次运行写入JSON。结果报告与来源的支持范围仍需由作者审阅，不把数值通过当成观测支持。
