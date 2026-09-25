# 带电阈值的对称性保护与破缺预算

**五条背景轨迹及270点独立高精度势核验已完成。** 全局N=1超对称启发的配对谱可以抵消重场m⁴势，而电磁阈值同号相加；时钟自身F项只留下很小残余。额外共同soft质量仍可破坏慢时钟，最大预定case出现反转。可见区破缺传递、横向稳定及引力完成尚未解决，没有新增观测支持。

- [中文结果与公式](reports/results_cn.md)
- [时钟/α响应图](figures/protected_clock_response.png)／[PDF](figures/protected_clock_response.pdf)
- [保护机制与soft预算图](figures/protection_and_soft_budget.png)／[PDF](figures/protection_and_soft_budget.pdf)
- [预定协议](protocol.md)、[实际输入](results/inputs.json)、[主检查与结果](results/summary.json)、[完整轨迹](results/trajectories.csv)、[理论预算](results/budgets.json)
- [独立高精度报告](reports/independent_precision_cn.md)、[270点基准](results/independent_precision.json)、[主从比较及30点阈值检查](results/independent_comparison.json)
- [解析与原始来源审计](../../reports/protected_threshold_theory_20260925.md)、[来源登记](../../reports/protected_threshold_theory_20260925_sources.json)、[材料身份](sources.json)

本轮采用常χ动能；不混用原固定盒R模型。质量函数和平方指数是输入。100 GeV例仅用于层级诊断；全部谱、热史和观测条件尚待验证。SUSY配对、已有量子阈值公式及Kähler匹配不是本项目的新发现。

仓库根目录复现，依赖Python、NumPy、SciPy、Matplotlib、mpmath、SymPy：

```bash
timeout 180 python experiments/protected_threshold_v01/code/run_protected.py
timeout 180 python experiments/protected_threshold_v01/code/independent_precision.py
python experiments/protected_threshold_v01/code/plot_results.py
python reports/protected_threshold_theory_20260925_checks.py
```

输入取自同仓库已归档的`charged_threshold_v01`，不依赖外部私有文件，不重跑早期背景或PM。重新运行会覆盖本目录相应产物，比较历史运行前请保留副本。正式正文保留原版，本轮先归档为条件理论材料。

首轮输入曾继承父实验的八例清单；已改为本轮实际五例后重新执行，轨迹CSV字节不变。[元数据修订记录](results/metadata_revision.json)和首轮日志/检查均保留；此修订没有改变方程、参数或初态。
