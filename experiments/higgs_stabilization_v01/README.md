# Higgs剩余模稳定与局域端点匹配

本阶段检验上一轮Higgs链中最后一个未稳定手征模。**原场内容的相位提升算符次数约为 $5.90\times10^{60}$；新增规范群能产生全重自由谱，却使真实局域端点的最低阶交叉接触精确为零。** 有限辅助场也未在明确限定的刚性时钟关闭分支中恢复该势。候选的限制与可能的其他完成方式分别保留；没有新观测或完整超引力结论。

- [中文结果报告](reports/results_cn.md)：推导、质量表、源的区别、有限辅助场及下一项任务。
- [四联图PNG](figures/stabilization_tradeoff.png)／[可导出PDF](figures/stabilization_tradeoff.pdf)。
- [计算协议](protocol.md)、[主结果](results/summary.json)、[独立比较](reports/independent_results_cn.md)。
- [材料及哈希清单](materials_manifest.json)、[归档审计](results/archive_audit.json)。
- 独立专题：[相位算符](../../reports/higgs_phase_operators_20260926.md)、[全部自由重谱](../../reports/higgs_extra_gauge_20260926.md)、[局域源](../../reports/higgs_local_currents_20260926.md)、[有限辅助场](../../reports/higgs_finite_auxiliary_20260926.md)。

从仓库根目录复现，依赖Python、NumPy、SciPy、SymPy、mpmath和Matplotlib：

```bash
OPENBLAS_NUM_THREADS=1 python experiments/higgs_stabilization_v01/code/run_stabilization.py
OPENBLAS_NUM_THREADS=1 python experiments/higgs_stabilization_v01/code/plot_results.py
OPENBLAS_NUM_THREADS=1 python experiments/higgs_stabilization_v01/code/independent_check.py
python experiments/higgs_stabilization_v01/code/audit_archive.py
```

专题检查脚本与机器结果在仓库 `reports/higgs_*_20260926*`。审计检查归档一致性，不代表物理验证；本阶段没有改动历史实验和86页正式稿。旧输入在清单中固定，已发表材料的通读状态和引用阅读范围不由内部AI检查代替。
