# no-scale 几何保护：稳定模量与量子反馈的深入检验

2026-09-25。承接上一轮隐藏部门传递问题，比较理想 no-scale、超势稳定和四次 Kähler 几何稳定三个明确作用量。**四次几何可以在局部稳定新增模量，同时保留对危险带电分裂的抑制；计算新增重模量的有限量子反馈及一次真空能调节后，仍有条件小反馈示例。**

这个进展限于指定重整化条件下已算出的局部贡献。完整同阶匹配、几何参数的技术自然性、真实物质及新宇宙背景仍待建立。没有新增观测拟合，也没有证明逆对数平方律或现实物理常数随时间变化。

- [中文结果、完整作用与边界](reports/results_cn.md)
- [反馈抑制与真空能调节图](figures/protection_and_radiative_vacuum.png)／[PDF](figures/protection_and_radiative_vacuum.pdf)
- [几何参数预算热图与相位图](figures/geometry_budget_and_phase.png)／[PDF](figures/geometry_budget_and_phase.pdf)
- [协议及主生产前解析补充](protocol.md)／[材料清单](materials_manifest.json)
- [理想 no-scale 推导](../../reports/noscale_exact_20260925.md)
- [超势稳定对照](../../reports/noscale_stabilization_20260925.md)
- [四次几何与真空能调节](../../reports/noscale_kahler_20260925.md)
- [量子接口、来源与技术自然性缺口](../../reports/noscale_quantum_sources_20260925.md)
- [独立高精度数值核查](reports/independent_results_cn.md)

## 结果范围

复用旧1025点shift轴轨迹作为参数截面。873组主扫描、5238组匹配方案诊断以及15375行代表曲线均已归档。没有把旧轨迹冒充新作用量已积分的解，86页正文保持原样。

以 $a=b=1,m_G=1$ eV、固定模量匹配能标为例，将标量圈图真空能对应的几何调节也计入后，已算贡献的力模之和相对于参考时钟力约为 $2.58\times10^{-21}$。在同一诊断中，增大 $m_G$ 会重新放大反馈，10%预算约在0.316 MeV触及；这不是粒子质量实验上限，也不是完整两圈约束。只用树级真空参数时得到约4.08 MeV，两者不同说明真空调节的连带影响不能忽略。

## 复现与文件

在仓库根目录执行，需要 Python、NumPy、SciPy、SymPy、mpmath、Matplotlib，确切版本见材料清单：

```bash
timeout 180 python experiments/noscale_protection_v01/code/run_protection.py
timeout 180 python experiments/noscale_protection_v01/code/plot_results.py
timeout 180 python reports/noscale_exact_20260925_checks.py
timeout 180 python reports/noscale_stabilization_20260925_checks.py
timeout 180 python reports/noscale_kahler_20260925_checks.py
timeout 180 python reports/noscale_quantum_sources_20260925_checks.py
```

命令覆盖对应输出；历史日志保留。高精度独立程序采用分批命令，详见其报告；它不读取或导入主生产实现。主生产66项、理想基准52项、超势稳定70项、四次几何35项及量子接口24项检查通过。独立数值汇总527/527通过，另有108个完整势精确点及54个原始超迹点的两档高精度复算。这些数字是实现核验，不是独立物理证据数。

`scan.csv` 为873组主方案；`scheme_scan.csv` 加三种固定匹配尺度和两种η选择，含5238组；`representative_curves.csv` 在同一行保留 `tree_eta_*` 与 `retuned_*` 两组结果；`budget_limits.csv` 为54个条件预算根。势和对无量纲χ/y的势导数以eV⁴表示，质量以eV、质量平方以eV²表示。`B_abs2` 为双线性分裂的模平方。

`selected_force_ratio` 为已算向量贡献相加后的模；`component_sum_ratio` 为各向量模之和，以避免用偶然抵消定义窗口。后者只上界本轮选出的有限贡献，不上界未知完整有效作用量。`above_threshold_AMSB_*` 只指100GeV带电阈值以上的规范接口，不是现实光子超伙伴预测。
