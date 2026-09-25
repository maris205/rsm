# 最小 U(1) 一圈匹配与条件常数响应

2026-09-25。沿用上一阶段的3条轴向背景，补齐指定重粒子谱的 KL/Konishi 规范匹配，并明确两种UV边界。**固定全纯高能耦合时，物理质量中的显式 Kähler 项被抵消；标准型的局部−22.023599 ppm读数变为匹配后的+23.116249 ppm，移位型约+17.877965 ppm保留。** 固定物理高能耦合是另一输入，仍给旧局部阈值曲线。

这是没有动力学轻带电物质的最小U(1)有效模型，记耦合为α_U。未加入电子／标准模型，不能把结果直接认作实际精细结构常数的预测。没有新增观测、参数拟合或背景积分，86页正文保持原样。

- [中文结果、推导与边界](reports/results_cn.md)
- [匹配结果及组成项图](figures/matched_gauge_response.png)／[PDF](figures/matched_gauge_response.pdf)
- [逆对数响应的近似边界图](figures/response_shape_limits.png)／[PDF](figures/response_shape_limits.pdf)
- [执行协议](protocol.md)／[材料与哈希清单](materials_manifest.json)
- [理论核查](../../reports/gauge_matching_theory_20260925.md)／[来源](../../reports/gauge_matching_theory_20260925_sources.json)
- [独立推导及18点高精度核查](../../reports/gauge_matching_independent_20260925.md)
- [独立3075点复算](reports/independent_results_cn.md)

## 复现

在仓库根目录执行，需要Python、NumPy、SymPy、mpmath、Matplotlib；确切版本见材料清单。父实验已有CSV和输入随仓库提供，无需重跑宇宙背景。

```bash
timeout 180 python experiments/gauge_matching_v01/code/run_matching.py
timeout 180 python experiments/gauge_matching_v01/code/independent_check.py
python reports/gauge_matching_theory_20260925_checks.py
python reports/gauge_matching_independent_20260925_checking.py
python experiments/gauge_matching_v01/code/plot_results.py
```

命令会覆盖对应生成文件；旧生产日志与协议作为历史记录保留。独立脚本不调用主匹配函数。主程序38项、理论22项、探索性独立125项和全轨迹独立149项检查全部通过；独立全曲线按峰值归一的最大差为 $4.01\times10^{-14}$。检查数不是独立物理证据数，完整定义以报告和JSON为准。

## 文件和单位

`results/matched_curves.csv` 保存3×1025行，模型参数和背景直接来自 `../sugra_shift_v01/results/`。`*_delta_g_inverse2` 是无量纲逆规范耦合的各项，乘4π得到逆α_U变化；`*_delta_alpha` 为相对今天的无量纲耦合变化，乘10⁶为ppm。`log_R` 是极小谱分裂的稳定对数，`split_delta_alpha` 是保留分裂相对简并谱的代数差。`shape_*` 仅是初端1、今天0的形状诊断，不是幅度拟合。

`results/summary.json` 保存端点、当前漂移（每年）、近似误差、负对照及全部主检查。`results/run_attempt_001.log` 保留首次语法错误，第二次生产日志为 `run_attempt_002.log`。`code/plot_results.py` 仅读取已归档结果，不改变或拟合曲线。来源PDF缓存位于忽略的build目录，版本、读取范围及哈希保存于来源记录。

本轮完成最小谱的局部一圈匹配；现实破缺传递、真实带电物质、完整重力有效作用量及质量律的算术来源仍待完成。
