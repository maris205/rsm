# 隐藏破缺传递：慢时钟的下一项一致性检验

2026-09-25。本轮在上一阶段的移位型超引力作用中加入一个明确的 nilpotent 隐藏场，计算它对时钟势和重带电谱的影响。**最简单的加性完成尚未通过检验：在本轮扫描和固定有限匹配条件下，带电量子力的小反馈区与计算 100 GeV 阈值所需的名义有效理论能标区不相交。**

一个值得保留的例外是：相位恰好为 π/2 时，横向重场移到局部谷底，树级指数势的形状仍能保留。但带电质量分裂仍在，不能仅靠选相位解决量子反馈。这个结果约束的是指定实现，并非对逆对数变化、所有超引力完成或所有隔离机制的排除。

- [中文结果、公式与适用边界](reports/results_cn.md)
- [受力预算与有效理论能标图](figures/hidden_force_and_domain.png)／[PDF](figures/hidden_force_and_domain.pdf)
- [特殊相位谷底及额外规范传递图](figures/heavy_valley_and_gauge_channel.png)／[PDF](figures/heavy_valley_and_gauge_channel.pdf)
- [执行协议](protocol.md)／[材料与哈希清单](materials_manifest.json)
- [完整作用与解析核查](../../reports/hidden_transfer_theory_20260925.md)／[一手来源](../../reports/hidden_transfer_theory_20260925_sources.json)
- [第二条独立解析路径](../../reports/hidden_transfer_independent_20260925.md)／[独立数值核查](reports/independent_results_cn.md)
- [规范传递分量及其限制](../../reports/gauge_leakage_20260925.md)

## 结果范围

复用 `sugra_shift_v01` 的 1025 点实轴背景，扫描 158 个隐藏参数值和 3 个相位，共 474 组；保存 21525 行代表曲线和 9 个局部重场谷。没有重新积分宇宙背景、拟合观测或生成新的 α(t) 预测，86 页正文保持原样。

内部 10% 力预算对应的隐藏超势参数约为 $m_G\lesssim2.0\times10^{-17}$ eV，而名义 $\sqrt{|f|}\ge100$ GeV 标志对应约 $m_G\gtrsim2.4\times10^{-6}$ eV。前者低于 nilpotent 理论已说明的重粒子计算域，只能视为形式预算；后者只是必要量纲标志，通过它也不等于完成了 UV 理论。这两者均不是观测约束，$m_G$ 也不等同于物理引力微子质量。

## 复现

在仓库根目录执行。需要 Python、NumPy、SciPy、SymPy、mpmath、Matplotlib；版本和输入哈希见材料清单。父实验输入随仓库提供，无需重跑旧宇宙轨迹。

```bash
timeout 180 python experiments/hidden_transfer_v01/code/run_transfer.py
timeout 180 python experiments/hidden_transfer_v01/code/independent_check.py
timeout 120 python reports/hidden_transfer_theory_20260925_checks.py
timeout 120 python reports/hidden_transfer_independent_20260925_checks.py
timeout 60 python reports/gauge_leakage_20260925_checks.py
timeout 180 python experiments/hidden_transfer_v01/code/plot_results.py
```

命令会覆盖对应生成文件；原始日志保留。主检查 35/35、独立数值 2015/2015、两条解析路径 29/29 与 24/24、规范分量核查 9/9 均通过。它们验证实现一致性，不能当作独立物理证据数或外部同行评审。独立数值脚本未读取、导入或调用主生产脚本。

`results/scan.csv` 保存每组的最大力比、最大值位置、谱正性、有效理论域标记及规范分量；`representative_curves.csv` 保存逐点势、势导数和带电谱；`valleys.csv` 保存完整势的局部极小值及渐近误差。质量以 eV、质量平方以 eV²、势与无量纲 χ 的势导数以 eV⁴ 表示。`beta=mG/H_ref`，`B_abs2_eV4` 表示双线性分裂系数的模平方。`VA` 是约定的量纲标志，`gauge_common_*` 只表示部分传递贡献。

下一步应选择一个明确的几何隔离或对称保护完成，同时检查树级交叉项与带电谱的共同软项、双线性分裂是否受抑制；本轮没有预先宣布这样的完成已成立。
