# 指定超引力部门、横向扰动与物理质量阈值

2026-09-25。接续[保护阈值实验](../protected_threshold_v01/README.md)，完成3个固定模型的9条双场背景、12组横向模（各2基解）及物理带电阈值。**移位型候选保留慢背景，横向质量为正但均匀初始位移只衰减0.01294%；标准型候选的物理质量修正使所算局部电磁质量阈值贡献方向翻转。** 完整规范耦合的同圈阶异常匹配尚未完成，不能把该贡献的符号当作完整超引力预测。100 GeV是全纯质量的演示参数，没有观测拟合。

- [中文结果与推导说明](reports/results_cn.md)
- [背景及电磁阈值图](figures/background_and_threshold.png)／[PDF](figures/background_and_threshold.pdf)
- [横向质量、位移及能量图](figures/transverse_modes.png)／[PDF](figures/transverse_modes.pdf)
- [冻结协议](protocol.md)与[材料、哈希及命令清单](materials_manifest.json)
- [理论核查](../../reports/sugra_shift_theory_20260925.md)与[原始文献范围](../../reports/sugra_shift_theory_20260925_sources.json)
- [独立背景复算](reports/independent_results_cn.md)／[方法](reports/independent_methods_cn.md)
- [320/400位粒子核验](reports/particle_precision_cn.md)

## 复现

在本仓库根目录执行，需要 Python、NumPy、SciPy、SymPy、mpmath 和 Matplotlib；已执行环境版本在清单及结果JSON中。父实验的已归档输入和参考轨迹随仓库提供，无须重跑父实验或下载观测数据。下面命令会覆盖本目录对应的生成结果；如需保留原运行时间／哈希，请先保存原产物。

```bash
timeout 180 python experiments/sugra_shift_v01/code/run_sugra.py
timeout 180 python experiments/sugra_shift_v01/code/independent_check.py
timeout 180 python experiments/sugra_shift_v01/code/check_particle_precision.py
python reports/sugra_shift_theory_20260925_checks.py
python experiments/sugra_shift_v01/code/plot_results.py
```

主检查149/149、理论符号23/23、独立背景／势435/435、粒子高精度43/43通过。主解与独立解分别使用 $N=\ln a$ 和宇宙时间，后者不调用主背景求解器。粒子独立核查直接计算超迹和完整导数，覆盖引力共同分裂的 $\chi$ 依赖。

## 结果文件与单位

| 文件 | 内容 |
|---|---|
| `results/inputs.json` | 冻结输入及父实验元数据 |
| `results/summary.json` | 9背景、12组模摘要及149项检查 |
| `results/trajectories.csv` | 每背景1025点，共9225行；阈值只在轴向计算，其他行对应列留空 |
| `results/linear_modes.csv` | 每模组1025点，共12300行；同时含两个基解及模能量 |
| `results/independent_checks.json` | 435项检查、方法、输入／代码／主CSV哈希 |
| `results/independent_*.csv` | 独立宇宙时间解在共同网格的值 |
| `results/particle_precision.json` | 18点高精度谱、势力及阈值核验 |
| `results/*log` | 正式运行日志及首次报告输出失败记录 |

背景列 `V,U,Vx,Vyy` 等以 $F_\chi^2H_{\rm ref}^2$ 为势单位；`u=chi-chi_i`，`qx=dchi/dN`，`wx=dot(chi)/H_ref`，`tau=H_ref*t`，`E=H/H_ref`，`f=dlnH/dN`。`s,d` 为 eV²，`h2=|B|²` 为 eV⁴；`U1,U1p` 为物理 eV⁴ 及其对无量纲 $\chi$ 的导数，不能与未恢复单位的 `Vx` 直接相除。`delta_alpha` 是仅保留局部质量阈值时的无量纲 $\alpha(z)/\alpha(0)-1$ 代理量，图中记为 $\Delta_{\rm th}$，不含完整超引力异常匹配。模文件的 `omega2` 单位为 $H_{\rm ref}^2$，能量单位相应为 $H_{\rm ref}^2$（省略共同正因子），`D0_N,D1_N` 为 $N$ 导数。

`inputs.json` 继承的 `b_em=1/3` 未被本轮使用；实际两复标量加一Dirac粒子的权重为 $(1/3,1/3,4/3)$。旧 `reference_code_sha256` 属于父输入沿革，本轮代码哈希存于 `summary.json` 和材料清单。

## 适用范围

这是超引力导出的指定标量部门在外部辐射、物质、Λ背景中的试验。有限局部重粒子一圈势只在树级轨迹上作预算；没有积分全量子反馈，没有完成隐藏／可见区破缺、全部引力／费米扰动或完整物质来源。$\chi^{-2}$ 质量律仍是输入，$\alpha$ 的完整阈值不是自动精确的 $1/\ln^2t$。正式86页论文及旧科学数据未改；当前结果用于下一轮理论整合。
