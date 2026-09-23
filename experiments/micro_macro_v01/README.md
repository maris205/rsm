# 局部非线性格点 → 连续场：首轮基准

**状态：已完成，保留一项最细参考精度提示。** 本轮将[研究主线](../../RESEARCH_DIRECTION.md)中的第一条微观—宏观连接落实为方程、独立数值实现和收敛检查。

核心结果是在指定平滑初态与有限时间内，基础空间误差约为二阶，加入首个格距修正后约为四阶；Verlet时间误差独立呈二阶，非线性第三空间谐波也收敛。它是已选择Hamilton量的数学基准，没有观测拟合或黎曼微观定律的确认。

- [结果报告](reports/results_cn.md)：数值表、主要结论、频谱适用范围及未通过的精度诊断。
- [数学推导](reports/derivation_and_limits.md)：正则归一化、连续极限、空间修正、离散相位和能量。
- [实验设定](protocol.md)：固定参数、初值、误差定义及对照方法。
- [简短论文段落 v0.3](../../draft_v03/micro_macro_bridge_v03.md)。

![微观格点与宏观连续场](figures/micro_macro_evolution.png)

## 复现

需要Python 3、NumPy、SciPy、Matplotlib。实际版本及代码／输入哈希记录在JSON中。在工作区 `dsc_world` 顺序执行：

```bash
python riemann_model/experiments/micro_macro_v01/code/run_lattice.py
python riemann_model/experiments/micro_macro_v01/code/run_reference.py
python riemann_model/experiments/micro_macro_v01/code/compare_and_plot.py
python riemann_model/experiments/micro_macro_v01/code/check_finest_reference.py
```

前两个脚本各自独立实现方程。所有命令会覆盖本实验中的对应生成文件；最后的精度审计本轮因一项门限未满足而以非零状态退出，并已写出全部诊断。复现时应阅读JSON中的逐项状态，不能将这一已记录的限制隐藏或改成通过。本轮主格点11例约11秒，参考11例的纯积分累计约4.5秒；这些是本环境下的小模型记录，不是普适的速度比较。

## 文件

| 文件 | 内容 |
| --- | --- |
| [lattice_summary.json](results/lattice_summary.json) | 11个Verlet案例、86项解析／实现检查。 |
| [reference_summary.json](results/reference_summary.json) | Fourier归一化、11个独立参考、22项核查及轨迹哈希。 |
| [comparison_summary.json](results/comparison_summary.json) | 空间二阶／四阶、时间二阶、谐波、低频占有与20项比较。 |
| [reference_precision_audit.json](results/reference_precision_audit.json) | 最细格点收紧容差：4项通过，1项1%参考门限未通过。 |
| [空间误差CSV](results/comparison_spatial.csv)、[时间误差CSV](results/comparison_temporal.csv)、[色散CSV](results/comparison_dispersion.csv) | 便于再分析的主要数值。 |
| [格点模态CSV](results/lattice_modes.csv)、[参考模态CSV](results/reference_modes.csv) | 各例1001时刻的模式系数与能量。 |
| `results/lattice_*.npz`、`results/reference_*.npz` | 原始场、速度、能量或Fourier系数轨迹。 |
| [图及图注](figures/captions.md) | 两张科研图的PNG、PDF与准确说明。 |

精度提示的具体含义：最细参考的全场最大RMS变化占修正误差1.092%，略超1%诊断门限；主要全时空间误差仅变化0.00781%。因此保留已观察到的近四阶收敛，同时不宣称全部精度要求均通过。数值截止$K=24$也不等于所有保留模式都属长波，适用性依赖本次初值的实际低频占有。
