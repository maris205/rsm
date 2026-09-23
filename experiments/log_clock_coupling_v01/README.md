# 逆对数平方响应与自主时钟：完整反馈实验

2026-09-22 · 已执行 · 无量纲候选模型，没有观测拟合

本轮把用户要求保留的$1/\ln^2$响应实际接入非线性格点，同时计算格点对时钟的反作用。指定正能量分支下，两者交换能量；弱反馈近似预设的对数时间曲线，强一些的反馈产生可计算的偏离。有限格点系统还满足一个条件晚时定理：响应相对于其常数极限的尾部保留逆对数平方尺度。

**黎曼物理来源仍未完成。** 平方响应及时间桥接是输入假设；本轮没有输入真实零点或算术残差，也未从这些结构导出真实粒子相互作用。新的内部$R$时钟不同于v0.2依赖膨胀背景的正则标量候选。

- [v0.4短稿](../../draft_v04/log_clock_bridge_v04.md)：核心想法、方程、条件结论。
- [结果报告](reports/results_cn.md)：全部对照、误差与局限。
- [设计记录](protocol.md)：参数、控制组、比较规则；不是公开预注册。
- [解析推导](reports/derivation_and_limits.md)：正则结构、能流、宏观方程、反馈及晚时界。
- [前作与算术桥接状态](reports/arithmetic_bridge_status.md)：具体引用了什么，哪些联系尚未建立。
- [响应与能量图](figures/log_clock_response.png)、[场演化与收敛图](figures/log_clock_field_bridge.png)、[图注](figures/captions.md)；同目录另有PDF。

独立参考下，$I=10,100,1000$时，完整场与外部规定对数驱动的全时空相对RMS差异分别为7.410%、0.8634%、0.08787%。这是所选参数下的模型轨迹差异，未经过拟合，不是观测准确率。改变$I$也改变初始时钟能量。

## 复现

从工作区根目录运行，需要Python、NumPy、SciPy及Matplotlib。脚本只写本实验结果与图；第一步读取上一轮固定系数结果作可选的$b=0$回归核验，第二步没有导入主实现。

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python riemann_model/experiments/log_clock_coupling_v01/code/run_coupled_lattice.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python riemann_model/experiments/log_clock_coupling_v01/code/run_reference.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python riemann_model/experiments/log_clock_coupling_v01/code/compare_and_plot.py
```

主积分8例，使用同时kick–drift–kick；独立DOP853参考12例，包括耦合格点、外部驱动、冻结系数、零效应控制和连续Fourier–Galerkin场。相同格点的独立求解用于时间误差，独立连续场用于空间误差。主结果不是通过重新对齐相位、改变振幅或拟合曲线得到的。

## 结果文件

| 文件 | 内容 |
| --- | --- |
| [lattice_summary.json](results/lattice_summary.json) | 8例主格点元数据、实际能量误差、88项初值／有效域／零控制检查。 |
| [reference_summary.json](results/reference_summary.json) | 12例独立参考、容差与谱截止核对、27项参考检查。 |
| [comparison_summary.json](results/comparison_summary.json) | 反馈和控制组、二阶空间／时间误差、29项比较诊断、输入SHA-256。 |
| [comparison_feedback.csv](results/comparison_feedback.csv) | 三种$I$的反馈、系数、能量与独立实现差异。 |
| [comparison_spatial.csv](results/comparison_spatial.csv)、[comparison_time.csv](results/comparison_time.csv) | 网格和时间步收敛表。 |
| [lattice_modes.csv](results/lattice_modes.csv)、[reference_modes.csv](results/reference_modes.csv) | 各例模式、时钟、能量和功的采样记录。 |
| `lattice_*.npz`、`reference_*.npz` | 完整轨迹；变量和傅里叶归一化见JSON。 |
| [artifact_checks.json](results/artifact_checks.json) | 本轮文件、公式、链接及原稿保全核验。 |

注意：主NPZ用`M2/Hfield/Hclock/Htotal`，参考NPZ用`mass_squared/field_energy/clock_energy/total_energy`。两者保存的`velocity`均为$\dot q$，正则格点动量为$p=a\dot q$。连续参考的波数从$-K$排列到$K$，格点傅里叶系数采用FFT顺序。外部与冻结对照没有可称为守恒总能量的“场＋自由时钟”列。

本轮的检查均通过；它们检验实现及误差预算，不证明物理真实性。v0.3额外最细参考审计的一项未过门限仍保留在原报告中，不随本轮结果改写。
