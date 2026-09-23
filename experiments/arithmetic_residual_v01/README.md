# 算术残差的动力学可辨识性

2026-09-22 · 已完成 · 1,593例主计算、25例时间加密、独立积分与空间敏感性核对

**本轮已实际接入真实数学零点的间距残差，但没有检出它在预设主指标上相对于替代输入的额外特殊性。** 八个真实块全部处于两类替代集合的中心95%范围；汇总经验双侧秩为0.40和0.36。这些秩不是经校准的$p$值，未检出也不证明等效或排除其他物理桥接。

核心逆对数平方背景继续保留，真实物理相互作用仍未从黎曼结构导出。当前额外桥接是残差如何进入恢复力、索引如何映射到内部时钟；本轮没有宇宙学或原子能级观测拟合。

- [v0.5短稿](../../draft_v05/arithmetic_response_v05.md)：建议先读。
- [完整结果](reports/results_cn.md)：所有块、主汇总、质量限制与误差。
- [冻结设计](protocol.md)及[冻结记录](results/protocol_freeze.json)：在主动力学结果之前保存，未按结果修改；不是公开预注册。
- [来源核查](reports/source_inventory.md)、[方法学复核](reports/methodology_review.md)、[实施记录](reports/implementation_notes.md)。
- [主比较图](figures/arithmetic_primary_comparison.png)、[响应与控制图](figures/arithmetic_response_and_controls.png)、[图注](figures/captions.md)。同目录另有PDF。

## 使用的输入与控制

前作固定公开提交的10,000个零点缓存，经Odlyzko官方表和本地高位参考交叉核对。1024个开发间距仅用于归一化，另1024个评价间距分成八块；每块真实输入、99个IAAFT和99个配对精确节点谱输入，再加一例纯包络。

**输入匹配并非全部通过。** IAAFT有4/792条未达到预定5%节点幅度谱门限，全部保留；窗化插值后的实际驱动谱也未严格匹配。因此未隔离纯高阶算术结构的独立作用。精确节点谱族是配对敏感性，不是独立重复实验。

## 复现

需要Python、NumPy、SciPy、Matplotlib。从工作区根目录运行：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python riemann_model/experiments/arithmetic_residual_v01/code/prepare_inputs.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python riemann_model/experiments/arithmetic_residual_v01/code/run_lattice.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python riemann_model/experiments/arithmetic_residual_v01/code/run_reference.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python riemann_model/experiments/arithmetic_residual_v01/code/analyze_and_plot.py
```

第一步只有原始缓存缺失时下载已列明来源，另只读现有Clock的83点参考。主要输入和随机种子固定，文件SHA-256记录在JSON中。主积分约51秒（本机观测值，包含所选细步长计算）；独立高精度参考需更长。时间不是跨硬件性能承诺。

第三步使用收紧后的最终设置，历史`reference_initial_*`保留为首轮精度记录。没有必要为复现最终结论重新运行已知较松设置。重跑会覆盖对应生成结果，请将修改后的探索另存目录。

## 结果入口

| 文件 | 含义 |
| --- | --- |
| [input_summary.json](results/input_summary.json) | 来源、展开、固定分段与种子、质量失败清单。 |
| [input_quality.csv](results/input_quality.csv) | 每条输入的节点及实际驱动匹配诊断。 |
| [inputs.npz](data/inputs.npz) | 1593条节点序列、样条系数与案例标识；`case_block`采用0起始编号。 |
| [lattice_summary.json](results/lattice_summary.json) | 主计算、17项数值核验、步长精度与旧基线回归。 |
| [reference_summary.json](results/reference_summary.json) | 独立DOP853、22例参考及tight抽检、保留首轮精度失败。 |
| [comparison_summary.json](results/comparison_summary.json) | 冻结主比较、24项数值比较检查、配对及质量敏感性、哈希。 |
| [comparison_blocks.csv](results/comparison_blocks.csv) | 八块×两族完整主指标与经验秩。 |
| [lattice_cases.csv](results/lattice_cases.csv) | 全1593例主指标及次要描述。 |
| [artifact_checks.json](results/artifact_checks.json) | 当前产物、公式、链接、来源哈希及原稿保全核验。 |

NPZ轨迹按`(case,time)`存标量，所选完整场按`(selected_case,time,grid)`保存；输出为201个时刻。`velocity`是$\dot q$，正则动量为$a\dot q$。源码中的数组索引与图表中从1起始的块号不可混用。

数值检查通过只说明指定实现和所考察精度可用。当前输入质量失败、独立参考首轮失败、旧v0.3的额外精度提示均保留，未汇总成全项目所有检查通过。
