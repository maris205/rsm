# α时间律：模拟恢复与可辨识性实验 v0.1

**结果：在指定误差模型下能正确恢复注入的漂移幅度，但现有观测设计不能识别逆对数平方的指数。** 联合幅度信息几乎全部来自原子钟；即使人为注入强漂移，$s=1,2,3$仍近乎不可分。

这是采用真实红移与误差的合成响应实验，没有把模拟当成新增观测，也没有用实测α或钟漂移均值作注入中心。

- [中文结果报告](reports/recovery_results_cn.md)
- [预先记录的实验协议](protocol.md)
- [观测设计审计](reports/data_design_audit.md)
- [统计设计审查](reports/statistical_design_review.md)
- [统计结果独立复核](reports/statistical_result_review.md)
- [绝热与数值相位预算](reports/phase_budget_cn.md)
- [机器可读汇总](results/summary.json)、[全部133种情形](results/recovery.csv)、[形状距离](results/shape_geometry.csv)

## 图示

[幅度恢复与信息来源](figures/recovery_capability.png) · [指数混淆热图与理想精度预算](figures/shape_identifiability.png) · [完整图注与矢量PDF](figures/captions.md)

## 方法和范围

King目录293测量行、131个QSO，另加一次当前原子钟漂移标准误；ESPRESSO单吸收体作为独立分析分支。7种数据／干扰／协方差设计，每个设计用30,000次独立零模拟校准扫描阈值，再用20,000次评估重复。不同注入复用噪声，133种情形不是133份独立模拟证据。

比较常数与$s=1,2,3$；背景、参考尺度、输入误差与仪器偏移的条件写在协议中。用QR分解模拟精确线性高斯充分统计量，并与完整观测向量的最小二乘交叉核对。正负幅度、零信号、相关误差和故意使用错误协方差的结果全部保留。

## 复现

从工作区 `/root/autodl-tmp/dsc_world` 运行，需要Python、NumPy、SciPy、Matplotlib：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python riemann_model/experiments/alpha_recovery_v01/code/run_recovery.py
python riemann_model/experiments/alpha_recovery_v01/code/check_phase_budget.py
python riemann_model/experiments/alpha_recovery_v01/code/plot_recovery.py
```

脚本读取相邻的`riemann_constant/data/`，只写本实验目录；源文件哈希保存在汇总JSON。当前误差预算是有条件的设计计算，不是完整光谱或钟时间序列的重新归约。独立结果核验的复现方法见其报告。

[返回项目入口](../../README.md) · [最小公设](../../draft_v01/minimal_framework_v01.md)
