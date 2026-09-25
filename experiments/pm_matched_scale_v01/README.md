# 匹配粒子与力网格尺度的受控塌缩

这一阶段检查上一轮三维塌缩中明显的空间离散误差。连续补偿球、参考背景、初始幅度与实际动量保持不变，仅执行两组每边粒子数与力网格数相同的配置：64/64和128/128。两组均沿用上一轮三项步长限制同时减半的时间控制。

本轮不增加clock分支，不重新选择初始过密度，也不引入平滑长度或力幅修正。它解决的是数值表示是否接近已知连续球形解的问题；结果不属于新观测拟合或星系形成模拟。

**两例均已完成。** 64/64的事件误差、最大半径误差及最大同调散布分别为2.1907%、20.2627%、9.2393%；128/128分别为0.3166%、3.1266%、5.3980%。128例径向精度明显改善，但半径误差仍未低于预设3%，两档事件差−1.8340%也未过1%要求。三个逐例失败和一个跨例失败全部保留；终态横向／径向力比由6.49%增至19.86%，不能说所有三维指标同步改善。实现检查95项、独立连续参照40项及完整产物审计440项通过，各自范围与科学门限分开报告。

- [运行前固定的协议](protocol.md)及[协议SHA256](protocol.sha256)
- [上一轮完整六配置](../halo_cooling_v01/reports/results_cn.md)及[本轮所据的设计备忘](../halo_cooling_v01/reports/pm_next_stage_design.md)
- [独立审计规范](reports/independent_review.md)
- [连续密度、球壳与初始力诊断](reports/shell_reference_cn.md)
- [全部结果与未通过门限](reports/results_cn.md)
- [实际执行及文件身份](results/pm_execution_record.json)
- [完整末态的独立审计](reports/independent_artifact_audit.md)
- [论文第4章](../../manuscript_cn/chapters/04_mathematical_framework.md)

## 如何判断结果

固定标签粒子的中位缩放定义半径代理，代理密度首次达到200时停止。除了事件时刻，还要比较共同事件前的完整采样轨迹、固定标签的径向散布、实际包围质量、二阶形状和四阶角向量。事件误差2%、半径误差3%、同调散布10%及双档事件差1%是运行前声明的实验门限。它们不是一般天体模拟的通用精度保证。

在理想代理密度200处，64/128力网格对应的核心半径仅约2.74/5.47个格胞。因此即使事件改善，也不能直接证明高密度对象内部已被充分分辨。本轮没有独立平移、进一步时间减半或clock分辨率矩阵；百分级绝对精度也不能单独证明约0.015%的clock差值可分辨。

独立连续参照包括初始局部密度的质量积分和33个固定拉格朗日球壳的物理半径方程。它们不调用PM力来校准引力幅度，亦不以球对称计算替代三维稳定性检验。所有旧目录保持只读，新增代码、原始结果、失败项及完整末态都进入本目录。

## 存档范围

两例均保存双精度末态位置与动量，分别压缩为独立文件，便于独立重算并控制单个Git文件大小。轨迹JSON保存全部接受步；NPZ保存径向剖面及投影，新投影轴使用CIC节点坐标`j/128`。输入依赖仍在同一仓库中，不需要相邻论文项目或下载观测数据。

代码核验、保存文件审计和科学精度验收分别报告；前两者通过不能抵消后者失败。已有正式结果不得覆盖，也不能删除不符合预期的配置。图件只读取存档，绘图不会重新积分。

## 复现

以下命令从仓库根目录执行，依赖Python、NumPy、SciPy，绘图另需Matplotlib。重算应使用独立检出目录，将本阶段已有`results/`和派生审计另存备份后再运行；保留两个旧实验的输入和原始结果。不要将新产生的摘要、版本或时间戳与本次粒子产物混用。生产程序拒绝覆盖已有正式产物，也不会自动重试。

```bash
python -B experiments/pm_matched_scale_v01/code/validate_matched_runner.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python -B experiments/pm_matched_scale_v01/code/shell_reference.py
python -B experiments/pm_matched_scale_v01/code/run_matched_pm.py --parallel 2 --fft-workers 4
python -B experiments/pm_matched_scale_v01/code/audit_matched_artifacts.py
python -B experiments/pm_matched_scale_v01/code/write_results_report.py
python manuscript_cn/figures/make_chapter04_matched_pm.py
```

科学门限未通过时，生产程序会保存全部结果并以状态1退出；这与执行崩溃不同。四小时资源界限在接受步之间检查，随后保存终态，落盘可略超出该界限。各案例开始和结束、命令、环境、执行身份以及限制原因均进入归档。初态连续力诊断每例调用一次PM力；正式轨迹正常计算各自起点力，末态独立审计再各调用一次旧全FFT，三种计算用途分别计数。

只查看或重画已有数据，无需重跑前三步。后续[气体接口设计](reports/gas_interface_next_design.md)已列出质量分账、膨胀项、热化和冷却的检验顺序，尚未执行；不能将这份设计当作本阶段已产生气体或星系的证据。
