# 受控塌缩与原初气体冷却接口

本阶段把 `cosmic_bridge_v01` 的正则时钟背景和标准引力连接推进到一个质量、初态及终点明确的球形塌缩基准，再计算指定气体状态的瞬时冷却时间。结果用于检验实现和条件物理连接，不是观测拟合，也不是完整星系形成。

六组PM已全部完成，完整保存产物的独立审计876项通过；9项逐例科学精度失败和2项跨例失败全部保留。增加粒子数使事件误差由6.49%降到2.95%，但半径误差仍为26.5%，尚未通过连续精度检验。连续球示例的时钟事件延后1.0264 Myr；同一质量、红移的标准冷却微物理保持相同。

- [运行前固定的协议](protocol.md)及[SHA256](protocol.sha256)
- [数值结果、分辨率与适用边界](reports/results_cn.md)
- [独立物理复核](reports/physics_independent_review.md)
- [后续PM方法检验设计](reports/pm_next_stage_design.md)（仅设计，未执行）
- [论文第4章](../../manuscript_cn/chapters/04_mathematical_framework.md)

球形初幅只由参考方程选择一次，使首次密度比200的事件发生在尺度因子0.5；所有PM配置和clock沿用同一初态及实际动量。事件200不等于维里平衡。25个质量、7个红移的冷却网格采用历史H/He碰撞电离平衡速率，其特征温度和均匀气体密度是指定接口，尚未由气体动力学求出。

## 复现与材料

命令从仓库根目录运行，依赖Python、NumPy、SciPy；绘图另需Matplotlib。已存数据可直接读，绘图不需要重新积分。下面是完整重算顺序：宜在独立仓库副本中先将本实验的`results/`移作备份，再执行，以便保留当前证据；须保留相邻`cosmic_bridge_v01/results/`背景输入。新运行会生成新的时间戳和产物身份，不应与旧PM结果混用。所有背景输入均在相邻实验目录内，旧代码保持只读，不需要访问其他项目或下载CAMB。

```bash
python experiments/halo_cooling_v01/code/spherical_collapse.py
python experiments/halo_cooling_v01/code/verify_sphere_profile.py
python experiments/halo_cooling_v01/code/run_cooling.py --events experiments/halo_cooling_v01/results/sphere_summary.json
python experiments/halo_cooling_v01/code/validate_refined_pm.py
python experiments/halo_cooling_v01/code/run_spherical_pm.py
python experiments/halo_cooling_v01/code/validate_physics_independent.py
python experiments/halo_cooling_v01/code/validate_pm_artifacts.py
python experiments/halo_cooling_v01/code/write_results_report.py
python manuscript_cn/figures/make_chapter04_halo_cooling.py
```

PM包含六个固定配置，其中细力网格为256³，计算开销显著大于独立球形方程。每个作业最多4个FFT线程，最多同时运行两个。保存全部轨迹、剖面、投影和诊断；两个主例另存完整双精度末态粒子坐标与动量。其余配置的独立存档核对限于保存统计和运行记录，不声称可从它们重建完整末态。

`sources/`归档的Grackle官方平衡工具仅作为速率验证参照，不参与生产冷却值的生成。归档保留原始文件、[完整许可](sources/GRACKLE_LICENSE)及[来源和哈希](sources/grackle_source_metadata.json)，固定提交`af7939494ce65007887ada7b98d1813df6843346`。`rates='cen'`是本实验明确选择的历史解析速率集，既非完整Grackle化学求解器，也非当前默认微物理。

## 结果解释

实现验证通过与物理近似适用是两件事。局部电离平衡弛豫慢于冷却的格点作为近似警示保留；不能为了得到星系形成的正结果而移除这些点或给冷却率增加任意时钟因子。在相同质量和红移下，两背景的标准冷却微物理必须给出相同结果。只有背景年龄和形成历史不同，才能在不同事件处取到不同的条件读数。

新接口沿用逆对数响应的研究动机，但平方指数、真实原子响应和引力的黎曼来源依然没有完成。标准引力、原子速率及补偿初态是本轮的物理输入；气体守恒方程、非平衡化学、热化、恒星形成和反馈属于后续工作。

已有PM结果不会被默认覆盖；仅在代码、协议、背景及球形校准文件均未改变时，才可用`--resume`复用已完成案例；重新生成球形校准后应作完整新运行。科学精度门限未通过时，程序保存完整结果并以非零状态退出，这不等于进程崩溃或结果丢失。存档投影数组的CIC节点为`j/N`，`projection_axis`是图像像素中心坐标，二者相差半个分析格胞，不能在物理相位分析中混用。图4.5不使用这些投影轴。
