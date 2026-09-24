# 第8章图注与来源

两图只重绘已冻结的结果表，未增加拟合、随机试验、动力学运行或物理观测。复现：

```bash
python manuscript_cn/figures/make_chapter08_figures.py
```

## 图8.1：合成恢复中的幅度检出与形状识别

[PNG](chapter08_recovery_cn.png) · [PDF](chapter08_recovery_cn.pdf) · [SVG](chapter08_recovery_cn.svg)

**图8.1　既定观测设计下的合成恢复。** 使用 King 目录的293条记录、今日原子钟漂移误差，以及 Keck、VLT 各一个自由仪器偏移；背景、参考尺度和噪声模型固定。**(a)** 预先固定指数 $s=2$，用双侧5%检验拒绝零漂移。横轴为正向注入的 $D_0/\sigma_{\rm clock}=0,1,2,5$，其中 $\sigma_{\rm clock}=2.5\times10^{-19}\,\mathrm{yr}^{-1}$。四点的拒绝率依次为4.795%、16.960%、51.760%和99.885%。每点使用20,000次评价重复，误差棒为约95%的蒙特卡洛误差；叉号为结果表中同一注入点的解析功效，虚线为名义5%拒绝率。零注入点表示误报检查，不能称为真实信号的检出功效。**(b)** 分别注入 $s=1,2,3$，固定 $D_0=1.25\times10^{-18}\,\mathrm{yr}^{-1}$，比较常数与三个时间响应候选的 AIC 选择频率。三行使用相同噪声实现，四列频率均为0.015%、50.135%、0.020%、49.830%；保留常数候选，未对三个变化候选重新归一化。近共线的时间响应无法由幅度检出自动区分，中间指数很少成为最优者不构成其物理劣势。选择频率不是理论成立概率，AIC 胜出也不是幅度检出。五倍钟误差的合成注入不表示与今日实测均值相容。

来源：`experiments/alpha_recovery_v01/results/recovery.csv` 的 `king_clock_offsets` 行；`summary.json` 用于核对设计与重复次数。`shape_geometry.csv` 额外校核五倍钟误差注入下 $s=2$ 对相邻指数的最小无噪声剩余 $\Delta\chi^2=3.4086648410\times10^{-9}$，不作为另一项检出统计。

## 图8.2：算术输入与条件替代族的终点比较

[PNG](chapter08_arithmetic_cn.png) · [PDF](chapter08_arithmetic_cn.pdf) · [SVG](chapter08_arithmetic_cn.svg)

**图8.2　固定终点的能量交换及其条件替代范围。** 主指标 $T=[E_R(20)-E_R(0)]/\mathcal H(0)$，纯包络基线 $T_0=0.04513942199602207$，横轴为 $10^4D=10^4(T-T_0)$。黑色菱形为八个真实零点残差块的读出，蓝、绿两族各含99条替代输入；细条为经验2.5%—97.5%分位范围，粗条为四分位范围，圆点为中位数。**这些范围是条件替代分布的描述，不是参数置信区间。** 黄色底纹的汇总行先按固定替代序号跨八块等权平均，再对99个平均值求分位范围；真实均值为 $0.04504179099443443$。所有真实块及汇总均值均处于两个对应的中央95%范围内。精确节点谱族由 IAAFT 的同一相位投影得到，是配对检查；不能将两族算作两次独立重复。四个节点质量失败保留在原定主比较中。插值、窗化之后的连续驱动及其导数功率谱未严格匹配，因而未建立黎曼特异的动力学响应，也未证明与替代输入等价。图中模型时间、能量及输入均未标定为宇宙年龄或物理观测量。

来源：`experiments/arithmetic_residual_v01/results/comparison_blocks.csv` 与 `comparison_summary.json`。图形仅作减去同一基线、乘 $10^4$ 的单位变换；没有改动原定主指标、块选择、终点或替代集合。经验尾部秩不画作校准 $p$ 值。

## 输入字节清单

下列路径以仓库根目录为基准。重绘脚本首先验证 SHA-256，保证没有隐式更换数据或方案。

| 输入路径 | SHA-256 |
|---|---|
| `experiments/alpha_recovery_v01/results/recovery.csv` | `70c5749a40e1c3e03563ff6d5e8c9a137bdef45d53154f355d0cf5eda1f057a7` |
| `experiments/alpha_recovery_v01/results/summary.json` | `c3806224d57c034218093da24e80c8c76a0f39791078ec61eaba3b08a2e25a6a` |
| `experiments/alpha_recovery_v01/results/shape_geometry.csv` | `e15b14f65cb721a202bef546ce03949d8b001f5989ead8cb975064363209e3bf` |
| `experiments/alpha_recovery_v01/protocol.md` | `3e0b6ecff2c73a9e830a4b7b481d07753e5953b86572c8c02ac56c73a1d1e100` |
| `experiments/arithmetic_residual_v01/results/comparison_blocks.csv` | `b1234f062c46209efb0acbf452c3ec65036b66b5d0b24994abbed26c8dc8e7e8` |
| `experiments/arithmetic_residual_v01/results/comparison_summary.json` | `278a9ae61df1d4fba46b4db4f240c5adc34e139213a5bd976b15904942d2f610` |
| `experiments/arithmetic_residual_v01/protocol.md` | `952804fa8aaf7ffa2092449010dc8a8c5d045ed319df1a7812a993ee565e9d56` |
