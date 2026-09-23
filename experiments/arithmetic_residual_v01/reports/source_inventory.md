# 算术输入来源与前作代码核查

核查日期：2026-09-22。范围为本工作区的 `riemann_engine`、`dsc_theory`、`Cosmic-Chaos-Alpha`、`riemann_clock`、`riemann_hubble` 和 `riemann_model` 既有来源审计，以及零点前作公开仓库和 Odlyzko 官方数据页。没有改动旧项目，没有重跑前作的大规模映射拟合，也没有搜索工作区以外的项目。

本轮选定的主输入是**前作公开仓库中前 10,000 个数学零点的 float64 缓存**，而不是原拟合模型生成的特征相位、人工对数包络或实验室的带误差零点估计。缓存可以用独立公布的数学表交叉核对；它本身不构成物理观测。

## 1. 输入库存与精度的实际来源

| 输入 | 数量与索引 | 可核实来源和精度 | 本轮用途 |
|---|---|---|---|
| 前作 `python/riemann_10k_true.npy` | 10,000 个，数组形状 `(10000,)`，按零点序号 1–10,000 严格递增 | 固定公开提交；80,128 字节；dtype 为 float64。生成笔记本调用 `float(mpmath.zetazero(n).imag)`，没有显式设定 `mp.dps`，因此不宣称高于存储精度或区间认证。 | 主算术输入 |
| 本地 `riemann_clock/data/processed/mathematical_reference_zeros.csv` | 83 个，索引 1–81、4200、4201 | 表记录 mpmath 1.3.0、35 位十进制；生成代码 `analyze_zero_precision.py` 设 `mp.mp.dps=35`，校验代码 `validate_analysis.py` 设 50 位。README 明确它不是区间认证结果。 | 已完成 83 点数值交叉检查 |
| Odlyzko 官方 `zeros1` | 前 100,000 个零点 | 官方目录标明误差在 $3\times10^{-9}$ 以内；不能把它称为 35 位表。 | 独立公开来源，供本轮输入流程逐项交叉核对前 10,000 个 |

前作缓存固定在提交 `e3419dda3d5515afd91a5dd6e8d8b8398f359dde`：

- [缓存文件](https://github.com/maris205/riemann_logistic/blob/e3419dda3d5515afd91a5dd6e8d8b8398f359dde/python/riemann_10k_true.npy)；[直接读取地址](https://raw.githubusercontent.com/maris205/riemann_logistic/e3419dda3d5515afd91a5dd6e8d8b8398f359dde/python/riemann_10k_true.npy)。
- [生成笔记本](https://github.com/maris205/riemann_logistic/blob/e3419dda3d5515afd91a5dd6e8d8b8398f359dde/backup/mvp2/p3_riemann_10k_true.ipynb)，第一个代码单元循环 `range(1, 10001)` 后保存数组。
- `backup/mvp2/riemann_10k_true.npy` 与 `python/riemann_10k_true.npy` 在本次读取中逐文件 SHA256 相同。
- 两份缓存的 SHA256 均为 `6bb18ef724b76280ee91e56c9ad0669a8026c996bd969c710b3a708f28e57545`。

本次内存读取确认首项约为 14.13472514，末项约为 9877.78265401，全体相邻差正。与上述本地 83 个参考值转 float64 后比较，最大绝对差为 **0.0**。这验证了这些索引上的数值一致性；同用 mpmath 的参考不能代替独立算法或完整区间认证。本报告没有重算全部 10,000 个零点。

[Odlyzko 官方目录](https://www-users.cse.umn.edu/~odlyzko/zeta_tables/)给出前 100,000 个零点表及精度声明；[该表](https://www-users.cse.umn.edu/~odlyzko/zeta_tables/zeros1)与 mpmath 缓存具有不同的发布来源。主分析采用哪一个数组、逐项误差及文件哈希，应由本轮实际输入清单保存。本报告只核实来源声明，不提前声称随后交叉核验的结果。

本地 [Clock 数据说明](../../../../riemann_clock/data/README.md)、[35 位参考表](../../../../riemann_clock/data/processed/mathematical_reference_zeros.csv)、[生成代码](../../../../riemann_clock/code/analyze_zero_precision.py)与[更高精度校验代码](../../../../riemann_clock/code/validate_analysis.py)提供可复查的出处。其余检查目录中未找到更长的直接数学零点表；`Cosmic-Chaos-Alpha` 的零点热力图属于可观测精度设想，不能当作真实算术输入。

## 2. 前作真实输入与非自治驱动是两个对象

[出版稿](../../../../riemann_hubble/docs/mca-31-00193.pdf)的 Data Availability Statement 指向 `maris205/riemann_logistic`。论文明确把数学零点作为拟合目标。公开代码中至少有以下两类实现：

- [task14_rmt_unfolded.py](https://github.com/maris205/riemann_logistic/blob/e3419dda3d5515afd91a5dd6e8d8b8398f359dde/task14_rmt_unfolded.py) 的第 22–24 行设 `mpmath.mp.dps=15`，计算前 100 个零点并转 float64；第 33–40 行用平滑计数公式展开。它为本轮“数学零点—去趋势间距”提供直接的前作方法接口。
- [task12_13_experiments_v2.py](https://github.com/maris205/riemann_logistic/blob/e3419dda3d5515afd91a5dd6e8d8b8398f359dde/task12_13_experiments_v2.py) 将真实零点与转移矩阵特征相位分开；映射中的 `i` 是轨道迭代步，驱动为 `u_temp + k_opt/log(i+offset)**2`。热身步与总步数也属于轨道采样设置，并不是零点序号，更不是宇宙时间。
- [较早的 1000 点脚本](https://github.com/maris205/riemann_logistic/blob/e3419dda3d5515afd91a5dd6e8d8b8398f359dde/python/macro_1000_scale_find_1d_v1.py) 显式设 `mpmath.mp.dps=25` 后计算前 1000 点。不同脚本的精度设置不能移植成 10,000 点缓存的生成保证。

出版稿主式中的 $n$ 是离散映射步数／尺度参数化，谱比较另有零点序号 $j$。新文稿宜坚持使用不同符号：轨道步 $n$、数学零点序号 $j$、内部时钟 $R$、演化时间 $\tau$。把任意两者相连都需要另写桥接假设。

原非自治系数的 $k_1/\ln^2 n+k_2/\ln^3 n$ 是指定的平滑时间表，其中系数通过零点目标校准。仅把这个函数采样、减去另一个平滑函数，**不能声称得到了真实算术残差**。前作未直接输出一个已证明决定物理耦合的零点残差；本轮从真实数学零点生成残差并接入格点，是新的有效模型定义。相关边界已见[上一轮算术桥接说明](../../log_clock_coupling_v01/reports/arithmetic_bridge_status.md)。

## 3. 最小可复现的去趋势方法

固定使用正虚部 $\gamma_j$，定义无待拟合参数的平滑计数函数

$$
\overline N(\gamma)=\frac{\gamma}{2\pi}\ln\frac{\gamma}{2\pi}
-\frac{\gamma}{2\pi}+\frac78,
$$

及展开后的间距残差

$$
s_j=\overline N(\gamma_{j+1})-\overline N(\gamma_j),\qquad
r_j=s_j-1.
$$

这与前作 `task14_rmt_unfolded.py` 的展开步骤一致。只有相邻零点对可生成一个间距；10,000 个零点最多给出 9,999 个这样的残差，不能称 10,000 个独立间距。去趋势的目的是剥离已知平均密度；剩余起伏仍是数学序列的结构，不是宇宙时间中的实测涨落。

如需中心化或方差归一化，固定以训练段估计其参数，再原样用于验证／留出段；不能各段独立缩放后声称是未见数据预测。抽取时间驱动还需要预先固定序号到时钟的映射、平滑插值和幅度。幅度、滤波长度、时钟缩放不能根据留出结果反复调优。

对照至少区分平滑包络、真实算术残差与统计匹配替代残差。若替代序列利用整段真实数据的幅度谱或排序分布，它属于条件于整段统计量的替代检验；不能把它描述为严格只见训练数据的未来预测。相同势函数可以把任意输入传到输出，因此“真实输入与零输入产生不同轨迹”只能证明模型响应，不能单独证明黎曼结构优于一般相关输入，更不能证明真实物理相互作用。

## 4. 当前可支持的结论

已有一个小体积、有固定版本和生成代码的真实数学零点输入，可以立即开始本轮算术残差对照；不必重跑原论文 $10^{10}$ 步映射拟合。缓存的有限精度、展开定义、索引—时钟桥接和对照构造需与结果一同公开。此步骤推进的是“算术结构在给定候选动力学中是否有可区分的作用”；**尚未从黎曼结构导出真实物理相互作用**，也没有加入新的天文观测证据。
