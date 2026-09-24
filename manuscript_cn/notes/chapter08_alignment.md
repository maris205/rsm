# 第8章与原稿的对应、来源及检查

## 原结构对应

对应原文件：[Chapter 8](../../chapters/09_chapter_08_simulation_case_studies.md)。保留3个主节、每节“设置—执行—结果”三个小节，共9个三级小节；为便于引用，原反复使用的1、2、3改为8.1.1—8.3.3。原章没有插图，本次补充两张实际结果图。

| 原文位置 | 保留的研究问题 | 当前处理 |
| --- | --- | --- |
| 8.1.1 α设置 | 逆对数时间响应如何与观测联系。 | 用固定背景、尺度、当前漂移归一化替代从第一个奇异tick起算；平方指数为主候选，1与3为对照。 |
| 8.1.2 α执行 | 数据怎样进入比较。 | 使用293条记录的设计、压缩钟误差与七种合成场景；分别报告零校准、评估重复、相关性和偏移。 |
| 8.1.3 α结果 | 能否恢复幅度和辨认形式。 | 用实际恢复率、覆盖率、形状距离及钟信息占比替代未经支持的零调参R²=0.974与漂移证明。 |
| 8.2.1 膨胀设置 | 内部慢变量能否影响宏观背景。 | 区分固定盒正则R、正则χ指数势背景与Hubble现象学；Planck和SH0ES均为今天H0推断。 |
| 8.2.2 膨胀执行 | 如何检验动力学和数值实现。 | 接入固定系数长波收敛、完整时钟反馈、独立DOP853／Radau背景求解及闭合条件。 |
| 8.2.3 膨胀结果 | 背景延拓带来什么可比较后果。 | 同今日漂移和同时间历史下约16.7%的传递差异；背景恢复与真实距离拟合分开，未声称张力消失。 |
| 8.3.1 实验室设置 | 怎样定义灵敏度及误差。 | 区分无量纲频率精度、年逆量纲漂移和数值相位误差；给出有效振子绝热参数，未建立物理噪声底。 |
| 8.3.2 精度执行 | 能否排除数值误差并识别输入机制。 | 用实际1,593例算术残差对照及收敛核验落实执行顺序；保留节点匹配、连续驱动及历史精度限制。 |
| 8.3.3 未来结果 | 什么结果可以拒绝指定候选。 | 给出钟时间序列的协方差预算及与高红移检验的分工；自由幅度含零时未检出不能排除全族，未承诺10年或10⁻²⁰必发现。 |

本章沿用假设的逆对数平方响应，没有把已有数值成果升级成真实物理相互作用的导出；数学零点、模型时钟与宇宙年龄继续分别定义。前作[7]、[8]提供动机及来源，不被当作物理同构证明。本文没有增加新观测、拟合、动力学扫描、计时或实验平台。

## 本章新增内容

正文新增6个编号式、3张表和2张矢量图；参考文献继续使用既有51条。

- 式(8.1)：继承已有α模型的当前漂移归一化，明确年与时间传递的量纲。
- 式(8.2)：定义已有形状几何指标，即固定真均值与替代候选重拟合后的协方差加权残差；不解释成已校准显著性。
- 式(8.3)：沿用自主背景候选的当前漂移归一化传递。
- 式(8.4)：冻结有效振子的辛Euler相位关系；状态一阶与频率偏差二阶分别说明。
- 式(8.5)：沿用算术主指标，以共同纯包络基线平移；原主汇总为平均T，平均D与其等价。
- 式(8.6)：未来钟序列的线性广义最小二乘预算，是设计表达，不冒充实际十年序列分析。

表8.1为反馈比较，表8.2为背景候选间的同归一化传递，表8.3分开数值核验和节点输入匹配。没有增加共同物理参数的联合拟合；不同候选及项目不相加为独立支持。

## 来源与复现边界

本章实验源均在写作前的固定仓库提交[2657a17](https://github.com/maris205/rsm/tree/2657a17a34439dd2a9bf0747f608305c04fa533d)内已有。以下为当前仓库中的直接入口；相关科学文件未在本次写作中修改。

| 实验 | 协议与报告 | 主要结果 |
| --- | --- | --- |
| α合成恢复 | [协议](../../experiments/alpha_recovery_v01/protocol.md)、[结果](../../experiments/alpha_recovery_v01/reports/recovery_results_cn.md)、[独立统计核对](../../experiments/alpha_recovery_v01/reports/statistical_result_review.md) | [recovery.csv](../../experiments/alpha_recovery_v01/results/recovery.csv)、[summary.json](../../experiments/alpha_recovery_v01/results/summary.json)、[shape_geometry.csv](../../experiments/alpha_recovery_v01/results/shape_geometry.csv)。 |
| 数值相位预算 | [报告](../../experiments/alpha_recovery_v01/reports/phase_budget_cn.md) | [phase_budget.json](../../experiments/alpha_recovery_v01/results/phase_budget.json)。 |
| 固定系数长波基准 | [结果](../../experiments/micro_macro_v01/reports/results_cn.md) | [空间／时间汇总](../../experiments/micro_macro_v01/results/comparison_summary.json)、[保留的严格精度失败](../../experiments/micro_macro_v01/results/reference_precision_audit.json)。 |
| 完整时钟反馈 | [协议](../../experiments/log_clock_coupling_v01/protocol.md)、[结果](../../experiments/log_clock_coupling_v01/reports/results_cn.md) | [反馈表](../../experiments/log_clock_coupling_v01/results/comparison_feedback.csv)、[汇总](../../experiments/log_clock_coupling_v01/results/comparison_summary.json)。 |
| v0.2自主背景 | [协议](../../experiments/action_completion_v02/protocol.md)、[数值结果](../../experiments/action_completion_v02/reports/numerical_results_cn.md)、[独立核验](../../experiments/action_completion_v02/reports/independent_background_review.md) | [背景表](../../experiments/action_completion_v02/results/background_summary.csv)、[独立比较](../../experiments/action_completion_v02/results/independent_background_checks.json)。 |
| 算术残差控制 | [协议](../../experiments/arithmetic_residual_v01/protocol.md)、[结果](../../experiments/arithmetic_residual_v01/reports/results_cn.md) | [各块统计](../../experiments/arithmetic_residual_v01/results/comparison_blocks.csv)、[主汇总](../../experiments/arithmetic_residual_v01/results/comparison_summary.json)。 |
| Hubble接口 | [第5章结果来源快照](../../reports/chapter05_sources/hubble/results_source.md)、[来源说明](../../reports/chapter05_sources/README.md) | 复用第5章数值；完整校准代码、数据及48次恢复档案的公开复现缺口仍未补齐。 |

α设计中相同视线相关ρ=0.25只是敏感性假定。每设计30,000零校准加20,000评估；同注入和不同指数复用噪声，133条件不是独立证据。主99.9998159%钟信息对应“King+clock+两偏移”，没有与第5章无偏移观测结果的99.9976%混淆。正向注入的功效数字与重绘图一致；覆盖率为预先固定真指数下的区间，不是模型选择后区间。

α相位预算含s=1、2、3，D0=0、±2.5×10⁻¹⁹、±1.25×10⁻¹⁸年⁻¹，κ=±7，共30条件。4,201点网格的绝热指标是给定ν0=5×10¹⁴Hz下的局部尺度比，不证明完整原子动力学或数值长时稳定。恢复计算使用解析响应，不逐载波积分。

反馈表的空间／时间RMS和积分误差均以相同对照归一化；变化I时也改变初始钟能。固定系数近四阶空间结果不属于完整反馈求解器；完整反馈的已有空间比较仍为二阶。主时钟的正则R动能与v0.2正则χ动能在换变量后不等价。背景的ε是示例，今天H0用Λ闭合固定，未预测H0；独立Radau沿用主闭合根，此限制继续明确。

算术所有1,593例的初始总能量同为1.3218983914666935，初钟能均为0.5，窗在初点为零。图分位数从T减去共同T0，非不同初能归一化的产物。25例时间敏感性最大ΔD/IQR为3.9731×10⁻⁶（非百分数）；空间64→128最大值为0.00260637，即0.260637%。八块不重叠未证明独立，两族配对不是独立复现，中心范围不是参数置信区间，0.40／0.36不是已校准p值。保留四个节点质量失败及连续驱动谱匹配不足，后者没有预设合格门限，不写成预定门限失败。

## 图件与静态核查

[重绘脚本](../figures/make_chapter08_figures.py)只读取已存档的表格并校验7项输入／协议SHA-256，未增加科学计算。[独立图注](../figures/captions_chapter08.md)记录数据过滤、单位变换、经验分位定义及全部哈希。

图8.1左图误差棒为1.96倍固定指数Monte Carlo标准误，约95%；右图保留常数列的3×4 AIC频率，三行均为0.015%、50.135%、0.020%、49.830%，档案确实相同而非复制错误。图8.2黑菱形、细95%范围、粗IQR、圆点中位数与正文一致。两图均已目视检查，中文、公式和符号清晰，无字体替换字符。

本章科学检查由原生协作代理分别核对α／精度、背景／反馈、算术／对照；其性质为作者侧辅助核查，不是外部同行评审。核查后补充绝热指标及示意载频、正向注入、固定系数b=0、背景示例范围和共同初能，未改科学结果。

前轮[摘要及第0—5章审阅](../../reports/manuscript_cn_body_review_20260923.md)的整体修订继续待办。本章的实际案例汇总不代表该轮建议全部处理。下一步按原设计进入第9章结论，完成主框架的中文逐章版本，再集中处理整体组织与重点突破。


## 构建交付

本章完成时累计单栏PDF为67页：第8章在第57—63页，图8.1／8.2在第58／62页，表8.1／8.2／8.3在第59／60／62页，参考文献在第64—67页。89个公式标签、51个参考锚点、56个目录项、337个有效PDF链接和14幅矢量图通过检查；40项构建输入哈希与当前文件一致，无缺字或溢出。图注与图、表题与表体同页，新章和参考页已逐页目视检查。旧章起页及旧图落页保持不变，本轮未修改构建器。本轮六份导航／正文／说明Markdown的204个本地链接另行检查通过，三个主节与九个小节的编号和顺序均已核对。

交付PDF的SHA-256为`3cc8bfccfd63c4b777136c8d89220a820d5333e71b1e9a06867578b13b334caf`。后续累计稿随新增章节更新，此处页数和哈希记录本次交付。
