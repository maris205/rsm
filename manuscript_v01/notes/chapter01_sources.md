# 第一章来源与结果定位

2026-09-23。此文件是逐章写作的资料说明，不属于论文正文。本轮新增章节组织与第一章，未重新运行实验。外部材料按下列实际阅读范围使用；不是对全部文献或前作证明的独立审计。

## 外部来源

编号对应[文献表](../references.md)。

| 文献 | 已核对的位置 | 用于第一章的范围 |
| --- | --- | --- |
| [1] Berry–Keating 1999 | 作者全文第236—239页，特别是第238页关于自伴算符与Hilbert–Pólya思想的段落；SIAM正式题名、卷页与DOI。 | 谱、动力学与算术的联系及研究动机。没有把其候选算符研究作为本模型的推导。 |
| [2] Montgomery 1973 | 作者大学托管全文第181页定理的RH假设及范围；第183—184页猜想、Dyson所指出的随机矩阵成对关联。 | 区分条件定理与更广猜想，避免把全套GUE统计写成已证同一性。 |
| [3] Odlyzko 1987 | 原论文首页摘要（大学托管副本）、作者文献页及出版商Crossref卷期页码记录。 | 零点间距的数值统计吻合与谱高度依赖；不是对某种原子系统的逐项能级识别。 |
| [4] Wang 2026，素数 | 本项目[前作核查](../../reports/sources_v01_audit.md)第1节；零点出版稿第25页Appendix A；本轮再次核对出版商登记元数据。 | 归属式陈述作者研究及报告的结果。出版商全文页面受403限制，本轮未新增对全文证明的独立审计。 |
| [5] Wang 2026，零点 | 同一核查；本地出版版第1页摘要、第3—4页§2.2、第18、22—25页的限制与Appendix A；本轮再次核对出版商登记元数据。 | 非自治二次映射、逆对数驱动与经验谱对应，明确校准、留出及局部GUE边界。 |
| [6] Hairer–Lubich–Wanner 2003 | 作者原始摘要与Cambridge出版条目；Crossref正式期刊卷页及DOI。 | 标准几何数值方法背景；不把辛性、Verlet方法或既有结构理论列为本项目原创。 |

新增经典来源的直接入口：[Montgomery全文](https://websites.umich.edu/~hlm/paircor1.pdf)、[Odlyzko原论文副本](https://web.williams.edu/Mathematics/sjmiller/public_html/ntprob19/handouts/computational/Odlyzko_DistrSpacingsZerosZeta.pdf)、[Berry–Keating全文](https://michaelberryphysics.wordpress.com/wp-content/uploads/2013/06/berry307.pdf)、[Hairer等作者摘要](https://www.unige.ch/~hairer/preprints/gniverlet.html)。

## 本项目结果

| 第一章表述 | 现有依据 | 必须一同保留的条件 |
| --- | --- | --- |
| 内部时钟与格点共享Hamilton量，允许双向反馈 | [v0.4短稿](../../draft_v04/log_clock_bridge_v04.md)§1—2及[推导](../../experiments/log_clock_coupling_v01/reports/derivation_and_limits.md) | 候选耦合；均匀时钟不等于已建立的局域相对论场。 |
| 逆对数平方晚时主导行为 | 同上§3及完整推导 | 纯包络、固定有限格点、有限初始能量、正回复系数下界和正初速度等分支条件；平方指数预设，未由结论选出。 |
| 长波方程、修正及误差比较 | [v0.3结果](../../experiments/micro_macro_v01/reports/results_cn.md)、[v0.4结果](../../experiments/log_clock_coupling_v01/reports/results_cn.md) | 平滑初态、有限时段、分开的空间与时间误差；保留v0.3额外精度门限未通过。 |
| 八个真实块位于替代集合中心95%范围 | [v0.5结果](../../experiments/arithmetic_residual_v01/reports/results_cn.md)、[比较摘要](../../experiments/arithmetic_residual_v01/results/comparison_summary.json) | 预设主指标、固定耦合、两种替代族；经验秩不是校准$p$值，未完成等效性或检出能力分析。 |
| 实际驱动及其导数的匹配不足 | [方法学复核](../../experiments/arithmetic_residual_v01/reports/methodology_review.md) | 节点谱匹配不自动保留窗化插值后的谱；不得声称已经隔离高阶算术结构。 |

第一章以结果概述为限。详细参数、失败记录、误差与统计解释将在第5章及附录完整交代，现有原始记录保持可追溯。
