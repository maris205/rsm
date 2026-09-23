# 第1章：原文对应、观测更新与小论文来源

2026-09-23。新稿为[第1章：系统异常](../chapters/01_system_anomaly.md)，对应原[Chapter 1: System Anomaly](../../chapters/02_chapter_01_system_anomaly.md)。保留原章提要及1.1—1.5顺序，以“两朵乌云”引出物理动机，最终过渡到第2章的算术来源。

## 原章结构与修订

| 原部分 | 新稿处理 | 补入的内容 |
| --- | --- | --- |
| Abstract | 保留膨胀率与真空能两条主线；区分观测张力与理论自然性问题。 | 明确内部响应的假设地位与四个项目的检验作用。 |
| 1.1 Webb's “Death Sentence” | 写成观测分歧及测量链，删除未经证实的Kelvin引语和“判死刑”结论。 | 2025年JWST/HST造父变星、CCHP TRGB及2026年H0DN结果；说明不同路径及数据重叠。 |
| 1.2 Server Clock Drift | 保留时钟响应想法，定义$H(t)$、$H_0$及原子钟读数映射。 | 区分坐标重参数化与相对物理响应；加入DESI最新距离进展、校准拟合与条件标量恢复。 |
| 1.3 Floating Point Truncation Error | 保留微观—宏观分层思路，改成真空能与有效引力问题。 | 硬截止示意、重整化总量、稳定性；区分物理截止与求解精度；标明当前经典格点的适用范围。 |
| 1.4 Excessive Uptime | 保留工程诊断的提问方式，写成共同慢变量与观测接口。 | 明确逆对数平方响应位置；纳入α、原生光谱和CMB结果，汇总为表1.1，保留先常规解释再函数比较的顺序。 |
| 1.5 Conclusion | 保留转向“源代码”章的逻辑，但具体指向数学输入与动力学联系。 | 总结已有工作提供的约束，不把四个项目相加为统一机制的显著性。 |

本章新增式(1.1)—(1.5)与一张接口表。原第1章没有正文图片；本轮沿用前置图A/B，不用假想观测曲线替代数据图。专门的数据图、声学峰和数值结果留待第5、8章按对应材料展开。

## 实际核查与引用范围

新增[参考文献14—24](../references.md#r14)，原1—13编号保留。外部来源核查截至2026-09-23，为与本章有关的定向核查，不是穷尽式文献综述。

- JWST造父变星：核对作者预印本`2509.01667`摘要、样本及联合$H_0$定义；未写成独立于HST的新增高斯先验。
- CCHP：使用`2408.06153v3`及ApJ出版信息，70.39一项明确对应HST+JWST的TRGB路径。2025年勘误修改附录B图B1，不影响所用主结果。勘误出版商全文访问受限，其内容由[该期刊论文的公开全文镜像](https://www.researchgate.net/publication/397458696_Erratum_Status_Report_on_the_Chicago-Carnegie_Hubble_Program_CCHP_Measurement_of_the_Hubble_Constant_Using_the_Hubble_and_James_Webb_Space_Telescopes_2025_ApJ_985_203)核对，不将镜像作为另一篇独立研究。
- H0DN：核对作者预印本与A&A书目，区分2025年首稿与2026年期刊出版；联合网络值不当作独立SH0ES结果。
- DESI：核对DR2 II及2026年DR2 IV作者预印本。后者含Lyα全形状的Alcock–Paczyński信息，不泛称所有信息均为BAO；所列暗能量偏好明确依赖模型和数据组合。
- 真空能：核对Weinberg出版商全文相关部分及Padilla讲义关于辐射不稳定性的章节。Padilla的arXiv提交记录是2015年，HTML转换显示的2026年日期不作为新论文年份。式(1.3)是标准单自由度截止估计，式(1.4)是说明重整化分工的示意式，不是本项目新增的真空能计算。

## 四个小项目怎样进入本章

| 项目 | 本章使用的记录 | 核对重点 |
| --- | --- | --- |
| Hubble | [原报告快照与来源](../../reports/chapter01_sources/README.md) | 三族$\Delta\chi^2$、时钟与正性条件、校准来源、48次背景恢复。源报告尚未在原项目HEAD中提交，因此保存了字节一致的有限文本快照。 |
| Constant | [固定版本研究报告](https://github.com/maris205/riemann_fine_structure_constant/blob/3827e38f1d7e633cbd1956a8c34fa363b83b54aa/reports/research_report_zh.md) | 293条旧测量与原子钟联合结果；$\Gamma$依赖项目尺度约定；零相容与指数难辨认。 |
| Clock | [固定版本光谱诊断](https://github.com/maris205/riemann_clock/blob/5b2b6b77bd544a812fb9140a9586189448f0a40f/experiments/highz_feasibility_2026-09-20/reports/publication_followup_results_cn.md)及[实验设计](https://github.com/maris205/riemann_clock/blob/5b2b6b77bd544a812fb9140a9586189448f0a40f/experiments/prediction_test_2026-09-21/readme.md) | 17/32数据数量、同线约−59 m/s负对照、四个条件位移与两个失配诊断；均值、额外不确定度和能级分布不混同。 |
| CMB | [固定版本结果](https://github.com/maris205/riemann_cmb/blob/f42db51748003fc83be064dade9bb52de5eb8319/readme_cn.md) | 2048实现、训练/保留峰的相反变化、自由初始场退化；不将下载ACT/SPT产品写成已经完成拟合。 |

本轮没有新增观测、数值实验或重新优化模型。两篇2026年已发表算术前作继续以文献7、8引用；四个物理方向按各自工作稿状态列示。正文介绍其结论与接口，详细参数、数据和公式保留在原报告以及后续专章。

## 复核与累计阅读版

由未参与正文起草的代理对照原章及本地四项目报告，只读核对了公式、数字与结论范围，未发现需要修改的实质性科学问题。这属于同一模型系列的协作检查，不代表外部同行评审。观测来源与数值分析记录分别核对，没有为本章重新计算显著性。

[累计单栏PDF](../paper.pdf)共17页：摘要与两幅矢量图在第1—3页，目录在第4页，第0章在第5—10页，第1章在第11—15页，24条参考文献在第16—17页。12个公式标签、24个参考锚点和112个PDF链接通过构建检查，未见缺字或版心溢出；第1章五页的公式和表1.1已目视检查。另核对14个相关Markdown文件的本地链接、全部引用锚点及报告快照哈希。

本轮修复了累计稿中$\hbar$的字体兼容问题：在生成的临时模板中使用已加载的AMS数学字形，章节公式和共用原模板不改。最终构建的输入哈希与当前章节、参考文献一致，报告位于未纳入版本控制的`build/manuscript_cn_full/build_report.json`。原三页摘要版继续保留。
