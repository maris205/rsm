# 黎曼标准模型：中文逐章修订

本稿按原文 *The Riemann Standard Model: A Dialogue Between Continuous and Discrete Worldviews* 的总体构建和第0—9章顺序逐章细化。保留原文的研究问题、架构设想和章节分工；在相应位置区分假设、条件推论、数值结果与物理解释。原小节将在进入对应章节时逐项核对。

“黎曼标准模型”（RSM）沿用原稿自拟名称，指这里提出的假设框架，不表示已经建立与粒子物理标准模型同等地位的理论。工作副标题为“基于离散同步演化的计算宇宙假说”。

**当前完成：** [中文摘要](abstract.md)、[第0章介绍](chapters/00_introduction.md)及两张[总览图](figures/captions.md)。[累计阅读版](paper.md)和[单栏PDF（11页）](paper.pdf)包含这些内容、目录及13条参考文献。原[摘要与图的三页版](front_matter.pdf)继续保留。第1—9章尚未在本目录逐章改写；下一步是第1章，不把章节计划计为已完成正文。

[摘要与原版的对应说明](notes/abstract_alignment.md)记录保留的想法、修订的结论强度、逆对数指数变化及原图处理安排。

[第0章对应说明](notes/chapter00_alignment.md)逐项说明原I—VIII部分如何融入已发表前作、内部时钟与格点模型、四类物理接口和现有数值结果。本轮继续采用中文，没有新增科学计算。

## 与此前七章稿的关系

[此前七章技术稿](../manuscript_v01/paper.md)及其[PDF](../manuscript_v01/paper.pdf)保持可读，作为已完成推导和计算的来源。其内部时钟、格点—连续场、算术残差对照将纳入原文相应章节；其七章结构不替代当前原文的十章结构。原文存档、阶段报告、原始图及数值结果均保留。

## 原文与中文稿逐章对应

以下中文章名保留原主题，具体题名和小节措辞随该章修订细化。每次先对照原章，再补入已发表前作、已有推导或计算；未完成的物理对应保留为待解决问题。

| 原文位置与主题 | 中文章名／主题 | 对应的现有材料与本轮状态 |
| --- | --- | --- |
| [前置材料／Abstract](../chapters/00_front_matter.md) | 标题、摘要与总体示意 | 本轮已写摘要、重画两张总览图。 |
| [第0章 Introduction](../chapters/01_chapter_00_introduction.md) | 介绍：从算术动力学到物理框架 | **已完成[中文稿](chapters/00_introduction.md)**；保留I—VIII顺序，纳入两篇已发表前作、7个核心公式及现有结果。 |
| [第1章 System Anomaly](../chapters/02_chapter_01_system_anomaly.md) | 系统异常：物理问题与研究动机 | **下一章**；保留常数、膨胀与真空问题，进入该章时更新来源和观测状态。 |
| [第2章 The Source Code](../chapters/03_chapter_02_source_code.md) | 源代码：算术结构与非自治动力学 | 已发表前作、算术输入核查、索引和时间的区别；保留逆对数平方核心假设。 |
| [第3章 Cosmic Lockstep Architecture](../chapters/04_chapter_03_cosmic_lockstep.md) | 形式化：宇宙同步架构 | 离散状态、同步规则、逻辑／空间描述；区分数值时间步与物理时间假设。 |
| [第4章 Mathematical Framework](../chapters/05_chapter_04_mathematical_framework.md) | 数学框架：辛演化与离散—连续映射 | 七章技术稿第2—4章、附录A／C；内部时钟反馈、能量收支与条件连续描述。 |
| [第5章 Observational Evidence](../chapters/06_chapter_05_observational_evidence.md) | 观测证据与约束：结构参数演化的检验 | 已有常数、原子钟与宇宙学项目；将响应假设、拟合能力和观测发现分别表述。 |
| [第6章 Predictions & Interpretations](../chapters/07_chapter_06_predictions_interpretations.md) | 预测与解释 | 对应原量子、真空、引力与波粒设想，逐项说明需要的机制和可能检验。 |
| [第7章 Riemann Engine Design v1.0](../chapters/08_chapter_07_engine_design.md) | 黎曼引擎设计 v1.0 | ECS、状态更新及现有计算代码；工程架构和物理本体假说分别定义。 |
| [第8章 Simulation Case Studies](../chapters/09_chapter_08_simulation_case_studies.md) | 模拟案例与系统验证 | 七章技术稿第5章、附录B及各实验报告；保留误差、未检出和失败记录。 |
| [第9章 Conclusion](../chapters/10_chapter_09_conclusion.md) | 结论：从几何描述到计算构想 | 汇总各章所得，原远期展望按可检验性逐项处理。 |
| [References](../chapters/11_references.md) | 参考文献 | 已建立[中文稿参考表](references.md)，按当前实际引用纳入13条；原表与技术稿书目继续作为来源。 |
| [Acknowledgments](../chapters/12_acknowledgments.md) | 致谢 | 完稿时按实际贡献与作者信息处理。 |

## 构建

从仓库根目录执行：

```bash
python manuscript_cn/figures/make_overview_figures.py
python manuscript_cn/build_front_matter.py
python manuscript_cn/build_manuscript.py
```

第一步重画总览图，不运行科学实验；第二步生成摘要与两图的三页版；第三步合并已完成章节及参考文献，生成累计Markdown和单栏PDF。正文修改入口为`abstract.md`、`chapters/`和`references.md`，图注在`figures/captions.md`。排版沿用已有单栏模板，依赖Python、Matplotlib、Pandoc、XeTeX及已安装字体；PDF检查使用PyMuPDF。中间文件分别进入未纳入版本控制的`build/manuscript_cn/`与`build/manuscript_cn_full/`。
