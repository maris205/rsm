# 黎曼结构启发的假设性动力学框架

当前写作方向：**按原文第0—9章构建中文逐章修订稿，已完成摘要、第0—5章及十一幅矢量图**，入口见[manuscript_cn](manuscript_cn/README.md)。此前摘要、七章正文及三个附录的完整稿作为技术材料保留；其结构不替代原文十章。计算已完成至真实零点残差与替代输入的条件比较v0.5。写作更新日期为2026-09-23；最近一轮计算完成于2026-09-22。当前成果是明确的条件模型与可复现计算，尚未获得常数演化的观测确认。

**论文定位：** 以黎曼／算术结构为动机，以逆对数平方响应为核心假设，通过明确的 Hamilton 模型研究微观自由度到宏观有效响应的连接。写作中始终区分输入假设、条件推论和待验证的物理对应；内部一致性的依据限于已推导、已核验的部分，不把“尚未发现矛盾”当作物理成立的证据。当前未检出特殊算术响应及已知数值、对照局限均属于需要报告的结果。

**正文审阅：** 已完成[摘要及第0—5章首轮审阅](reports/manuscript_cn_body_review_20260923.md)，基准为`ac39f8c`。本轮未定位核心数学链中的实质错误；优先问题是Hubble校准扩展的公开复现入口，其次是贡献与候选模型关系、最小方法说明、符号及章节重复。报告区分修订建议与后续研究，第0—5章的整体修订尚未执行；新增第6—7章已注意贡献分级、候选关系和实现边界，没有重跑科学计算。

项目仓库：[maris205/rsm](https://github.com/maris205/rsm)。使用SSH克隆后，下文新列的命令从仓库根目录执行：

```bash
git clone git@github.com:maris205/rsm.git
cd rsm
```

论文、原稿、实验代码、输入和已有结果一并保存；`build/`、Python缓存及Notebook检查点不纳入版本控制。两个早期跨项目数据依赖已整理成本地副本，见[复现路径与来源说明](reports/repository_portability.md)。部分历史报告仍保留原工作区的相邻项目引用，属于背景资料；已有运行的哈希和审计记录保留其原始版本含义。

原稿题为 *The Riemann Standard Model: A Dialogue Between Continuous and Discrete Worldviews*。
这里保留它作为研究构想的原始记录；“Standard Model”目前是原稿自拟名称，并不表示已建立公认或经过观测确认的标准模型。

当前新稿：

- **[按原文逐章修订：中文累计稿](manuscript_cn/paper.md)**／[PDF（60页）](manuscript_cn/paper.pdf)：已完成摘要、[第0章介绍](manuscript_cn/chapters/00_introduction.md)、[第1章物理动机](manuscript_cn/chapters/01_system_anomaly.md)、[第2章算术来源](manuscript_cn/chapters/02_source_code.md)、[第3章同步架构](manuscript_cn/chapters/03_cosmic_lockstep.md)、[第4章数学框架](manuscript_cn/chapters/04_mathematical_framework.md)、[第5章观测约束](manuscript_cn/chapters/05_observational_evidence.md)、[第6章预测与解释](manuscript_cn/chapters/06_predictions_interpretations.md)、[第7章引擎设计](manuscript_cn/chapters/07_engine_design.md)、十二幅矢量图和51条参考文献。第7章保留原7.1—7.5及10个三级小节，对照实际代码给出驱动、同步积分、读出和验证流程，纳入已有线性终场性能基准，明确实现范围；已有观测结果、未检出与失败记录继续保留。对应说明见[第0章](manuscript_cn/notes/chapter00_alignment.md)、[第1章](manuscript_cn/notes/chapter01_alignment.md)、[第2章](manuscript_cn/notes/chapter02_alignment.md)、[第3章](manuscript_cn/notes/chapter03_alignment.md)、[第4章](manuscript_cn/notes/chapter04_alignment.md)、[第5章](manuscript_cn/notes/chapter05_alignment.md)、[第6章](manuscript_cn/notes/chapter06_alignment.md)、[第7章](manuscript_cn/notes/chapter07_alignment.md)及[章节总表](manuscript_cn/README.md)。下一步进入第8章模拟案例与系统验证。
- **[完整论文初稿：单栏 PDF](manuscript_v01/paper.pdf)**（22页）与[Markdown连续阅读版](manuscript_v01/paper.md)：已整合摘要、七章、三个附录及三张结果图；也可按[章节目录](manuscript_v01/README.md)阅读。[参考文献](manuscript_v01/references.md)包含两篇已发表前作。以下v0.1—v0.5作为分阶段来源保留。
- **[研究主线备忘：从微观非线性动力学到宏观理论](RESEARCH_DIRECTION.md)**：记录用户的长期方向，对照原第4、7章等已有设想，明确辛结构、长波极限与统计粗粒化的不同作用。
- **[假设性框架与算术输入检验 v0.5](draft_v05/arithmetic_response_v05.md)**：最近一轮计算的阶段短稿，开篇说明假设、条件结果和物理对应的区别。实际使用前作零点缓存，计算1593例；八个真实残差块的主指标均落入所构造替代集合的中心范围，未检出特殊宏观响应。节点及实际驱动匹配的局限明确保留。[实验入口](experiments/arithmetic_residual_v01/README.md)、[完整结果](experiments/arithmetic_residual_v01/reports/results_cn.md)、[主比较图](experiments/arithmetic_residual_v01/figures/arithmetic_primary_comparison.png)和[响应热力图](experiments/arithmetic_residual_v01/figures/arithmetic_response_and_controls.png)。
- **[逆对数驱动与反馈 v0.4](draft_v04/log_clock_bridge_v04.md)**：把核心平方对数响应实际接入同一Hamilton体系，计算反作用、能流与连续场，并给出有条件的晚时尾律。[实验入口](experiments/log_clock_coupling_v01/README.md)、[数值报告](experiments/log_clock_coupling_v01/reports/results_cn.md)、[响应图](experiments/log_clock_coupling_v01/figures/log_clock_response.png)和[场演化图](experiments/log_clock_coupling_v01/figures/log_clock_field_bridge.png)。**真实黎曼物理来源仍待推导**，见该阶段的[算术桥接状态](experiments/log_clock_coupling_v01/reports/arithmetic_bridge_status.md)。
- **[格点到连续场 v0.3](draft_v03/micro_macro_bridge_v03.md)**：首个非线性微观—宏观桥接实例；[实验入口](experiments/micro_macro_v01/README.md)、[结果报告](experiments/micro_macro_v01/reports/results_cn.md)、[时空演化图](experiments/micro_macro_v01/figures/micro_macro_evolution.png)及[收敛图](experiments/micro_macro_v01/figures/micro_macro_convergence.png)。基础空间误差约二阶，首个格距修正后约四阶；保留最细参考的一项精度提示。
- **[驱动场与条件结论 v0.2](draft_v02/action_framework_v02.md)**：指数势标量＋电磁响应的作用量、场方程和能量收支；逆对数平方的成立条件、晚期偏离及与非自治探针的衔接。它与v0.4的自由内部时钟是两个不同候选。
- **[自主背景实验](experiments/action_completion_v02/README.md)**：正向求解尘埃＋Λ＋标量。在示例中，固定当前漂移后的高红移传递与精确对数假设相差约16.7%；这是预测之差，没有新增观测拟合。附[图](experiments/action_completion_v02/figures/background_completion.png)、[数值报告](experiments/action_completion_v02/reports/numerical_results_cn.md)及[独立复核](experiments/action_completion_v02/reports/independent_background_review.md)。
- v0.2的[变分与能量推导](reports/action_completion_v02_derivation.md)、[独立数学核对](reports/action_completion_v02_independent_review.md)、[原始文献与限制](reports/action_completion_v02_sources.md)。
- **[最小公设 v0.1](draft_v01/minimal_framework_v01.md)**：五条公设、明确的非自治更新、数学条件，以及α与原子钟的联合观测关系。
- **[变量、假设与证据表](draft_v01/variables_and_claims_v01.md)**：定义与量纲、旧符号修正、哪些是推论、哪些是新增假设。
- [数学推导与边界](reports/minimal_dynamics_v01_review.md)、[α接口核对](reports/alpha_interface_v01_review.md)、[前作及引用审计](reports/sources_v01_audit.md)。
- **[模拟恢复实验](experiments/alpha_recovery_v01/README.md)**：7种设计、133种注入情形；[结果报告](experiments/alpha_recovery_v01/reports/recovery_results_cn.md)及[科研图](experiments/alpha_recovery_v01/figures/captions.md)。结果是幅度可按预期恢复，但当前设计不能辨认平方指数，联合幅度信息几乎全来自原子钟。

原稿和整体计划：

1. **[章节目录](chapters/README.md)**：英文原文，按第0—9章拆分，另有前置材料、参考文献和致谢。
2. **[初读梳理](reports/first_reading_cn.md)**：整体思路、各章用途和首先需要解决的一致性问题。
3. **[下一阶段工作计划](WORK_PLAN.md)**：先做最小定义与可检验模型，再逐步连接四个观测方向。

补充材料：

- [已有四个项目的证据与接口](reports/prior_projects_evidence.md)
- [数学与物理基础的初步检查](reports/foundations_first_pass.md)
- [转换范围、已知原稿问题与校验说明](reports/conversion_notes_cn.md)
- [机器可读转换清单](reports/conversion_manifest.json)
- [原始 Word 文档](ori_paper/merged_document_v2_with_latex.docx)

## 已完成的整理

共生成 **13个章节文件**，保留134个原始标题、335处Word公式（315处行内、20处独立公式）、5张表和9处正文图片引用。图片字节与Word内嵌资源一致；原文件哈希未变。原稿的英文表述、公式、旧拟合结论和推测均未在转换时改写。科学意见单列在报告中，便于以后逐条审定。

v0.1已另建新稿，保留逆对数平方的核心设想，明确时间桥接、有效探针动力学和电磁响应的假设地位。小矩阵核验28项通过，检查辛性、冻结稳定域、数值相位、外功和局部依赖；它们是数学实现检查，没有新增观测拟合，也没有证明完整统一理论。既有项目提供模型、工具和约束，尚未共同验证同一个物理机制。

第一轮模拟沿用真实观测的红移与误差，只生成合成响应。主实现183项检查及独立统计路径207项检查通过，附绝热与数值相位预算。检出注入信号不代表在真实数据中发现变化；原始观测均值未被用作模拟真值。

v0.2进一步给出一个自主演化的驱动候选；在单一流体标度支上得到对数场轨迹，在加入Λ后计算实际偏离。29项主数值检查及24项独立输出比较通过。响应函数的平方指数、真实物质源与微观原子响应仍待说明，数值示例也未检验第五力、辐射背景或量子修正。v0.1的精确时间律和v0.2的动力学延拓是有共同极限的两种候选，不能把后者结果直接套到旧拟合。

v0.3以固定参数的局部非线性正则格点为基准，分别求解格点与连续场，并区分空间差异、时间误差和高频失效。86项格点检查、22项参考检查、20项比较诊断通过；另一次最细网格精度审计4项通过、1项未通过，主误差对收紧容差仅变化0.00781%。这一成果支持指定初值与有限时段的条件桥接，尚未从算术规则导出真实微观Hamilton量。

v0.4实际接入逆对数平方响应，完整计算时钟—格点反馈，并与规定外部驱动和冻结系数比较；88项主核验、27项独立参考检查、29项比较诊断通过，另有独立数学与比较代码复核。所选正能量分支下，场向时钟释放能量，有限时段偏离外部对数律，但固定有限格点的晚时主导尾律可以保留。响应幂仍是预设输入，实际零点／算术残差尚未输入；数值检查不等于物理确认，也不覆盖上一轮已记录的精度边界。

v0.5随后实际加入零点间距残差：1024个评价间距、八个固定块，每块99个IAAFT及99个配对精确节点谱对照。终点能量交换汇总经验秩为0.40和0.36，没有提供特殊算术作用的证据；不是等效性结论，也不是物理$p$值。主数值17项、最终参考5类、比较24项核验通过，四个节点匹配失败及初始参考的两项精度失败原样保留。窗化插值后的实际驱动谱未严格匹配，下一轮需要改进这一控制并检查检出能力。

## 目录

```text
riemann_model/
├── ori_paper/       原始Word，保持原样
├── chapters/        按章Markdown及原图assets/
├── manuscript_cn/   按原文第0—9章逐章修订；当前已写摘要、第0—7章及十二幅图
├── manuscript_v01/  完整初稿（Markdown/PDF）、逐章正文及合并与排版脚本
├── draft_v01/       最小公设短稿、变量与主张表
├── draft_v02/       驱动场作用量、条件结论与数值结果
├── draft_v03/       非线性微观格点到连续场的条件桥接
├── draft_v04/       逆对数响应、时钟—格点反馈及条件晚时结果
├── draft_v05/       真实算术残差的条件比较与未检出结果
├── experiments/     注入恢复、自主背景、格点与驱动反馈的独立计算
├── reports/         初读、推导核对、来源与数值核验记录
├── code/            转换脚本与最小动力学核验
├── build/           转换中间文件，可重新生成
└── WORK_PLAN.md     待开展的研究顺序与验收条件
```

## 重新转换

需要Python 3及Pandoc。本次使用Pandoc 2.9.2.1；无额外Python包依赖。在仓库根目录运行：

```bash
python code/split_original_docx.py
```

脚本会重新生成13个原稿章节和转换报告，覆盖这些生成文件。后续实质修订应放入独立的新稿目录，保留本目录中的原稿转换版用于对照。新的章节目录提供可用的文件链接；原稿前置材料中的旧Word页码和目录字段仅作存档。

## 最小动力学核验

需要Python 3和NumPy。在仓库根目录运行：

```bash
python code/check_minimal_dynamics_v01.py
```

输出见[核验JSON](reports/minimal_dynamics_v01_checks.json)。这些检查不包含真实对数驱动的长期稳定证明、原子结构计算或观测似然。
