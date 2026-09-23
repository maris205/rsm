# 黎曼结构启发的逆对数平方动力学：内部时钟与微观—宏观连接

英文工作题名：*Riemann-Inspired Inverse-Log-Square Dynamics: An Internal-Clock Framework and Micro–Macro Connections*

2026-09-23开始逐章撰写，已形成中文Markdown正文与单栏PDF；英文版本待中文整稿审阅后推进。本文定位为假设性动力学框架及其条件检验。原始章节和v0.1—v0.5阶段短稿作为来源保留。本目录是连贯新稿，版本号不表示新增计算。

**当前已完成中文完整初稿：[单栏PDF](paper.pdf)（22页）与[Markdown连续阅读版](paper.md)。** 包含摘要、七章正文、三个附录、三张已有结果图及八篇参考文献。也可按下表逐章阅读。本版完成论述与证据的整合，未新增数值实验；后续集中改进关键物理与识别问题。

## 贯穿全文的问题

给定黎曼结构启发的逆对数平方响应，能否构造具有完整反馈的时钟—格点系统，从中得到有条件的连续场与晚时行为，并检验真实零点残差是否带来统计替代输入之外的响应？

主模型采用v0.4的内部坐标$R$、正则动量$P_R$与非线性格点；v0.5增加算术残差。v0.2的正则标量场与宇宙学背景属于不同候选，在讨论或附录中比较，不能将两者的方程和结论混用。平方指数、算术耦合与物理时间对应的假设地位贯穿全文。

## 章节结构与写作顺序

| 章 | 标题 | 要完成的论证 | 主要资料 | 状态 |
| --- | --- | --- | --- | --- |
| 1 | 引言：从算术谱到时间依赖动力学 | 交代谱统计背景、两篇前作、为何引入内部时钟，以及本文能回答的问题。 | [正文](chapters/01_introduction.md)、[参考文献](references.md) | 已写初稿 |
| 2 | 变量、尺度与基本假设 | 分开零点序号、映射步数、内部时钟、演化时间和宇宙时间；定义正则自由度、边界、响应函数及有效域。 | [正文](chapters/02_variables_and_assumptions.md)；[阶段变量表](../draft_v01/variables_and_claims_v01.md) | 已写初稿 |
| 3 | 自主时钟—格点动力学 | 给出共同Hamilton量、双向反馈与能流；推导纯包络分支的有界性及晚时逆对数平方尾律，列明假设。 | [正文](chapters/03_autonomous_dynamics.md)；[原推导](../experiments/log_clock_coupling_v01/reports/derivation_and_limits.md) | 已写初稿 |
| 4 | 从局部格点到连续有效场 | 推导长波方程和首个格距修正；区分形式一致性、解的收敛与统计粗粒化，说明均匀时钟的全局耦合。 | [正文](chapters/04_continuum_limit.md)；[v0.3](../draft_v03/micro_macro_bridge_v03.md) | 已写初稿 |
| 5 | 数值结果与算术输入对照 | 依次呈现固定系数基准、完整反馈、真实残差与替代输入；报告误差、未检出结果与匹配局限。 | [正文及图](chapters/05_numerical_tests.md)；[残差报告](../experiments/arithmetic_residual_v01/reports/results_cn.md) | 已写初稿 |
| 6 | 物理解释与可检验路线 | 说明连接原子、常数和宇宙学还缺哪些机制；给出能够限制或排除候选的检验，区分现有结果与未来目标。 | [正文](chapters/06_physical_tests.md)；[工作计划](../WORK_PLAN.md) | 已写初稿 |
| 7 | 结论 | 集中陈述模型构造、条件连接和本轮识别边界；保留后续物理推导的具体问题。 | [正文](chapters/07_conclusion.md) | 已写初稿 |

已附：[A，有限格点晚时命题的完整证明](chapters/appendix_a_asymptotics.md)；[B，数值复现、替代输入及精度与匹配失败记录](chapters/appendix_b_reproducibility.md)；[C，不同驱动候选与条件观测接口](chapters/appendix_c_physical_interface.md)。影响结论的限制同时保留在主文。

## 本轮资料与后续衔接

[第一章的来源与结果定位](notes/chapter01_sources.md)记录原始文献位置和已有实验依据，供写作复核使用，不属于论文正文。当前文献表只收录正文实际使用的条目；后续章节增加来源时再补充。

[完整初稿核对说明](notes/full_draft_status.md)记录本轮写作、公式与链接检查及独立代理复核的具体范围。

后续优先处理第6章列出的实际驱动匹配与检出能力，再研究局域驱动和真实物质接口。先保留这一完整初稿，新的正负结果均通过具体章节修订进入论文。

修改逐章源文件或文献后，从仓库根目录运行以下命令更新连续阅读版：

```bash
python manuscript_v01/assemble.py
```

`paper.md`由源文件生成；应修改各章、`abstract.md`或`references.md`，再重新合并。

生成单栏PDF：

```bash
python manuscript_v01/render_pdf.py
```

排版需Python 3、Pandoc、XeTeX及模板使用的LaTeX宏包、AR PL SungtiL GB／DejaVu Sans Mono／Tinos字体，以及用于成品检查的PyMuPDF。本环境没有`xelatex`命令，脚本使用已有`xetex`在`build/manuscript_typeset/`生成局部格式；不安装系统软件。所有中间文件和构建报告也写入该目录。成品包含可点击目录、文献锚点与仓库资料链接；[排版模板](typeset/paper.tex)和[脚本](render_pdf.py)均随仓库保存。
