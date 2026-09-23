# 数学与物理基础：初步内部分析

本记录仅针对第 2、3、4、6、9 章的首次阅读，为下一阶段工作定位；不属于外部同行评审，也未改写章节。以下区分可直接验算的问题、尚缺定义的联系和工程隐喻，不据此判定整个研究方向。

值得保留的核心是：研究具有算术驱动参数的离散非自治动力系统，明确状态空间、更新规则及观测量，再检验其统计性质、稳定性和局部传播。这能形成可证明或反驳的数学问题。[§4.1.2 的辛条件](/root/autodl-tmp/dsc_world/riemann_model/chapters/05_chapter_04_mathematical_framework.md:17)和[§4.2.3 的近邻更新](/root/autodl-tmp/dsc_world/riemann_model/chapters/05_chapter_04_mathematical_framework.md:52)提供了起点；但 ECS、指针、渲染和系统老化，目前主要承担类比作用。

1. **算术索引尚未成为物理时间，且漂移方向写反。** [§2.4](/root/autodl-tmp/dsc_world/riemann_model/chapters/03_chapter_02_source_code.md:61)没有给出从整数或零点索引到宇宙时间的可检验映射，也没有从谱统计相似到原子能级、再到无量纲常数 \(\alpha\) 的算符关系。对正文公式，若 \(k>0,n>1\)，则 \(u_n=u_c-k/(\ln n)^2\) 单调增加并从下方趋于 \(u_c\)；连续延拓的导数为 \(2k/[n(\ln n)^3]>0\)。这与[§2.5](/root/autodl-tmp/dsc_world/riemann_model/chapters/03_chapter_02_source_code.md:97)“越来越偏离临界点”相反。可以衰减的是修正项，不能据此称 \(u_n\) 本身衰减。

2. **辛性不足以保证指定能量守恒或轨道有界。** [§4.1.3](/root/autodl-tmp/dsc_world/riemann_model/chapters/05_chapter_04_mathematical_framework.md:30)把体积守恒提升为能量及稳定性保证。直接反例为 \(M(q,p)=(2q,p/2)\)：它保持 \(dq\wedge dp\)，但对 \(H=(q^2+p^2)/2\)，有 \(H\circ M-H=3q^2/2-3p^2/8\)，且 \(q_0\ne0\) 时轨道无界。该反例针对“仅凭辛性”的推论；额外稳定性条件仍可研究。[§4.4.2](/root/autodl-tmp/dsc_world/riemann_model/chapters/05_chapter_04_mathematical_framework.md:90)又允许截断破坏辛性，需说明哪些结论仍有效。

3. **全局帧号不等于瞬时全局通信。** [§3.1](/root/autodl-tmp/dsc_world/riemann_model/chapters/04_chapter_03_cosmic_lockstep.md:9)的统一更新序号可以与局部规则并存；按 §4.2.3，传播依赖经过 \(m\) 步最多扩展 \(m\) 条邻接边。若[§3.3](/root/autodl-tmp/dsc_world/riemann_model/chapters/04_chapter_03_cosmic_lockstep.md:51)的零步更新能改变远处可读状态，就超出这一依赖域。须明确逻辑层是否影响可观测量、如何恢复相对论时钟关系，而不能由 Lockstep 名称推出因果一致性。

4. **共享指针没有建立纠缠测量规则。** [§6.2](/root/autodl-tmp/dsc_world/riemann_model/chapters/07_chapter_06_predictions_interpretations.md:31)未定义测量概率、不同测量方向或可控操作；[§9.3](/root/autodl-tmp/dsc_world/riemann_model/chapters/10_chapter_09_conclusion.md:39)却进一步声称可实现零延迟通信。在标准量子形式下，对 A 的局部完全正保迹操作且不筛选测量结果，有 \(\rho'_B=\operatorname{Tr}_A[(\mathcal E_A\otimes I)\rho_{AB}]=\rho_B\)，因此 B 的独立统计不变。需明确拟议模型修改了哪项前提，并给出新预测；指针地址和 \(O(1)\) 算法复杂度本身不提供物理传输时间。[Ghirardi 对无信号证明的重述，§4](https://arxiv.org/html/1305.2305)。

5. **普朗克单位的恒等式尚非光速推导。** [§4.2.3](/root/autodl-tmp/dsc_world/riemann_model/chapters/05_chapter_04_mathematical_framework.md:52)用 \(c=l_P/t_P\) 解释光速，但通常定义本已含 \(c\)：\(l_P=\sqrt{\hbar G/c^3}\)、\(t_P=\sqrt{\hbar G/c^5}\)。因此该比值是恒等式；若要独立推导，须先给出不以目标常数定标的长度、时间及测量对应。[NIST 2022 CODATA 表，第 1 页](https://physics.nist.gov/cuu/pdf/all.pdf)。另在 [§3.3](/root/autodl-tmp/dsc_world/riemann_model/chapters/04_chapter_03_cosmic_lockstep.md:51)，若 \(n\) 无量纲且 \(d_P\) 为长度，\(\Delta d_P/\Delta n\) 与速度的单位不一致。

6. **截断噪声没有推出暗能量，更没有推出无限能源。** [§6.3](/root/autodl-tmp/dsc_world/riemann_model/chapters/07_chapter_06_predictions_interpretations.md:53)把能量密度、哈密顿量差的范数和无量纲机器精度相等，缺少单位、体积归一化与系数。[§9.4](/root/autodl-tmp/dsc_world/riemann_model/chapters/10_chapter_09_conclusion.md:57)也未给出包含驱动场、导体、复位成本和输出功的能量收支；“垃圾回收失效”不能代替动力学。无限能量只能记为未建立的推测，并须与前文守恒主张核对。

建议下一步先制作符号与假设表，明确上述桥接关系，再挑选一个最小更新模型验证；当前材料不足以支持将常数漂移、超光速通信或无限能源表述为已推导结果。
