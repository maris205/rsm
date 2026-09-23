# 最小作用量 v0.2：一手来源与适用边界

核查日期：2026-09-22。对象是候选标量场 $\chi$、指数势 $V_s e^{-2\chi}$ 及电磁规范动能函数 $B(\chi)$，不是完整文献综述或新的实验约束分析。本报告只核读以下四篇原始研究论文的相关部分；不采用论文发表时的观测数值作为 2026 年约束。

## 1. 四篇来源各支持什么

| 原始论文 | 实际核读位置 | 可以支持的范围 | 不能据此声称 |
|---|---|---|---|
| E. J. Copeland, A. R. Liddle, D. Wands, *Exponential potentials and cosmological scaling solutions*, **Phys. Rev. D 57, 4686–4690 (1998)**；[作者原稿](https://arxiv.org/pdf/gr-qc/9711068)、[DOI](https://doi.org/10.1103/PhysRevD.57.4686)。 | §II，式 (1)–(10)，Table I；§III，式 (11)。 | 在平直 FRW、正指数势、正则标量及常数状态方程的未耦合流体条件下，存在已知的标量主导解与跟踪流体的 scaling 解；后者的存在及吸引子条件为 $\lambda^2>3\gamma$，且 $\Omega_\phi=3\gamma/\lambda^2$。 | 条件不能直接移植到电磁源显著、背景跨时代变化或有额外耦合的模型；旧核合成数值界不等于今日界。 |
| H. B. Sandvik, J. D. Barrow, J. Magueijo, *A Simple Varying-alpha Cosmology*, **Phys. Rev. Lett. 88, 031302 (2002)**；[作者原稿](https://arxiv.org/pdf/astro-ph/0107512)、[DOI](https://doi.org/10.1103/PhysRevLett.88.031302)。 | 第 1–2 页，规范场变量重定义；式 (1)–(3) 及式 (3) 后对电磁源的说明。 | 通过标量依赖的电磁动能系数表示可变 $\alpha$，是已有的协变场论机制；电磁能量也会反过来驱动标量。自由纯辐射的平均 $F_{\mu\nu}F^{\mu\nu}=0$ 不等于物质内电磁贡献为零。 | 本候选的特定 $B(\chi)$、指数势或平方对数律不是该论文的结论；其当年拟合不能作为本候选通过观测的证据。 |
| G. Dvali, M. Zaldarriaga, *Changing alpha With Time: Implications For Fifth-Force-Type Experiments And Quintessence*, **Phys. Rev. Lett. 88, 091303 (2002)**；[作者原稿](https://arxiv.org/pdf/hep-ph/0108217)、[DOI](https://doi.org/10.1103/PhysRevLett.88.091303)。 | §2，式 (2)–(7)。 | 在所讨论的轻标量有效场论中，核子质量的 $\alpha$ 依赖产生标量—物质耦合，因成分不同而引出第五力及等效原理检验。 | 没有提供本候选的现代排除区域；不能把简单核子估计替代真实材料的质量与核结合能响应、场质量和环境解。 |
| T. Banks, M. Dine, M. R. Douglas, *Time-Varying alpha and Particle Physics*, **Phys. Rev. Lett. 88, 131301 (2002)**；[作者原稿](https://arxiv.org/pdf/hep-ph/0112059)、[DOI](https://doi.org/10.1103/PhysRevLett.88.131301)。 | §1，式 (1.1)–(1.7) 及后续真空能讨论。 | 在一般有效场论与所列自然性假设下，可变耦合的量子修正会改变标量势和真空能，保护宇宙尺度的轻标量需要额外机制或调参说明。 | 其截止尺度估计不是无条件的实验上限，也不是对所有紫外理论的排除；本候选的辐射稳定性仍须具体计算。 |

这些来源说明选用的是成熟的数学和场论构件。现阶段可能形成的贡献是**明确的条件组合、近似成立范围和联合观测关系**，不是首次提出指数势、可变 $\alpha$、自治系统的背景约化或标量场能量交换。

## 2. 对候选模型的独立代数对应

以下是本报告对候选式的代入核对，不是上述论文已经证明了本模型。使用 $c=\hbar=1$、度规 $(-+++)$，$M_{\rm Pl}=(8\pi G)^{-1/2}$ 为约化普朗克质量。令 $\chi$ 无量纲、$F_\chi>0$ 有质量量纲，正则场为 $\phi=F_\chi\chi$。标量部分为

$$
\mathcal L_\chi=-\frac{F_\chi^2}{2}(\partial\chi)^2-V_s e^{-2\chi},
\qquad \lambda=\frac{2M_{\rm Pl}}{F_\chi},\qquad V_s>0.
$$

忽略非引力源，且背景严格满足 $a\propto t^p$ 时，代入 $\chi=\ln(t/t_*)$ 得

$$
V_s t_*^2=\frac{F_\chi^2}{2}(3p-1),\qquad
\rho_\chi=\frac{3pF_\chi^2}{2t^2},\qquad
\Omega_\chi=\frac{F_\chi^2}{2pM_{\rm Pl}^2}.
$$

因此需要 $p>1/3$。对跟踪流体的分支，$p=2/(3\gamma)$，$\gamma=1+w$，

$$
\Omega_\chi=\frac{3\gamma F_\chi^2}{4M_{\rm Pl}^2},
\qquad \lambda^2>3\gamma.
$$

这是与第一篇论文 scaling 分支的参数对应；等号处不能直接沿用严格稳定结论。若采用标量主导分支，必须改用相应 Friedmann 条件，不能一面忽略标量反作用，一面把该分支当作背景的全部来源。

**同一纯指数势不产生全宇宙历史上任意背景的精确对数解。** 直接代入一般 $H(t)$ 要求 $3H(t)t-1$ 为常数。辐射期 $p=1/2$ 与物质期 $p=2/3$ 要求不同的 $V_s t_*^2$；同一 $V_s$、同一 $t_*$ 不能在两个时期均满足上述精确解。过渡期及晚期加速必须求解完整背景方程，或明确引入新的势形假设并重新检验。$t_*$ 在这条精确支上受势归一化关系约束，也未由此被确定为普朗克时间。

## 3. 规范归一化及被忽略的源

取

$$
\mathcal L_{\rm EM}=-\frac14 B(\chi)\mathcal F_{\mu\nu}\mathcal F^{\mu\nu},
\quad
B(\chi)=\left[1+\beta(\chi^{-2}-\chi_0^{-2})\right]^{-1},
\quad B(\chi_0)=1.
$$

在局部缓变背景上，正则规范场归一化给出 $e_{\rm eff}=e_0/\sqrt B$，故 $\alpha/\alpha_0=B^{-1}$。这是候选耦合的代数结果；$B>0$ 和 $F_\chi^2>0$ 只检查经典动能符号，不能证明所有扰动、物质响应与量子修正均可接受。空间依赖的规范场重定义还会产生导数项，不能把局部缓变归一化误写成全空间恒定重定义。相关场论构造的原始例子见 [Sandvik 等式 (1)–(3)](https://arxiv.org/pdf/astro-ph/0107512)。

该 $B$ 是**选择的响应函数**。指数势能产生对数场轨迹，但不选择 $\chi^{-2}$，也不把 $s=2$ 从可选的 $\chi^{-s}$ 中唯一挑出。平方指数的数论动机和物理作用量假设仍需分开陈述。

由上述作用量变分，均匀背景方程的源写为

$$
F_\chi^2(\ddot\chi+3H\dot\chi)+V_{,\chi}
=\mathcal S_\chi,
\quad
\mathcal S_\chi=-\frac14 B_{,\chi}\langle\mathcal F_{\mu\nu}\mathcal F^{\mu\nu}\rangle
+\left\langle\frac{\partial\mathcal L_m}{\partial\chi}\right\rangle.
$$

微观电磁描述与把结合能收入 $m_A(\chi)$ 的有效物质描述必须一致匹配，不能把同一结合能重复计算。即使宏观平均自由辐射源消失，原子核与物质的质量响应仍可能存在。应实际比较 $|\mathcal S_\chi|$ 与 $|V_{,\chi}|$ 或背景惯性／摩擦尺度；$\Omega_\chi\ll1$ 只约束对引力背景的贡献，不保证该源小。轻标量导致的成分响应机制见 [Dvali 与 Zaldarriaga §2](https://arxiv.org/pdf/hep-ph/0108217)。

## 4. 一个可继续检验的条件关系

定义正则归一化电磁耦合及今天的漂移

$$
d_{e0}=M_{\rm Pl}\left.\frac{\partial\ln\alpha}{\partial\phi}\right|_0
=-\frac{2M_{\rm Pl}\beta}{F_\chi\chi_0^3},
\qquad D_0=\left.\frac{\dot\alpha}{\alpha}\right|_0.
$$

**仅在今天也处于上述精确对数支的条件下**，

$$
d_{e0}=\frac{M_{\rm Pl}}{F_\chi}t_0D_0,
\qquad
\Omega_\chi d_{e0}^{\,2}=\frac{3\gamma}{4}(t_0D_0)^2.
$$

这是本候选消去参数的关系，没有输入新数据。它显示：在固定非零漂移下，缩小 $F_\chi$ 降低背景能量占比，会增大正则归一化电磁耦合。两种效果不能视作可彼此独立调小。实际晚期背景偏离 scaling 时，应使用 $D_0=(\partial_\chi\ln\alpha)_0\dot\chi_0$ 的一般关系重新计算。

若只有 $\alpha$ 这一通道，物体的标量电荷可写成 $q_A=d_{e0}Q_A$，其中 $Q_A=\partial\ln m_A/\partial\ln\alpha$。把它变成实验预言，还需要源与测试材料的 $Q_A$、标量的有效质量与作用程、边界条件及环境中的场解；不能从钟漂移一项直接宣布通过等效原理检验。[成分相关力的来源](https://arxiv.org/pdf/hep-ph/0108217)

## 5. v0.2 可以写下的结论与仍需完成的计算

若完整系统的作用量不显含坐标时间，背景场 $\bar\chi(t)$ 仍能随初值与动力学演化。将某条背景解代入探针哈密顿量后，探针就成为非自治系统；其能量变化来自与驱动场的交换。此处是标准背景约化的直接解释，不是违反守恒的新增机制。在膨胀时空中，应检验协变总能动量守恒及分量间交换，不能据此声称普遍存在一个全局恒定的宇宙总能量。

目前可以支持的表述是：**给定势、响应函数和近似条件，构造出一条通往平方对数响应的作用量路径，并得到应联合检验的背景—漂移—耦合关系。** 暂不能写成平方对数律被唯一推导、黎曼结构已获物理证明或全部实验界已满足。

接下来的必需项具有具体计算对象：

1. 完整 $H(t)$ 下的场轨迹、初值依赖和过渡误差，以及源项对轨迹的改动。
2. 指定材料与环境的标量电荷、作用程和等效原理预言；随后使用对应实验的一手最新结果给出允许参数区。
3. 明确有效场论截止、重整化条件和保护机制，估计 $\Delta V(\chi)$ 及 $\Delta m_\phi^2$。把势写成指数形式本身不保护它免受量子修正；这一问题在 [Banks、Dine、Douglas 的原论文](https://arxiv.org/pdf/hep-ph/0112059) 中已被提出。
4. 重新评估新增的 $F_\chi,V_s$、初值和响应参数是否实际减少自由度；形式上加入作用量不会自动改善现有数据对指数的辨识能力。

本轮没有计算第五力排除界、辐射修正或完整宇宙学似然，也没有新增真实观测证据。
