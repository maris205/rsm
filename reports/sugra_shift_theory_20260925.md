# 单时钟超引力候选：移位型 Kähler 势与横向稳定核查

日期：2026-09-25。材料范围：对本轮指定的一个中性手征场、一个带电手征对和一个规范多重态做解析审查；不加入新的稳定子，不拟合观测，不重写旧实验。以下是**条件模型的推导**。正横向曲率不是自然界已经选择该模型的证据，也不代表完成全部超引力宇宙学。

本轮得到一个清楚的候选：保持上一轮超势，改用沿时钟实方向具有移位对称性的 Kähler 势，可以保持指数时钟势的形状，并使另一实标量方向获得正树级质量。代价是明确选择了一种新的引力完成；带电粒子的物理质量和质量分裂也必须按照该完成重新计算。

## 1. 新增的结构与来源边界

令 $M_p$ 为约化普朗克质量，$F=F_\chi$，$\varepsilon=F^2/M_p^2$，并取

$$
Z=\frac F{\sqrt2}(\chi+iy),\qquad
W=W_0(Z)+M(Z)Q_+Q_-,
$$
$$
W_0=-\frac F{\sqrt2}\sqrt{U_i}
e^{-\sqrt2(Z-Z_i)/F},\quad Z_i=F\chi_i/\sqrt2,
\qquad
M(Z)=m_\infty\sqrt{1+\frac{\xi}{(\sqrt2 Z/F)^2}}.
\tag{S1}
$$

全纯平方根只在正实轴附近、避开零点和支点的既定分支中使用。质量律、超对称粒子谱和此非多项式超势都是新增 EFT 假设，尚未从黎曼动力学推导。

比较两种 Kähler 势（带电场的规范指数在这里略写）：

$$
K_{\rm can}=Z\bar Z+|Q_+|^2+|Q_-|^2,
\qquad
K_{\rm sh}=-\frac12(Z-\bar Z)^2+|Q_+|^2+|Q_-|^2.
\tag{S2}
$$

二者 $K_{Z\bar Z}=1$，实标量动能均为
$-F^2[(\partial\chi)^2+(\partial y)^2]/2$。这里只完成常 $\chi$ 动能分支，不能直接移植到原来常 $R$ 动能分支。

使用标准 $N=1$ 超引力 F 项

$$
V_F=e^{K/M_p^2}\left[K^{i\bar j}D_iW D_{\bar j}\bar W-
\frac{3|W|^2}{M_p^2}\right],\qquad
D_iW=W_i+K_iW/M_p^2.
\tag{S3}
$$

其来源及 Kähler 不变组合可核对 [Martin, *A Supersymmetry Primer*, §7.6，式 (7.6.17)—(7.6.24)](https://arxiv.org/html/hep-ph/9709356v7)。[Chiang 与 Murayama, *Building Supergravity Quintessence Model*, §2](https://arxiv.org/html/1808.02279v1)提供移位型 Kähler 构造的既有背景；该文的隐藏破缺区和势不同，因此本报告的具体势、稳定性与数值量级均由 (S1)—(S3) 重新推导，不借用其模型结论。

特别要分清：$K_{\rm sh}=K_{\rm can}+f+\bar f$，$f=-Z^2/2$，但真正的 Kähler 变换还要求

$$
W_{\rm sh}=e^{-f/M_p^2}W_{\rm can}
=e^{Z^2/(2M_p^2)}W_{\rm can}.
\tag{S4}
$$

本轮保持 $W$ 不变，所以 (S2) 给出两种不同的物理模型，并非把同一模型换了一个好看的规范。$W_0$ 和 $M$ 都打破了实方向的移位对称性；不能把所选 $K$ 说成全作用量的精确移位保护。

## 2. 全势、梯度与横向标量

在 $Q_\pm=0$，$U(\chi)=U_i e^{-2(\chi-\chi_i)}$，直接代入 (S3) 得

$$
V_{\rm sh}(\chi,y)=U e^{\varepsilon y^2}
\left(1-\frac32\varepsilon+\varepsilon^2y^2\right),
\tag{S5}
$$
$$
V_{\rm can}(\chi,y)=U e^{\varepsilon(\chi^2+y^2)/2}
\left[\left(1-\frac{\varepsilon\chi}{2}\right)^2
+\frac{\varepsilon^2y^2}{4}-\frac32\varepsilon\right].
\tag{S6}
$$

因此移位型在 $y=0$ 保留原指数形状，仅乘常数 $1-3\varepsilon/2$。它的精确梯度是

$$
V_{{\rm sh},\chi}=-2V_{\rm sh},\qquad
V_{{\rm sh},y}=2\varepsilon yUe^{\varepsilon y^2}
\left(1-\frac\varepsilon2+\varepsilon^2 y^2\right).
\tag{S7}
$$

设 $a_c=1-\varepsilon\chi/2$，在实轴上

$$
V_{{\rm can},\chi}|_0
=Ue^{\varepsilon\chi^2/2}
\left[(\varepsilon\chi-2)(a_c^2-3\varepsilon/2)
-\varepsilon a_c\right].
\tag{S8}
$$

以正则横向场 $Fy$ 定义质量：

$$
m_{y,\rm sh}^2=\frac{V_{{\rm sh},yy}|_0}{F^2}
=\frac{2U}{M_p^2}\left(1-\frac\varepsilon2\right),
\qquad
m_{y,\rm can}^2=\frac{Ue^{\varepsilon\chi^2/2}}{M_p^2}
(a_c^2-\varepsilon).
\tag{S9}
$$

在本轮 $\varepsilon=10^{-4}$，移位型的势为正且横向质量为正。更一般地，$0<\varepsilon<2/3$ 足以使 (S5) 在实 $y$ 上为正且以 $y=0$ 为横向最低点。正质量不意味着很快收敛：$m_y/H$ 必须沿真实背景计算，轻横向场可以保持很长时间的偏移。

场空间平直，实轴是无转弯的不变轨道。其线性横向扰动在树级满足

$$
\delta\ddot y_k+3H\delta\dot y_k+
\left(\frac{k^2}{a^2}+m_y^2\right)\delta y_k=0.
\tag{S10}
$$

这里 $V_y=\dot y=0$，所以该横向扰动在线性阶不产生背景应力的一阶扰动。式 (S10) 是这一玻色扇区的检验，不是全部纵向/引力/费米扰动的证明，也不代表任意远离实轴的初态均已排查。上一轮全局圈图的极小负横向项不能代替本轮完整超引力圈图；本轮应重新做量级预算，不能把旧 $\delta K$ 原样当作引力完成。

## 3. 重新计算局部带电质量谱

这一节保留 $K_0(Z,\bar Z)+|Q_+|^2+|Q_-|^2$，且 $K_{0,Z\bar Z}=1$。记 $w=W_0$，$A=D_Zw$，$C=M_Z+K_ZM/M_p^2$。把 (S3) 直接展开到 $Q$ 二阶，得到

$$
V=V_0+(s+d)(|Q_+|^2+|Q_-|^2)
+(\mathcal B Q_+Q_-+{\rm c.c.})+O(Q^4),
\tag{S11}
$$
$$
s=e^{K_0/M_p^2}|M|^2,\qquad
m_{3/2}^2=e^{K_0/M_p^2}\frac{|w|^2}{M_p^4},\qquad
d=m_{3/2}^2+\frac{V_0}{M_p^2},
\tag{S12}
$$
$$
\mathcal B=e^{K_0/M_p^2}
\left[A^*\left(M_Z+\frac{K_Z M}{M_p^2}\right)
-\frac{w^* M}{M_p^2}\right].
\tag{S13}
$$

混合项最后一项的负号来自两个 $|D_QW|^2$ 交叉项与 $-3|W|^2/M_p^2$ 的合并，不能只保留全局 $W_Z^*M_Z$。规范 D 势在 $Q=0$ 周围为四阶，不改变这里的二阶谱。于是

$$
m_\pm^2=s+d\pm|\mathcal B|,\qquad m_\psi^2=s.
\tag{S14}
$$

谱中包含两只复标量和一只 Dirac 费米子。需检查 $s+d>|\mathcal B|$、$H,m_{\rm probe}\ll\sqrt s$ 和绝热条件。本轮没有真实带电重粒子的背景占据。

一般软质量与物理双线性质量的结构可核对 [Brignole、Ibáñez 与 Muñoz, *Soft supersymmetry-breaking terms from supergravity and superstring models*, 式 (9)、(11)、(13)](https://arxiv.org/html/hep-ph/9707209v3)。该文讨论真空匹配；这里的瞬时系数是从 (S3) 独立展开，再在慢滚、重场极限使用，未把一般滚动背景强行当作驻定真空。特别是不能把仅在零真空能条件下的 $m_{3/2}^2=F_iF^i/(3M_p^2)$ 用于本轮 $V_0>0$ 背景；(S12) 是本模型使用的质量系数定义。

为便于实现，令

$$
r=\frac{\partial_\chi M}{M}
=-\frac{\xi}{\chi(\chi^2+\xi)},\qquad k_c=\varepsilon\chi^2/2.
$$

两支在实轴上的表达式分别为

| 量 | 移位型 $K_{\rm sh}$ | 标准型 $K_{\rm can}$ |
|---|---|---|
| $s$ | $M^2$ | $e^{k_c}M^2$ |
| $d$ | $(1-\varepsilon)U/M_p^2$ | $e^{k_c}U(a_c^2-\varepsilon)/M_p^2$ |
| $\mathcal B$ | $\sqrt{2U}M(r+\varepsilon/2)/F$ | $e^{k_c}\sqrt{2U}M[a_c(r+\varepsilon\chi/2)+\varepsilon/2]/F$ |
| $\partial_\chi\ln s$ | $2r$ | $\varepsilon\chi+2r$ |

所以本轮的局部质量阈值必须使用物理 $s$。相同全纯 $M(Z)$ 不保证两种引力模型的这一项贡献方向相同；完整超引力 $\alpha$ 尚须合并下节所述的同圈阶规范匹配，不能只由本表判定。并且 $d$ 随 $\chi$ 变化，不能沿用上一轮“共同软质量为常数”的导数程序。

## 4. 费米子、圈图与未计算项

带电方向满足 $D_{Q_\pm}W=0$、$\dot Q_\pm=0$；其二阶质量块没有与 $Z$ 或规范费米子的背景混合，也不进入本模型的 Goldstino 方向。因此 $\sqrt s=e^{K_0/(2M_p^2)}|M|$ 是重场局部阈值所用的 Dirac 质量。相反，$Z$ 费米子处在破缺方向，滚动背景还涉及运动学破缺与引力微子纵向模式；不能把它简单当作另一个独立的重费米子加进 CW 超迹。本轮未完成该费米子/引力微子扰动系统，也未给出引力微子产生或稳定性结论。

**局部质量阈值不是完整超引力规范耦合。** 本轮的 $\alpha$ 读数只保留固定匹配约定下、由所指定带电粒子物理质量产生的局部阈值贡献。全纯 Wilsonian 规范耦合与物理耦合的匹配还包含超Weyl、Kähler及Konishi异常相关项，参见 [Kaplunovsky–Louis, *Field Dependent Gauge Couplings in Locally Supersymmetric Effective Quantum Field Theories*, v2，§3.1—3.3，尤其式 (3.4)、(3.7)、(3.13)、(3.26)](https://arxiv.org/pdf/hep-th/9402005v2)。这些可在同一圈阶出现，不能由 $H/m\ll1$ 自动忽略；相关有限匹配、隐藏区与UV边界本轮均未完成。因此标准型候选在当前计算中的符号翻转只属于这一项局部阈值贡献，不是完整量子超引力 $\alpha$ 必然反向变化的结论。此次仅补充这一适用范围，不将尚未计算的异常项假定为零，也不假定其一定抵消或一定保留当前符号。

对于仅含所指定重带电粒子的、固定有限减法方案的局部平直空间势核，仍可定义

$$
V_1^{\rm heavy}=\frac{f(s+d+|\mathcal B|)+f(s+d-|\mathcal B|)-2f(s)}{32\pi^2},
\quad f(x)=x^2[\ln(x/\mu^2)-3/2].
\tag{S15}
$$

它的首项是

$$
V_1^{\rm heavy}
=\frac{sd}{8\pi^2}(L-1)+\frac{d^2+|\mathcal B|^2}{16\pi^2}L+\cdots,
\quad L=\ln(s/\mu^2).
\tag{S16}
$$

必须连同 $s',d',\mathcal B'$ 求力；固定点常数减法不能删掉力。主阶 $s^2$ 仍相消，但该重场谱满足

$$
\operatorname{Str}\mathcal M^2=4d,\qquad
\operatorname{Str}\mathcal M^4=8sd+4d^2+4|\mathcal B|^2.
\tag{S17}
$$

因此有限阈值很小不等于所有紫外敏感项已经消失。还未包括曲时空局部项、全部超引力圈图、匹配反项和可见区的破缺传递。曲率项可有 $m^2\mathcal R/(16\pi^2)$ 的维数规模，可能涉及普朗克项和导数作用量的重整化；$H^2/s\ll1$ 说明重场的局部展开可用，**不单独保证这些项相对于极小的时钟势也可忽略**。若展示该量级预算，应标为没有计算精确系数的 EFT 诊断，不能写成完整一圈超引力结果。

原平直全局模型的 $\delta K$ 不能未经匹配直接移植，更不能既用其 $U/G$ 势，又重复加入相同阶的完整 CW 修正。若所选重场的势预算远低于双精度背景求解误差，合理产物是树级轨迹加独立高精度量级预算，不声称积分器分辨出了该微小修正。

## 5. 对宇宙学完成度的准确表述

把本树级玻色扇区放入 Einstein–FLRW 背景并外加辐射、尘埃和 $\Lambda$，能做一个定义清楚的背景测试。它还不是所有宇宙成分的局部超对称完成。尤其：

- (S12) 中的 $V_0$ 是由本 $K,W$ 得到的超引力势，不能未经模型推导把外加流体密度或 $\Lambda$ 放进 $d$。
- 若将来用隐藏超场产生 $\Lambda$ 或现实破缺尺度，其 F 项、交叉项、重整化和可见区耦合会改变本轮谱；当前结果不能替代那一步。
- 相同初始 $\chi$ 和物理速度 $\dot\chi$ 在两种势下给不同的初始 $H$，因此比较时应记录 $\dot\chi/H_{\rm ref}$，不要把相同 $d\chi/dN$ 误当相同物理速度。
- 此构造仍未推出 $1/\chi^2$ 质量律、$\chi\sim\ln t$ 的适用区间、真实电磁物质的反馈或所有粒子实验限制；也没有新增观测支持。

本轮最有根据的结论是：**存在一个简单的树级超引力候选，同时保持时钟指数势的形状并给出正横向标量质量；局部带电质量阈值贡献仍可计算，但需使用新的物理质量和质量分裂。** 可见区隔离、隐藏破缺区、引力微子动力学以及包括同圈阶异常项在内的完整规范匹配仍是独立问题。

复核材料：[23项符号恒等式脚本](sugra_shift_theory_20260925_checks.py)、[结果](sugra_shift_theory_20260925_checks.json)、[来源身份与读取范围](sugra_shift_theory_20260925_sources.json)。执行命令为 `python reports/sugra_shift_theory_20260925_checks.py`；23/23通过，退出码0。检查覆盖从 $D_ZW$ 得到全势、带电二阶展开、横向质量、Kähler变换和物理质量导数。它们是模型内部代数核验，不是外部同行评审或物理真实性认证。
