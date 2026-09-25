# 最小带电标量完成：作用量、量子阈值与宇宙反馈核查

日期：2026-09-25。作者侧独立解析复核；没有运行宇宙生产积分，没有新增观测拟合，未更改正文或旧实验。本记录延续[常数—共轭时钟备忘录](constant_clock_bridge_20260925/README.md)，区分结构动机、物质匹配假设和可以推导的结果。

**结论：单个带电复标量足以把指定的质量变化传递到低能精细结构常数，但同一阈值产生的有效势通常远大于慢变宇宙时钟的尺度。只保留电磁响应，会遗漏主导反作用。** 此完成适合做明确的量子反馈对照；它还不足以自然解释常数漂移。

## 1. 作用量及适用范围

采用自然单位、度规 $(-+++)$，$\chi=\ln(R/R_*)>0$ 为无量纲变量。新增电荷 $Q=1$ 的复标量 $X$：

$$
S=\int\sqrt{-g}\,d^4x\left[
\frac{M_{\rm Pl}^2}{2}\mathcal R
-\frac{F_\chi^2K(\chi)}2(\partial\chi)^2-U_{\rm tree}(\chi)
-|DX|^2-s(\chi)|X|^2-\frac{\lambda_X}{2}|X|^4
-\frac14F_{\mu\nu}F^{\mu\nu}\right]+S_{\rm other}+S_{\rm ct},
$$

$$
D_\mu X=(\nabla_\mu-ieQA_\mu)X,\qquad
s(\chi)=m_X^2(\chi)=m_\infty^2(1+\xi/\chi^2).
$$

要求 $F_\chi^2K>0$、$s>0$、$\lambda_X\geq0$，沿 $X=0$ 的未破缺电磁支计算。非多项式质量函数是EFT输入，不是已证明的紫外完成。现有实格点没有导出电荷、复标量及其物理质量；这些是新增假设。$\lambda_X$ 不出现在该支的重场二次算符中，本次重场一圈匹配不等于已经包括轻场圈图和所有高阶反馈。

局域展开要求 $H/m_X\ll1$、$|\dot m_X|/m_X^2\ll1$、相关外动量 $p/m_X\ll1$；全区间没有阈值交叉。高能耦合边界和其他阈值保持固定。重整化能标 $\mu$ 不是宇宙时间的倒数。

以下采用低曲率的局部平直一圈近似。完整曲时空量子EFT还需登记 $\mathcal R|X|^2$、$m_X^2\mathcal R$ 和曲率反项，不能将下面的CW势与动能项称为完整曲时空有效作用量。标量曲率耦合与引力反项的明确原始例见 Markkanen、Tranberg [arXiv:1303.0180v2，式(2.5)—(2.6)、附录E](https://arxiv.org/pdf/1303.0180v2)。$m_X\ll M_{\rm Pl}$ 时诱导引力项可按 $m_X^2/(16\pi^2M_{\rm Pl}^2)$ 等估量；这不消除有效势的尺度问题。

## 2. 原始来源与一圈系数

核对来源为 Henning、Lu、Murayama，*How to use the Standard Model effective field theory*，arXiv:1412.1837v2（2015-07-08）。[原始PDF](https://arxiv.org/pdf/1412.1837v2) 式(2.54)、其后的约定及式(2.55)给出复标量因子 $c_s=1$、质量势／二导数项及规范beta系数。定位为PDF第28、30页（印刷页26、28）。以下转成本文约定，并独立检查二导数归一化。

单复标量贡献为

$$
\beta_e=\frac{Q^2e^3}{48\pi^2},\qquad
\frac{d\alpha}{d\ln\mu}=\frac{b_{\rm em}\alpha^2}{2\pi},\qquad
b_{\rm em}=Q^2/3=1/3.
$$

$\overline{\rm MS}$ 重场一圈势及动能修正为

$$
U_1=\frac{s^2}{32\pi^2}\left(\ln\frac{s}{\mu^2}-\frac32\right),\qquad
G=F_\chi^2K+\Delta G,\qquad
\Delta G=\frac{(s_{,\chi})^2}{96\pi^2s}>0.
$$

固定 $\mu$，记 $L=\ln(s/\mu^2)$，直接求导：

$$
U_{1,\chi}=\frac{s s_{,\chi}}{16\pi^2}(L-1),\qquad
U_{1,\chi\chi}=\frac{s_{,\chi}^2L+s s_{,\chi\chi}(L-1)}{16\pi^2},
$$

$$
s_{,\chi}=-2m_\infty^2\xi\chi^{-3},\qquad
s_{,\chi\chi}=6m_\infty^2\xi\chi^{-4}.
$$

二导数系数的独立推导：任意固定 $s>0$ 附近，单复标量的欧氏作用量为 $\Gamma_E={\rm Tr}\ln(-\partial_E^2+s+\delta s)$，二次项为 $-\delta s(p)\delta s(-p)I(p)/2$。泡积分满足

$$
I(p)-I(0)=-\frac{p^2}{16\pi^2s}\int_0^1x(1-x)dx+O(p^4)
=-\frac{p^2}{96\pi^2s}+O(p^4).
$$

于是作用量含 $+(\partial\delta s)^2/(192\pi^2s)$，相应动能系数为上式 $\Delta G$；回到 $(-+++)$ 写为 $-\Delta G(\partial\chi)^2/2$。实标量行列式会多一个 $1/2$。此处仅保留一圈、二导数阶，不代表完整的非局域作用量。

参考点 $B(\chi_0)=1$、固定参考电荷 $\alpha_0=e_0^2/(4\pi)$ 时，阈值给

$$
B(\chi)=1-\frac{\alpha_0}{12\pi}\ln\frac{s(\chi)}{s(\chi_0)},\qquad
\alpha(\chi)=\frac{\alpha_0}{B(\chi)},\qquad B>0.
$$

保留 $1/B$ 是匹配式的代数求逆，不等于已完成高圈重求和。小变化的响应与同一时间参数下的漂移为

$$
\frac{\Delta\alpha}{\alpha_0}\simeq
\frac{\alpha_0m_\infty^2\xi}{12\pi s_0}(\chi^{-2}-\chi_0^{-2}),\qquad
\frac{\dot\alpha}{\alpha}=-\frac{\alpha m_\infty^2\xi}{6\pi s\chi^3}\dot\chi.
$$

$\xi>0$、$\dot\chi>0$ 时该支的 $\alpha$ 下降。平方指数来自指定质量律，圈积分传递了它，没有选择指数或导出宇宙历史。

## 3. FLRW闭合与交换符号

定义 $V=U_{\rm tree}+U_1+U_{\rm ct}$。齐次宇宙时间 $t$ 下

$$
G(\ddot\chi+3H\dot\chi)+\frac{G_{,\chi}}2\dot\chi^2+V_{,\chi}+S_\chi=0,
$$

$$
\rho_\chi=G\dot\chi^2/2+V,\qquad p_\chi=G\dot\chi^2/2-V,\qquad
\dot\rho_\chi+3H(\rho_\chi+p_\chi)=-S_\chi\dot\chi.
$$

无真实 $X$ 占据及经典电磁源时可设 $S_\chi=0$；真空势 $U_1$ 仍存在。若保留非相对论真实粒子，$\rho_X=n_Xm_X$、$\dot n_X+3Hn_X=0$，则

$$
S_X=n_Xm_{X,\chi}=\rho_X\partial_\chi\ln m_X,\qquad
\dot\rho_X+3H\rho_X=+S_X\dot\chi.
$$

真空圈与真实粒子占据可同时存在，但不得重复计入同一真空行列式；非绝热粒子产生不在本次近似内。一般有压强粒子源还应采用 $(\rho_X-3p_X)\partial_\chi\ln m_X$。

显式电磁源为

$$
S_{\rm EM}=\frac{B_{,\chi}}4\langle F^2\rangle,\qquad
\dot\rho_{\rm EM}+4H\rho_{\rm EM}=+S_{\rm EM}\dot\chi.
$$

最后一式要求各向同性平均、无电流作功；单一相干电磁场通常不适合作严格各向同性背景。局部正交标架有 $F^2=2(|\mathbf B_{\rm phys}|^2-|\mathbf E_{\rm phys}|^2)$，因此当 $B_{,\chi}>0$，纯电源把标量推向增大的 $\chi$，纯磁源反向。自由辐射等分可给平均 $F^2=0$，但这不消除物质电磁结合能的源；已计入物质质量的结合能不能再次重复计数。存在电流作功时，物质端与电磁端分别加 $+\langle\mathbf E\cdot\mathbf j\rangle$ 和其负值。

引力方程为

$$
3M_{\rm Pl}^2H^2=\rho_\chi+\rho_m+\rho_r+\rho_o,\qquad
-2M_{\rm Pl}^2\dot H=G\dot\chi^2+\rho_m+4\rho_r/3+\rho_o+p_o.
$$

所有交换项相消后得到总连续方程。真空常数可放入 $V$ 或 $\rho_o$，不得双计数。固定 $H(t)$ 而不回代标量能量的计算须称为测试场近似。

令 $N=\ln a$、$v=d\chi/dN$，数值方程可写成

$$
\frac{dv}{dN}=-\left(3+\frac{d\ln H}{dN}\right)v
-\frac{G_{,\chi}}{2G}v^2-\frac{V_{,\chi}+S_\chi}{GH^2},\qquad
H^2=\frac{\rho_m+\rho_r+V+\rho_o}{3M_{\rm Pl}^2-Gv^2/2}.
$$

后式要求所选物理解支及分母条件成立，不允许任意指定 $v(N)$ 后声称获得自主解。

## 4. 两个动能分支不同

- $K=1$：树级正则场为 $\phi=F_\chi\chi$。给定指数势可在适当标度支出现 $\chi\propto\ln t$；混合宇宙成分下仍须解实际轨迹。
- $K=e^{2(\chi-\chi_{\rm ref})}$：令 $R=R_{\rm ref}e^{\chi-\chi_{\rm ref}}$，树级动能恰为 $F_\chi^2\dot R^2/(2R_{\rm ref}^2)$，保留原 $R$ 的常系数动能。标量方程多出 $\dot\chi^2$ 联络项，不能照搬 $K=1$ 的轨迹。

第二支无势、无源时满足 $\ddot R+3H\dot R=0$；固定物质主导测试背景给 $R=R_\infty-C/t$，不是 $R\propto t$。坐标恒等式不替代时间对应。

任意 $G>0$ 可局部定义 $d\phi=\sqrt G\,d\chi$，正则势曲率为

$$
V_{,\phi\phi}=\frac{V_{,\chi\chi}}G-\frac{G_{,\chi}V_{,\chi}}{2G^2}.
$$

这尚不是含物质／引力扰动的所有模的完整质量谱。

## 5. 辐射尺度与局部减法

下面是指定方案中**没有保护或精细抵消时的量级诊断**，不是方案无关的排除定理。固定 $\mu=m_0$，参考点附近

$$
\Delta U_1\simeq-\frac{m_0^4}{8\pi^2}\Delta\ln m_X,\qquad
\frac{\Delta\alpha}{\alpha_0}\simeq\frac{\alpha_0}{6\pi}\Delta\ln m_X,
$$

$$
|\Delta U_1|\simeq\frac{3m_0^4}{4\pi\alpha_0}
\left|\frac{\Delta\alpha}{\alpha_0}\right|.
$$

以项目 $H=67.4$ km s$^{-1}$ Mpc$^{-1}$、$M_{\rm Pl}=2.435\times10^{27}$ eV、$\alpha_0^{-1}=137.035999084$ 作参考标尺，并非本次重估这些量。$3M_{\rm Pl}^2H^2=3.6768\times10^{-11}$ eV$^4$。假设响应为 $10^{-6}$，则 $\Delta\ln m_X\simeq0.002583$：

| 尺度示例 $m_0$ | $|\Delta U_1|$（eV$^4$） | 相对 $3M_{\rm Pl}^2H^2$ |
|---|---:|---:|
| 1 eV | $3.27\times10^{-5}$ | $8.90\times10^5$ |
| 1 MeV | $3.27\times10^{19}$ | $8.90\times10^{29}$ |
| 100 GeV | $3.27\times10^{39}$ | $8.90\times10^{49}$ |

这些是展示 $m^4$ 放大的量纲示例，不是获实验许可的粒子候选；1 eV粒子也不满足所有光谱测量的重阈值条件。

另取 $m_0=100$ GeV、$\chi_0=5$、$\xi=1$、$F_\chi=M_{\rm Pl}$、$K=1$：$\Delta G/F_\chi^2\simeq4.21\times10^{-40}$，而未减法势给 $|U_1''/F_\chi^2|^{1/2}/H\simeq2.18\times10^{25}$。该点曲率为负；这里报告的是其绝对尺度，并非稳定振子质量。“动能修正小”不等于“量子反作用小”。

局部减法

$$
U_1^{\rm sub}=U_1-U_1(\chi_0)-U_1'(\chi_0)(\chi-\chi_0)
-\tfrac12U_1''(\chi_0)(\chi-\chi_0)^2
$$

只是把参考点的值、斜率与曲率吸收入指定有限反项。高阶残项仍受 $m^4$ 放大，下一圈通常需要重新调整。它可作为明确标注的调谐EFT对照，不能称为自然性问题的解决。令 $\mu=m(\chi)$ 却冻结其余运行参数也不能消除反馈；一致RG改进须包含所有参数运行。

保护对称性、粒子谱抵消或其他机制可能改变判断，但本次未构造这些机制。减小 $\xi$ 会同时减小电磁响应，增大 $F_\chi$ 会改变动力学及能量预算，均不能被当作无代价的预测。

## 6. 对本次数值设计的审查

根任务的[固定试验](../experiments/charged_threshold_v01/protocol.md)在 $z=4.2$ 初态比较 $K=1$ 与 $K=e^{2(\chi-\chi_i)}$，使用同一树势和宇宙常数；$m^2=m_i^2g(\chi)$，$g=(1+\xi/\chi^2)/(1+\xi/\chi_i^2)$、$\xi=\chi_i^2$、$\mu=m_i$，只扣初始CW常数，不扣导数。两分支初始 $K=1$，所以同一 $d\chi/d\ln a$ 确实对应相同树级初始动能。

此方案可做条件实验，但需保持以下边界：

1. 零圈强度是关闭量子反馈的树级／探针对照，不能把 $m_i=0$ 后仍称为积分掉重场的阈值模型。
2. 初始常数减法后 $g^2(\ln g-3/2)+3/2$ 可以为负；必须保留Friedmann物理解支、$\chi>0$ 及有效理论域的终止，不能为跑到 $z=0$ 再修改势或宇宙常数。
3. 本次有限圈强度实际对应 $m_i=0.000788768$、$0.0024943$、$0.00788768$ eV。即使 $m_i\gg H$，仍不能当光学测量下的重阈值。这只能是形式低能响应与反馈诊断，不是可行新粒子或光学 $\alpha$ 预测。
4. eV、MeV、GeV阈值可先报告势、力与动能修正的量级；极端刚性下没有必要强行跑ODE，也不能把未算轨迹写成计算结果。

可以报告明确的物质匹配、阈值传递、反馈和实际算出的条件轨迹；不能报告平方指数已被量子理论选择、宇宙历史必然是精确 $1/\ln^2t$，或数学检查通过等于观测支持。

已只读复核[数值实现](../experiments/charged_threshold_v01/code/run_threshold.py)的质量函数、CW减法、动能度量、FLRW方程、能量功积分和阈值响应，未发现与本记录冲突的显式公式。实际轨迹与proper-time独立复算见实验目录，本记录没有重复执行生产模拟。各例都以今天的 $\alpha_0$ 归一化，是共同参考校准；固定UV边界指各例内部的比较，不额外宣称不同模型的数值UV匹配常数相同。

本记录由AI协助定位原始来源、推导和量纲核对，另由协作代理独立核对FLRW交换符号；属于项目内部复核，不是外部同行评审。10项代数恒等式检查全部通过，见[脚本](charged_threshold_action_audit_20260925_checks.py)和[结果](charged_threshold_action_audit_20260925_checks.json)。原始来源版本和文件身份见[来源记录](charged_threshold_action_audit_20260925_sources.json)。
