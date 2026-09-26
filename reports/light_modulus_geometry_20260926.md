# 保留剩余 Higgs 手征模：准确商空间动能与端点接入

日期：2026-09-26。基线提交 `fbfae9851ce0c38452d18d2ec8d70b1d840e6c6d`。本报告是 [light_modulus_v01](../experiments/light_modulus_v01/README.md) 的独立几何核验，不修改前轮归档和正式正文。

**原链保留一个完整的轻手征模，而非仅一个相位。它在左端的线性接入极小，在右端的接入为正常大小。准确的非线性动能可以求出；把该模固定在原点后只保留重矢量接触，不能构成完整低能理论。**

## Material Passport

- 问题：不添加新规范群，也不假定剩余模已经变重时，原 Higgs 链的低能标量几何和局部流是什么？
- 输入：[原场相位算符审计](higgs_phase_operators_20260926.md)、[局部物理流匹配](higgs_local_currents_20260926.md)。输入哈希和代码哈希记录在机器结果中。
- 假设：刚性超对称，规范化的原始 Kähler 势，固定乘积约束，零 FI 项，原整数电荷，树级两导数展开。乘积径向模式和重矢量采用领先超对称消元。
- 方法：解析约束求解、直接分量动能、相位速度的规范投影、精确有理矩阵恒等式、100／180 位计算。没有导入主实验实现。
- 输出：[独立脚本](light_modulus_geometry_20260926_checks.py)、[507项机器核验](light_modulus_geometry_20260926_checks.json)、[来源阅读范围](light_modulus_geometry_20260926_sources.json)。另一个代理独立审计动能及源项的因子和全局坐标范围。
- 证据性质：内部理论和实现核验，不是观测证据，也不是外部同行评议。没有统计检验或数据拟合。

重手征场的超势驻点和重矢量的 Kähler 驻点，是领先超对称有效作用的标准消元条件；辅助场同样要相对于重质量保持足够小。因此以下“准确”指**声明的约束商空间内准确**，不表示有限辅助场及有限径向质量修正全部为零。[Brizi、Gomez-Reino、Scrucca，arXiv:0904.0370v2，摘要](https://arxiv.org/abs/0904.0370v2)

## 1. 原电荷矩阵留下一个复模量

令 $a=0,\ldots,n$ 标记 $n+1$ 对相反电荷的手征场，原 $U(1)^n$ 电荷矩阵为

\[
Q_0=q e_0^T,\quad
Q_a=e_{a-1}^T-qe_a^T\quad(1\le a<n),\quad
Q_n=e_{n-1}^T.
\]

超势为

\[
W_{\rm link}=\lambda\sum_a Z_a(\Phi_a\widetilde\Phi_a-v^2).
\]

在产品约束成立、$Z_a=0$ 的支上，定义

\[
w=(1,-q,-q^2,\ldots,-q^n)^T,\qquad
\mathcal W=w^Tw=\sum_{a=0}^{n}q^{2a},\qquad u=w/\sqrt{\mathcal W}.
\]

注意从第二项开始均为负号，不是交替符号。$Q^Tu=0$，$u^Tu=1$，并且 $Q$ 满列秩。共 $n+1$ 个乘积约束和 $n$ 个复化规范商，给出

\[
2(n+1)-(n+1)-n=1
\]

个复维数。这包括一个实振幅方向和一个角方向。

## 2. 全非线性 D 平坦支及动能

定义有质量平方量纲的实变量 $z$，以及

\[
d_a=|\Phi_a|^2-|\widetilde\Phi_a|^2,\qquad
s_a=|\Phi_a|^2+|\widetilde\Phi_a|^2.
\]

$Q^Td=0$ 和 $|\Phi_a\widetilde\Phi_a|^2=v^4$ 共同给出

\[
\boxed{d_a=u_a z,\qquad
s_a(z)=\sqrt{4v^4+u_a^2z^2}.}
\tag{1}
\]

原点两场都等于 $v$。整条约束支可写为

\[
\Phi_a=v e^{r_a+i\varphi_a},\qquad
\widetilde\Phi_a=v e^{-r_a-i\varphi_a},\qquad
r_a=\frac12\operatorname{arsinh}\frac{u_a z}{2v^2}.
\]

令 $\theta=u^T\varphi$、$h(z)=\sum_a u_a^2/s_a(z)$。规范投影后的准确两导数标量动能为

\[
\boxed{\mathcal L_{\rm kin}=
\frac{h(z)}4(\partial z)^2+\frac1{h(z)}(\partial\theta)^2.}
\tag{2}
\]

推导中的因子可以直接由原场验证。单个场对的径向动能是 $s_a(\partial r_a)^2$，且 $dr_a/dz=u_a/(2s_a)$，所以径向系数为 $h/4$。相位速度 $\nu=\partial\varphi-QA$ 满足 $u^T\nu=\partial\theta$。在该约束下最小化 $\sum_a s_a\nu_a^2$，得到 $\nu_a=u_a\partial\theta/(s_a h)$，角向系数因此是 $1/h$。对所有有限实 $z$，两个系数都为正。

局部规范不变的全纯坐标取为

\[
\boxed{C=\sqrt2v\sum_a u_a\ln(\Phi_a/v).}
\tag{3}
\]

设 $R=\sum_a u_a r_a$，则 $C=\sqrt2v(R+i\theta)$，且

\[
R'(z)=\frac{h(z)}2,\qquad
K_{\rm quot}(z)=\sum_a s_a(z),\qquad
\boxed{K_{C\bar C}=\frac1{2v^2h(z)}}.
\tag{4}
\]

在原点，$h(0)=1/(2v^2)$，故 $K_{C\bar C}(0)=1$。若 $C=(\sigma+i a)/\sqrt2$，则

\[
\sigma=v\sum_a u_a\operatorname{arsinh}\frac{u_a z}{2v^2}
=\frac{z}{2v}+O(z^3/v^5),\qquad a=2v\theta.
\tag{5}
\]

**$\sigma=z/(2v)$ 只是原点的线性近似。** 全局沿实支规范化的弧长变量应为 $\sigma_{\rm can}(z)=\int_0^z\sqrt{h(z')/2}\,dz'$。它通常不同于全纯坐标的实分量；不能同时把全局规范化和全纯性当作免费条件。

式(3)是对数坐标，在全局覆盖空间上定义。因为 $w$ 是本原整数核，规范不变角具有 $\theta\sim\theta+2\pi/\sqrt{\mathcal W}$ 的周期。原点规范化角场 $a$ 的周期为 $2\pi f_a$，$f_a=2v/\sqrt{\mathcal W}$。该周期本身不等于已经证明的散射截断尺度，尤其当前并没有相位势。

## 3. 局部流同时接到重模式和保留的轻模式

保留明确的局部非最小耦合

\[
K_{\rm src}=c_I I d_0+c_X S_Xd_n,\qquad S_X=|X|^2,
\]

其中 $I,S_X$ 质量维数为2，$c_I,c_X$ 质量维数为$-2$。商空间的领先源项是

\[
K_{\rm eff}\supset z(c_Iu_0 I+c_Xu_n S_X).
\tag{6}
\]

在原点附近，它给出

\[
K_{\rm eff}\supset(g_I I+g_X S_X)\sigma+O(\sigma^3),\qquad
\boxed{g_I=2v c_Iu_0,\quad g_X=2v c_Xu_n.}
\tag{7}
\]

因此重矢量接触很小，并不意味着源对整个低能谱的接入都很小。$q=3,n=127$ 时，

\[
u_0=2.39896763795084\times10^{-61},\qquad
u_n=-0.942809041582063.
\]

相同源系数时，$g_I/g_X=-3^{-127}\simeq-2.54449\times10^{-61}$。若沿用复现旧站点归一化的 $c_X=3c_I$，该比值是 $-8.48163\times10^{-62}$。**右端 SUSY 破缺来源可以对这个剩余实方向产生正常大小的作用；不能因左端接入小就忽略它。**

重、轻投影分别为

\[
P=Q(Q^TQ)^{-1}Q^T=\mathbf1-uu^T,\qquad L=uu^T.
\]

| 端点 | 保留轻模的平方投影 $L_{aa}$ | 重子空间的平方投影 $P_{aa}$ |
|---|---:|---:|
| 左端 | $5.75505\times10^{-122}$ | $1-5.75505\times10^{-122}$ |
| 右端 | 约 $8/9$ | 约 $1/9$ |

右端轻／重平方投影之比约为8。两端之间 $P_{0n}=-u_0u_n$、$L_{0n}=u_0u_n$。这是完整模空间的投影恒等式，**不是说具有不同传播分母的轻、重物理交换在任意动量自动抵消**。

## 4. 有限背景上的重场接触也能保留

为了将轻模与重接触放在同一个作用中，令 $S=\operatorname{diag}(s_a)$，用 $b$ 表示源列向量，$b_0=c_I I$、$b_n=c_X S_X$，其余为零。固定全纯 $C$，消去原重规范方向给出

\[
\boxed{K_{\rm eff}=\sum_a s_a+z\,u^Tb
-\frac12 b^T\left(S-\frac{uu^T}{h}\right)b+O(b^3).}
\tag{8}
\]

其中

\[
S Q(Q^TSQ)^{-1}Q^TS=S-\frac{uu^T}{h}.
\]

这个半正定核是固定轻坐标后重方向的响应，零向量为 $S^{-1}u$。在原点，二次项恢复 $-v^2 b^T(\mathbf1-uu^T)b$。相应的归一化端点交叉系数为

\[
\tau_{\rm loc}=\frac{c_X}{c_I}\frac{-u_0u_n}{1-u_0^2}.
\]

对相同 $c$，它为 $2.26176837952281\times10^{-61}$；取 $c_X=3c_I$ 后为旧站点值 $6.78530513856844\times10^{-61}$。这些重接触仍然存在，但式(6)的轻模接入必须同时留下。

式(8)是所声明微观源的诱导项。若微观 Kähler 势另有显式 $IS_X$ 或其他二次源算符，必须额外加入，规范对称性并不自动禁止它们。

## 5. 相位保护及本报告不能证明的事项

现有产品超势和所有模平方源项都在 $\Phi_a\to e^{iu_a\alpha}\Phi_a$、$\widetilde\Phi_a\to e^{-iu_a\alpha}\widetilde\Phi_a$ 下不变。因此角方向在所声明刚性经典模型中没有势，不能把保留的轻手征模擅自换成一个已经有质量的实标量。

原矢量型 Higgs 对对残余全局流的混合规范三角异常也逐对抵消：

\[
u_a Q_{ai}Q_{aj}+(-u_a)(-Q_{ai})(-Q_{aj})=0.
\]

引力和纯全局三角异常同样相消。这限于列出的场内容，不排除新增物质、非微扰作用或量子引力效应。特别是它不能替代一个已计算的紫外保护机制。

有限 $F_X$ 可能使重径向模式和重实规范标量偏离本报告的约束支；完整势的稳定性要结合这些响应继续判断。这里没有证明固定乘积、D平坦商与全有限辅助场极小化完全等价，也没有导入未声明的稳定质量。

复现命令：`python reports/light_modulus_geometry_20260926_checks.py`。本轮 **507/507** 项一次通过，包括原场动能的直接微分、一般权重矩阵消元、规范速度投影、相位不变性和高精度端点数值；没有初始失败需要覆盖。一个独立代理另行确认了式(2)、式(4)、式(8)的因子，并强调固定全纯坐标及全局对数周期的范围。

**这一步建立了可继续计算的轻模低能作用，尚未建立稳定宇宙学解。没有加入或导出 $1/\ln^2 t$ 驱动，也没有从黎曼结构导出真实物理相互作用。**
