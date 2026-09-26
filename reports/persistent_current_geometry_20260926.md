# 持久共享电流：完整局部几何、恢复力与带电接触

2026-09-26。结论是：**把模量辅助场接入同一个局部电流，可以产生不随时钟势 \(U\) 消失的横向恢复力；同一作用量也产生持久的带电软质量。** 两者有固定相关系数，不能分别保留或删除。可接受例子需要另一个很小的输入耦合，本计算没有解释其来源。

## Material Passport

- 材料：明确声明的四维局部 Kähler 势与超势；不是新观测、数据拟合或完整重矢量量子匹配。
- 阅读范围：[上一轮说明](../experiments/clock_stabilization_v01/protocol.md)、[上一轮完整四次几何](clock_quartic_geometry_20260926.md)。未导入主实验脚本或旧核验脚本。
- 方法：独立构建完整 \(T,Z\) Kähler 矩阵、辅助场、标量势以及带电二次谱；100 位十进制精度。
- 输出：[独立脚本](persistent_current_geometry_20260926_checks.py)、[全部结果](persistent_current_geometry_20260926_checks.json)。**559/559** 项检查通过。
- 执行：首次质量/电流核验 507/507；随后扩展精确软质量一阶 jet 与四个物理 Hessian，得到 559/559。没有失败检查或容差放宽。放大测试使用独立 \(h_Q\ne h_Z\)，物理例使用 \(h_Q=h_Z\)。
- 研究性质：AI 辅助内部复核，未经外部同行评议。没有新增外部文献断言；下列公式均由声明的局部作用量直接推出。五维嵌入、重场全谱与圈、黎曼相互作用来源及逆对数指数仍未完成。

## 1. 固定一次的参考点与径向延拓

记 \(P=M_p^2\)、\(t=T+\bar T\)、\(v=-i(T-\bar T)\)、\(Z=F_{\rm hol}(\chi+iy)/\sqrt2\)，保留

\[
g=t+\eta(t-1)^2+(t-1)^4+v^4-\frac A{t^2},\qquad
h_Z=1+\frac A{t^3},\qquad W=W_c+w(Z)+M(Z)Q_+Q_-.
\]

星号表示关闭时钟后、调好有限真空能的原参考驻点。令

\[
p=g_t,\quad q=g_{T\bar T},\quad D=p^2-gq,\qquad
G_{T\bar T,*}=\frac{3PD_*}{g_*^2},
\]
\[
X_c=\sqrt{G_{T\bar T,*}}(T-T_*),\quad
k_c=\frac{h_{Z,*}}{g_*}\left[-\frac12(Z-\bar Z)^2\right]
=F_{\rm phys}^2y^2,\quad
S_c=\frac{h_{Q,*}}{g_*}(|Q_+|^2+|Q_-|^2).
\]

本轮局部候选明确为

\[
\boxed{\Omega=g-\frac{h_Zk+h_QS}{3P}
+\frac{g_*}{3P\Lambda^2}
\bigl[k_c+\sigma S_c+\tau|X_c|^2\bigr]^2,\qquad
K=-3P\ln\Omega.}
\]

这里 \(\Lambda=\Lambda_{\rm phys}\)，带星号的归一化和电流中心是固定常数。不能在每个 \(T,\chi\) 点重新定中心，否则计算的是不同作用量。该冻结径向延拓也是新 EFT 输入；它与上一轮使用 \(h_Z(t)\) 作为整个四次项系数的延拓，在涉及 \(T\) 导数的 \(O(AU)\) 项上不同。轴上恢复力及接触值的相关关系相同，某些离轴第一 jet 不同。

零真空能原参考点满足

\[
(t_*-1)t_*^5=-A,\qquad
\eta=-6(t_*-1)^2+\frac{3A}{t_*^4},\qquad q_*=q'_*=0.
\]

有限 \(\rho\) 时，定义物理 \(m_G=|W_c|/(Pg_*^{3/2})\)、\(r_\rho=\rho/(3Pm_G^2)\)，用

\[
q_*=\frac{r_\rho p_*^2}{g_*(1+r_\rho)},\qquad
q'_* =\frac{2p_*q_*}{g_*}-\frac{q_*^2}{p_*}.
\]

新电流平方在 \(X_c=Q=y=0\) 从四阶才开始，因此该参考点的 \(\Omega\) 前三阶 jet、原势值及一阶梯度均不变，**无需为该点再调一次 \(\eta\)**。它在离开 \(T_*\) 后改变轴上模量势，所以不能声称整条多场轴势完全不变。开启时钟后原模型已有的模量力也仍存在。

## 2. 不冻结模量辅助场的完整 Hessian

在参考点，完整两场辅助场为

\[
F^T=\frac{p\bar W}{P\sqrt gD},\qquad
F^Z=-\frac{\bar w_Z}{\sqrt g h_Z},\qquad
F_X=\sqrt{G_{T\bar T,*}}F^T,\quad
F_z=\sqrt{h_{Z,*}/g_*}F^Z.
\]

记 \(f=|F_X|^2\)、\(u=|F_z|^2=U\)、\(C=F_X\bar F_z\)。规范实坐标为 \(X_c=(x+ia)/\sqrt2\)、时钟横向坐标 \(b=F_{\rm phys}y\)。在原几何的同一参考点，新项对势 Hessian 的完整增量是

\[
\boxed{\Delta H=\frac1{\Lambda^2}
\begin{pmatrix}
4\tau^2f+2\tau u&0&-4\tau\operatorname{Im}C\\
0&4\tau^2f+2\tau u&4\tau\operatorname{Re}C\\
-4\tau\operatorname{Im}C&4\tau\operatorname{Re}C&4\tau f+12u
\end{pmatrix}.}
\]

因此

\[
\Delta m_y^2=\frac{12U+4\tau|F_X|^2}{\Lambda^2},\qquad
\Delta m_x^2=\Delta m_a^2=
\frac{4\tau^2|F_X|^2+2\tau U}{\Lambda^2}.
\]

这些是完整逆度规对辅助场的响应，未把 \(T\) 或 \(F^T\) 冻结成无响应的背景常数。\(\tau>0\) 时，新增三维 Hessian 为正；混合项不能漏掉，但不会使这个新增矩阵变成不定。具体地，其非平凡二阶主块行列式分子为

\[
16\tau^3 f^2+40\tau^2fu+24\tau u^2>0
\]

（至少一个辅助场非零）。这不自动证明原 Hessian 与新增项之和、完整滚动扰动或远处势全局稳定。

新项的前三阶 Kähler jet 为零，故参考点的联络不变；势梯度也不变。因此新增**协变势 Hessian**与上式新增坐标 Hessian一致。但完整滚动扰动还涉及场空间曲率、速度与引力约束，本轮没有求解它们。

在 \(U=0,q=0\) 时，\(f=3Pm_G^2\)，所以

\[
\boxed{\Delta m_{y,\rm persist}^2=\frac{12\tau Pm_G^2}{\Lambda^2}.}
\]

对有限 \(\rho\)，这个式子乘 \(p^2/D=1+r_\rho\)。它在时钟关闭后仍存在，是本候选与上一轮纯 \(k_c^2\) 稳定的区别。

## 3. 持久恢复力同时产生持久带电软质量

在 \(Q=0\) 定义 \(J=k_c+\tau|X_c|^2\)。带电度规的分子为

\[
H_Q=h_Q(t)-\frac{2\sigma h_{Q,*}}{\Lambda^2}J,\qquad
\mathcal Z_Q=H_Q/\Omega.
\]

由 \(\ell_Q=\ln H_Q\) 的完整 \(T,Z\) 混合 Hessian 得到

\[
\boxed{\Delta d_Q=\frac{2\sigma}{\Lambda^2}
\bigl(U+\tau|F_X|^2\bigr),\qquad
\Delta d_{Q,\rm persist}=\frac\sigma2\Delta m_{y,\rm persist}^2.}
\]

恢复力和软质量由同一个参数生成。令 \(\sigma=0\) 删除带电接触是额外耦合假设；不能在 \(\sigma\ne0\) 时只保留恢复力。\(\sigma<0\) 可给负软质量，实际带电质量是否为负还需与 \(M^2\) 及 \(B\) 一起判断，不能只看软项符号。

双线性轴值保持 \(\Delta B=0\)，但第一 jet 为

\[
\Delta B_y=-\frac{4i\sqrt2\sigma M F_{\rm phys}F_z}{\Lambda^2}
=\frac{4i\sqrt2\sigma M F_{\rm phys}\sqrt U}{\Lambda^2},
\]
\[
\Delta B_x=\frac{2\sqrt2\sigma\tau M F_X}{\Lambda^2},\qquad
\Delta B_a=-\frac{2i\sqrt2\sigma\tau M F_X}{\Lambda^2}.
\]

最后一个等号采用原有负实 \(w_r\) 约定。新 \(\tau\) 不增加 \(B_y\) 的一阶 jet，但产生模量方向的 \(B\) jet。仅轴值为零不足以删除谱导数及其圈干涉。

对于固定参考点的 \(w\propto e^{-\chi}\)，令

\[
r_w=\frac{F_{\rm phys}\sqrt{U/2}}{Pm_G},\quad
f=3P\frac{p^2}{D}m_G^2|e^{i\theta}-r_w|^2,
\]
\[
f_\chi=6P\frac{p^2}{D}m_G^2r_w(\cos\theta-r_w),\qquad
\boxed{\Delta d_{Q,\chi}=-\frac{4\sigma U}{\Lambda^2}
+\frac{2\sigma\tau}{\Lambda^2}f_\chi.}
\]

在本冻结径向延拓下，写 \(\ell'_0=\partial_t\ln h_Q\)，横向一阶 soft jet 的精确式是

\[
\boxed{\Delta d_{Q,y}=2\sqrt2 F_{\rm hol}\frac{2\sigma}{\Lambda^2}\frac{h_Z}{g}
\left[(1+\tau)\frac{h'_Z}{h_Z}+\ell'_0-\frac qp
+\frac{\tau p}{2g}\right]\operatorname{Im}(F^T\bar F^Z).}
\]

它使用完整辅助场响应

\[
(F^Z)_y=iF^Z-i\sqrt2F_{\rm hol}
\left(\frac{h'_Z}{h_Z}-\frac qp\right)F^T,
\]
\[
(F^T)_y=\frac{iF_{\rm hol}h_Z}{\sqrt2gG_{T\bar T,*}}
\left(2\frac{h'_Z}{h_Z}+\frac pg\right)F^Z.
\]

特别是 \(\tau=0,h_Z=h_Q,q=0\) 时，括号为 \(2h'_Z/h_Z\)，不会相消。近中心领先系数为 \(-24\sqrt2\sigma A F_{\rm phys}\sqrt U m_G\sin\theta/\Lambda^2\)，是上一轮不同径向延拓下该系数的两倍。该差异已由完整谱直接微分核对；旧归档模型不改写。

## 4. 模量四次项、正度规与有限参考域

在 \(y=Q=0\) 上，新几何项为

\[
\Delta g=c_4[(t-t_*)^2+v^2]^2,\qquad
c_4=\frac{g_*\tau^2G_{T\bar T,*}^2}{48P\Lambda^2}.
\]

它包含原来分离的 \((t-1)^4+v^4\) 没有的交叉项 \(2(t-t_*)^2v^2\)，并提高模量质量；不能把它仅当成时钟质量修补。该平方在其固定中心的二阶、三阶几何 jet 仍为零。

局部 EFT 还要求 \(\Omega>0\)、\(K_{T\bar T}>0\)、\(\det K_{I\bar J}>0\) 以及 \(H_Q>0\)。独立计算逐点检查完整度规，不把截断四次式延拓到任意场值。参考点及多个有限邻域点通过这些条件；这不是全场空间无鬼的证明。重场质量、全部规范化辅助场与 EFT 尺度分离也仍须额外检查。

## 5. 有物理归一化的非零 \(A\) 例子

取 \(M_p=2.435\times10^{27}\) eV、\(m_G=10^{-6}\) eV、\(m_{\rm KK}=100\) TeV、\(\Lambda=1\) TeV、\(\sigma=1\)，沿用原有限真空能。由已有有限引力系数得到 \(A=4.3359709873\times10^{-31}\)。定义

\[
r_H=\frac{\tau P}{A\Lambda^2}.
\]

忽略很小的有限真空能项，\(U=0\) 时总时钟质量领先为

\[
m_y^2\simeq12A m_G^2(r_H-2).
\]

令 \(r_H=4\)，即 **额外输入** \(\tau=2.9251519295\times10^{-61}\)，完整局部计算得到：

| 项目 | 数值 |
|---|---:|
| 原负曲率领先项 | \(-1.0406330369\times10^{-41}\) eV² |
| 持久新增曲率 | \(2.0812660739\times10^{-41}\) eV² |
| 完整关闭时钟局部曲率 | \(1.0406330369\times10^{-41}\) eV² |
| 伴随新增 \(d_Q\) | \(1.0406330369\times10^{-41}\) eV² |
| \(\lvert F_X\rvert^2\) | \(1.7787675\times10^{43}\) eV⁴ |

四个物理测试点（\(U=0,U_{\min}\)，相位 \(0,\pi/2\)）的完整三维局部势 Hessian 都为正。\(U=0\) 时参考点是驻点，因此其势 Hessian也等于协变势 Hessian；开启时钟后的结果仅是固定参考点局部切片诊断。

这个例子证明存在声明作用量内的持久局部恢复力，**没有解释为什么 \(\tau\) 应该如此小**。同一组其他参数若取 \(\tau=1\)，得到 \(\Delta m_y^2=7.11507\times10^{19}\) eV²、\(\Delta d_Q=3.557535\times10^{19}\) eV²；它已是完全不同的带电反馈与模量质量问题。把小耦合移到新电流，不能称为消除了小参数问题。

本报告不计算全部带电 CW、重矢量及模量量子修正；持久软质量对慢时钟的影响应由本轮另一个独立谱计算继续检查。尤其 \(U\) 趋近零时，“横向仍有正恢复力”与“辐射反馈仍小于原驱动力”是两个不同条件，前者不保证后者。
