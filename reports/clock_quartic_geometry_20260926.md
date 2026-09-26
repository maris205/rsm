# 时钟四次项与共享接触：完整 K/W 的独立局部推导

日期：2026-09-26。结论是：**一个明确加入的四维四次项能增加开时钟时的横向恢复力，并保持实轴势及其一阶导数；但同边界接触、偏离实轴后的势形和有效理论尺度必须一起检查。它不能稳定时钟关闭后的真空。**

## Material Passport

- 材料：本地理论模型的独立推导与可复现数值核验，不是观测数据。
- 输入范围：[上一轮实施说明](../experiments/sequestering_test_v01/protocol.md)、[上一轮完整几何推导](protection_power_correction_20260925.md)。本次没有导入主实验脚本或结果拟合函数。
- 方法：由完整两手征场的 Kähler 矩阵与辅助场计算势和带电谱；100 位十进制精度。
- 输出：[独立脚本](clock_quartic_geometry_20260926_checks.py)、[检查结果及局部谷底](clock_quartic_geometry_20260926_checks.json)。最终 **843/843** 项检查通过。
- 执行记录：第一次四次项检查为 667/667；随后扩展到共享接触后为 843/843。没有失败检查或放宽容差。人为放大的无量纲参数用于检查恒等式，不冒充物理观测点。
- 限定：AI 辅助研究与项目内部独立复核；未经外部同行评议。没有完成五维 UV、完整圈匹配、全场宇宙轨迹或逆对数平方指数的推导。

## 1. 模型、归一化与实轴恒等式

记 $P=M_p^2$、$t=T+\bar T$、$Z=F_{\rm hol}(\chi+iy)/\sqrt2$，

\[
k=-\frac12(Z-\bar Z)^2=F_{\rm hol}^2y^2,\qquad
k_{\rm eff}=k-\frac{k^2}{\Lambda_{\rm hol}^2},
\]
\[
\Omega=g(t,v)-\frac{h(t)k_{\rm eff}+h_Q(t)(|Q_+|^2+|Q_-|^2)}{3P},
\quad K=-3P\ln\Omega,
\]
\[
W=W_c+w(Z)+M(Z)Q_+Q_-,\qquad W_T=0.
\]

涉及具体一阶时钟 jet 时，沿用此前的全纯指数形式 $w=w_r(\chi)e^{-iy}$，$w_r<0$、$\partial_\chi w_r=-w_r$，因此 $U\propto e^{-2\chi}$；$W_c$ 可有一般相位。新增四次质量差本身不依赖这个特定指数形状。

这里 $h=h_Z$；四次项是**新增四维 Wilson 系数的假设**，并非上一轮有限五维引力系数自动导出的结果。参考点定义

\[
F_{\rm phys}^2=\frac{h}{g}F_{\rm hol}^2,\quad
\Lambda_{\rm phys}^2=\frac{h}{g}\Lambda_{\rm hol}^2,\quad
U_{\rm phys}=\frac{|w_Z|^2}{g^2h}.
\]

设 $p=g_t,q=g_{T\bar T},D=p^2-gq$，所有量取 $v=y=Q=0$。完整实轴式为

\[
V=U_{\rm phys}+\frac{3q|W|^2}{Pg^2D},\qquad
F^T=\frac{p\bar W}{P\sqrt gD},\qquad
F^Z=-\frac{\bar w_Z}{\sqrt g h}.
\]

四次项在这条轴上为零，它对 $\Omega$ 的首次改变量是 $O(y^4)$，对场度规的首次改变量是 $O(y^2)$。因此上述轴值、轴上的全部 $\chi$ 导数以及一阶横向力都不变；并未把一般相位下本来非零的横向力删除。

## 2. 新增恢复力及辅助场 jet

逆 Kähler 度规的二阶变化给出

\[
\Delta V_{yy}=\frac{12U_{\rm phys}F_{\rm hol}^2}{\Lambda_{\rm hol}^2},\qquad
\boxed{\Delta\!\left(\frac{V_{yy}}{G_{yy}}\right)
=\frac{12U_{\rm phys}}{\Lambda_{\rm phys}^2}.}
\]

这是固定参考 $T$ 切片上的精确二阶差，适用于有限 $q$、一般 $h,h_Q$ 和 $W_c$ 相位；它没有取 $q=0$ 才成立的限制。这里 $G_{yy}=F_{\rm phys}^2$。新增辅助场 jet 为

\[
\Delta F^Z=\frac{6F_{\rm hol}^2}{\Lambda_{\rm hol}^2}
F^Z_{\rm axis}y^2+O(y^3),\qquad
\Delta F^T=O(y^3).
\]

轴上的度规和一阶度规导数不变，所以 Christoffel 符号也不变；所有轴上势梯度不变。因此该算符对**轴上协变势 Hessian** 的新增 $yy$ 分量也等于 $\Delta V_{yy}$。然而滚动扰动还含场空间曲率项，而场空间 Riemann 张量受到该算符影响；不能仅凭势 Hessian 宣称完成多场宇宙稳定性分析。

最重要的边界是 $U_{\rm phys}=0$：新增质量严格为零。若原模型的 $-24A m_G^2$ 负项继续存在，该四次项**不能解决 clock-off 真空的失稳**。若将来 $U$ 趋于零、$A,m_G$ 不变，有限时段内通过的稳定条件也会再次失效。

## 3. 不共享带电接触时的谱变化

先对上面分离写法定义 $\ell_Q=\ln h_Q$、$c_Q=\partial_t^2\ell_Q$。带电共同软质量与双线性为

\[
d=\frac{2V}{3P}-c_Q|F^T|^2,
\]
\[
B=M_{\rm phys}\left[-\frac{\bar W}{P\Omega^{3/2}}
+F^T\left(\frac{\Omega_T}{\Omega}-2\partial_t\ell_Q\right)
+F^Z\left(\frac{\Omega_Z}{\Omega}-\frac{M_Z}{M}\right)\right].
\]

其中 $M_{\rm phys}=M_{\rm hol}/(\sqrt\Omega h_Q)$，一般离轴时可为复数。四次项不改变 $d,B$ 的轴值和第一横向导数，但会改变第二横向导数：

\[
\boxed{\Delta d_{yy}=\frac{2}{3P}\Delta V_{yy},\qquad
\Delta B_{yy}=\frac{12F_{\rm hol}^2}{\Lambda_{\rm hol}^2}B_{\rm clock,axis},}
\]
\[
B_{\rm clock,axis}=-M_{\rm phys}F^Z\frac{M_Z}{M}.
\]

带电超对称质量平方 $s=|M_{\rm phys}|^2$ 到二阶 $y$ jet 都不变。因此“轴上没有新力”不等于“带电圈不改变横向曲率”；后者必须使用新增的二阶谱 jet。

## 4. 共享电流平方产生的同边界接触

更明确的局部候选采用

\[
\Omega=g-\frac{h k+h_Q S-h[k+\sigma\gamma S]^2/\Lambda_{\rm hol}^2}{3P},
\qquad S=|Q_+|^2+|Q_-|^2,\quad \gamma=\frac{h_{Q,*}}{h_*}.
\]

$\gamma$ 是固定参考比值，不在对 $T$ 求导时重新变化；$\sigma$ 是明确保留的独立接触参数。$Q=0$ 后恢复上一节的四次时钟，但带电度规现在为

\[
\mathcal Z_Q=\frac{h_Q}{\Omega}(1-Ck),\qquad
C(t)=\frac{2\sigma\gamma h}{h_Q\Lambda_{\rm hol}^2},\quad
\ell_Q=\ln h_Q+\ln(1-Ck).
\]

使用完整的 $T,Z$ 混合 Hessian，得到参考轴上的精确新增值

\[
\boxed{\Delta d=C|F^Z|^2=\frac{2\sigma U_{\rm phys}}{\Lambda_{\rm phys}^2},
\qquad \Delta B=0.}
\]

对于既有 $U\propto e^{-2\chi}$，$\Delta d_\chi=-4\sigma U/\Lambda_{\rm phys}^2$，$\Delta B_\chi=0$。一般相位下的第一横向导数则为

\[
\boxed{\Delta d_y=2\sqrt2F_{\rm hol}C
\left(\partial_t\ln h_Q-\frac qp\right)
\operatorname{Im}(F^T\overline{F^Z}).}
\]

推导中关键的完整辅助场恒等式是

\[
(F^Z)_y=iF^Z-i\sqrt2F_{\rm hol}
\left(\partial_t\ln h-\frac qp\right)F^T.
\]

接触的混合 $T\bar Z$ Hessian 与这一辅助场响应一起消去了未受抑制的 $m_G$ 项。不能只微分 $C|F^Z|^2$ 而漏掉混合曲率。这个抵消把最终系数中的 $\partial_t\ln h_Z$ 换成 $\partial_t\ln h_Q$；它并不要求 $h_Z=h_Q$。

在 $w_r<0$ 的既有约定下，双线性第一 jet 是

\[
\boxed{\Delta B_y=-2i\sqrt2CF_{\rm hol}M_{\rm phys}F^Z
=\frac{4i\sqrt2\sigma M_{\rm phys}F_{\rm phys}\sqrt{U_{\rm phys}}}
{\Lambda_{\rm phys}^2}.}
\]

该项没有额外的 $A$ 因子，必须保留其与原 $B$ 的干涉。对 $q=0,h_Q=h_Z=1+A/t^3$ 的领先近中心式，

\[
\Delta d_y\simeq-\frac{12\sqrt2\sigma A F_{\rm phys}\sqrt U\,m_G\sin\theta}
{\Lambda_{\rm phys}^2}.
\]

在固定重整化约定下，只取带电 CW 的线性软质量及双线性干涉，可写为

\[
\Delta V_{1,Q}\supset\frac{s\,\Delta d}{8\pi^2}(L-1),\qquad
\Delta(V_{1,Q})_y\supset
\frac{s\,\Delta d_y}{8\pi^2}(L-1)
+\frac{2\operatorname{Re}(\bar B_0\Delta B_y)}{16\pi^2}L,
\quad L=\ln(s/\mu^2).
\]

这些是指定树级接触的一环选定响应；有限局部势、$d^2$ 高阶项和其他同阶部门仍需分别说明。若该树级接触本身来自一次圈匹配，原始理论中的环阶应再加一。共享平方中的 $S^2$ 同时给出带电自相互作用，不能在完整 UV/圈匹配时假装不存在。该局部平方结构本身不证明某个完整重矢量理论能够实现所有符号与系数。

## 5. 正度规范围

令 $O=\Omega$，下标 $t,y$ 表示实坐标偏导。完整两手征场度规的正定条件可直接检查

\[
O>0,\qquad K_{T\bar T}>0,\qquad
\det K_{I\bar J}>0,
\]
\[
\det K=\frac{9P^2}{2F_{\rm hol}^2O^3}
\left[O(O_{tt}O_{yy}-O_{ty}^2)-O_{tt}O_y^2-O_{yy}O_t^2
+2O_{ty}O_tO_y\right].
\]

在 $q=h'=h''=0$ 控制下，这简化为

\[
\det K=\frac{3Ph p^2}{O^3}
\left(1-\frac{6F_{\rm hol}^2y^2}{\Lambda_{\rm hol}^2}\right).
\]

因此除了 $O>0$，该控制的局部有效范围严格要求

\[
F_{\rm phys}^2y^2<\Lambda_{\rm phys}^2/6.
\]

有限 $q,h',h''$ 时使用完整行列式，不把简化边界误称精确结果。共享接触还要求 $1-Ck>0$；在参考归一化及 $\sigma=\pm1$ 控制中，上述时钟范围内这个带电条件也满足。数值检查覆盖正定范围内的多个离轴点，以及简化控制中越过边界后负行列式的例子。

## 6. 一般相位的谷底与尺度限制

写局部势

\[
V(\chi,y)=V_0(\chi)+J(\chi)y+\tfrac12H(\chi)y^2+\cdots.
\]

在 $H>0$ 且位移足够小时，$y_*\simeq-J/H$，$\Delta V_{\rm valley}\simeq-J^2/(2H)$。如果四次稳定贡献占主导，且仅取

\[
J\simeq-6\sqrt2A F_{\rm phys}\sqrt U\,m_G\sin\theta,\qquad
H\simeq12F_{\rm phys}^2U/\Lambda_{\rm phys}^2,
\]

那么领先谷底校正恰为

\[
\Delta V_{\rm valley}\simeq-3A^2m_G^2\Lambda_{\rm phys}^2\sin^2\theta,
\]

它在固定参数下与 $\chi$ 无关。但已有横向曲率、有限真空能和非线性离轴项使这种常数结论不再精确。独立脚本直接求出三个一般相位的真实局部 $y$ 驻点：放大测试中，精确谷底校正相对于该领先常数分别偏离约 **0.070%、0.239%、0.723%**。所有点具有正度规和正局部横向曲率；这些例子检验近似的限制，不是新宇宙轨迹。

尺度也不能略过。若用这项抵消 $-24A m_G^2$，需要

\[
\Lambda_{\rm phys}^2<\frac{U_{\min}}{2A m_G^2}.
\]

对上一轮 $m_G=1$ eV、$m_{\rm KK}=1$ TeV、$U_{\min}\simeq8.33\times10^{-16}$ eV$^4$ 的示例，右侧意味着 $\Lambda_{\rm phys}\lesssim3.1$ GeV。如果同时把 $\Lambda$ 当作有效理论截断尺度，它低于原有 100 GeV 带电质量，因此不能沿用原带电圈计算而宣布这一示例已获救。较小 $m_G$ 可改变这个层次判断，需要逐点同时检查；而 $\Lambda$ 是否确为物理截止，还需具体匹配模型确定。

可保留的结果是一个明确的条件机制及其相关谱响应：它增加开时钟的恢复力，有限区间内可能保持慢变化，且暴露了同边界接触和层次要求。它尚未给出渐近稳定真空、完整高能来源、真实物理相互作用或 $1/\ln^2t$ 的指数来源。
