# 全部相位被规范场吸收后的有限辅助场审计：局域端点的刚性分支

日期：2026-09-26。类型：独立解析与数值审计。对象：`higgs_stabilization_v01` 所比较的额外 $U(1)$ 候选；不导入主实验代码，不修改历史实验。

## 结论和范围

在下文明确规定的刚性超对称分支中，加入最后一个规范场并吸收剩余 Higgs 相位后，局域左端时钟流与局域右端超对称破缺流**不能产生所需的静态时钟恢复势**。该结论并非只从零时空动量的接触项消失推断：直接保留有限 $F_X$ 并对重场求驻点，时钟势仍然平坦。

解析结果适用于 $X=Z_a=0$、$W_{\rm clock}=0$、正定 Kähler 度量以及连续存在的重场驻点分支。下面的数值审计只验证这个分支中的重 Higgs 场驻点，未证明超对称破缺后完整的标量 Hessian 正定，未加入完整超引力、非零时钟辅助场、其他超势耦合、量子修正或一般高阶超导数算符。**不能把“这个分支没有恢复势”改写成“所有有限辅助场效应均为零”。**

## 1. Gaussian massive-vector 描述：有限 $F_X$ 的直接移位

取 $n+1$ 个 Higgs pairs、$n+1$ 个 gauge factors。包含 gauge-coupling ratio $\eta>0$ 的实 charge matrix 为

\[
Q_+=\begin{pmatrix}
q&0&\cdots&0&0\\
1&-q&\cdots&0&0\\
\vdots&&\ddots&&\vdots\\
0&\cdots&1&-q&0\\
0&\cdots&0&1&\eta
\end{pmatrix},\qquad H_+=Q_+^TQ_+.
\]

局域 clock source 取左端 Higgs row 的 $1/q$ 倍；hidden source 取完整的右端 Higgs row：

\[
a=(e_0,0),\qquad b=(e_{n-1},\eta),\qquad
 a^TH_+^{-1}b=0.
\]

在规定的 Gaussian Kähler 模型

\[
K=k+|X|^2+\frac12J^TH_+J+
\frac{\sqrt2}{M}\left[(a^TJ)k+(b^TJ)|X|^2\right],
\qquad W=fX,\quad k=y^2,\quad S=|f|^2
\]

中，$X=0$ 时 $K_{X\bar J}=K_{X\bar T}=0$，所以非零 $F_X$ 不会被误当作零。消去 auxiliaries 后，此刚性模型的静态 potential 为

\[
V(j,k)=\frac{S}{1+\sqrt2 b^Tj/M}
+\frac{M^2}{4}\left|H_+j+\frac{\sqrt2}{M}ak\right|^2.
\]

作有限移位

\[
u=j+\frac{\sqrt2}{M}H_+^{-1}ak.
\]

第二项变为 $M^2|H_+u|^2/4$；第一项的 denominator 变为

\[
1+\frac{\sqrt2}{M}b^Tu-
\frac{2}{M^2}b^TH_+^{-1}a\,k
=1+\frac{\sqrt2}{M}b^Tu.
\]

因此整个静态 potential 都与 $k$ 无关。这是有限 $S$ 的代数恒等式，不是 $S\to0$ 或 $p^2\to0$ 的辅助场截断。Heavy stationary branch 的能量及 clock 曲率均保持此性质；kinetic terms 可因场重定义和背景而改变。

## 2. 显式 canonical Higgs 模型也保留此结论

考虑在 chiral coordinates 中给出的局域 norm-difference currents：

\[
\begin{aligned}
K={}&k+|X|^2+\sum_{a=0}^{n}
\bigl(|\Phi_a|^2+|\widetilde\Phi_a|^2+|Z_a|^2\bigr)\\
&+c_I k\bigl(|\Phi_0|^2-|\widetilde\Phi_0|^2\bigr)
+c_X|X|^2\bigl(|\Phi_n|^2-|\widetilde\Phi_n|^2\bigr),\\
W={}&fX+\lambda\sum_{a=0}^{n}Z_a(\Phi_a\widetilde\Phi_a-v^2).
\end{aligned}
\]

带电模平方中的规范场因子已隐含。$c_I,c_X$ 的质量维数为 $-2$；$\eta$ 已包含在 $Q_+$ 中，下面的共同规范耦合为 $g$。对于 $T=(\chi+iy)/\sqrt2$，取 $k=-(T-\bar T)^2/2=y^2$。

定义 $d_a=|\Phi_a|^2-|\widetilde\Phi_a|^2$、$s_a=|\Phi_a|^2+|\widetilde\Phi_a|^2$。在 $X=Z_a=0$ 上，$X$ 度量 block 单独分离，各 $Z_a$ 的度量仍为 1；因而完整的该分支 F potential 为

\[
V_F=\frac{S}{1+c_Xd_n}
+\lambda^2\sum_a|\Phi_a\widetilde\Phi_a-v^2|^2.
\]

这里没有漏掉 inverse-metric mixing：脚本用独立 complex variables 显式检查了 $X$ row 和全部 $Z_a$ rows。其他 superpotential derivatives 在该分支为零。D potential 为

\[
\mu=Q_+^T[d+c_Iks_0e_0],\qquad V_D=\frac{g^2}{2}|\mu|^2.
\]

现在只对左端 pair 作实变换

\[
\Phi_0\mapsto e^\xi\Phi_0,\qquad
\widetilde\Phi_0\mapsto e^{-\xi}\widetilde\Phi_0.
\]

此变换保持左端 product，亦不改变右端 $d_n$，故严格保持 $V_F$。同时

\[
\partial_\xi\mu=2(s_0+c_Ikd_0)Q_{+,0}^T.
\]

只要 $|c_Ik|<1$、左端 pair 非零，则 $s_0+c_Ikd_0>0$。Heavy stationarity 要求

\[
Q_{+,0}\cdot\mu=0.
\]

在同一点对 clock variable 求偏导则得

\[
\partial_kV=c_Is_0g^2Q_{+,0}\cdot\mu=0.
\]

Envelope theorem 随即给出 $dV_{\rm stat}/dk=0$。在 $k=0$ 的 clock mass coefficient 为零；如果 stationary branch 在允许区间连续存在，它的整个静态能量都与 $k$ 无关。

这个论证的关键条件是只有右端 pair 进入 hidden source，以及左端 product-preserving direction 不进入其他 F interactions。若新增相互作用破坏这个条件，必须重新计算。

## 3. Branch 的显式构造与有限 $F_X$ 响应

对实 positive Higgs vevs，令 $p_a=\Phi_a\widetilde\Phi_a>0$。在 $|h|<1$、$h=c_Ik$ 时，左端 D source

\[
r_0=d_0+h\sqrt{d_0^2+4p_0^2}
\]

是从实数轴到实数轴的一一映射，因为它的 derivative 为 $1+hd_0/s_0>0$。其反解为

\[
d_0=\frac{r_0-h\sqrt{r_0^2+4p_0^2(1-h^2)}}{1-h^2}.
\]

取 $r_a=d_a$ for $a>0$，则 energy 完全可以写成

\[
V=\frac{S}{1+c_Xr_n}
+\lambda^2\sum_a(p_a-v^2)^2
+\frac{g^2}{2}r^T(Q_+Q_+^T)r.
\]

新坐标中的表达式不再包含 $k$。Product stationary solution 为 $p_a=v^2$。令 $A=Q_+Q_+^T$、

\[
\gamma=(A^{-1})_{nn}=\eta^{-2},\qquad t=c_Xr_n,
\]

positive branch 由

\[
t(1+t)^2=\frac{Sc_X^2}{g^2\eta^2},\qquad
r=\frac{Sc_X}{g^2(1+t)^2}A^{-1}e_n
\]

给出，其 energy 为

\[
V_{\rm stat}(S)=\frac{S}{1+t}+\frac{St}{2(1+t)^2}.
\]

所以**有限辅助场的自身响应没有消失**：$t\neq0$、真空能偏离 $S$，隐扇区度量也改变；消失的是它对这个局域时钟方向的静态势依赖。若之后规定 $S$ 随另一个背景场变化，这个非平凡的 $V_{\rm stat}(S)$ 仍可产生反馈，必须另行计算。

## 4. 独立数值检查

运行：

```bash
python reports/higgs_finite_auxiliary_20260926_checks.py
```

结果为 **678 / 678 项通过**。这些检查是内部一致性验证，不是 678 个独立物理证据。

- Gaussian finite-$S$ shift 与 source overlap：使用 rational SymPy matrices，$n=1,2,3,4$，$\eta=1/2,2$。
- Canonical Higgs complex-coordinate metric 和 superpotential derivative blocks：显式检查 $n=1$ 的八个 chiral fields；推导对任意 $n\ge1$ 成立。
- 108 个直接 nonlinear stationary solves：$n=1,2,4$，$\eta=0.5,1,2$，$S=0.01,0.5,4$，$k=0,0.01,0.1,0.5$，取 $v=1,g=0.7,\lambda=2,c_I=0.2,c_X=0.3$。每次从零 log amplitudes 开始；未把解析解作为初值。
- 最大 stationarity residual：$5.76\times10^{-15}$。在同一 $n,\eta,S$ 下改变 clock $k$，最大 energy span：$8.88\times10^{-16}$。
- 所有被检查点的 hidden metric 正值，clock/Higgs block 的最小 Schur complement 为 $0.6650>0$。这验证所采样域的 kinetic positivity，不能替代完整 scalar mass Hessian。
- 额外检查了 nonlinear coordinate map 的 45 个正负参数控制点、有限 $F_X$ branch equation 以及 product-preserving derivative。

机器可读结果：[higgs_finite_auxiliary_20260926_checks.json](higgs_finite_auxiliary_20260926_checks.json)。脚本在结果文件中保存 Python/library versions 与自身 SHA-256。初次执行即全部通过，没有放宽阈值或删除失败点。

## 5. 对候选机制的含义

额外 $U(1)$ 能在超对称自由谱中吃掉剩余相位，但这不自动完成所需的持续恢复。对于这里明确的**局域端点、刚性、$W=fX$** 实现，源的物理取法导致静态恢复势缺失，而且有限 $F_X$ 不能挽回这一项。

其他 clock F terms、supergravity curvature、有限 momentum、loop threshold 或额外 interactions 可以改变结果，但必须逐项匹配和检查。不能保留抽象旧 gauge-site source、同时把它未经检查地解释成新模型中同一个局域 Higgs current；两者已非同一 microscopic operator。此审计也不涉及从 Riemann structure 推导真实物理 interaction 或 $1/\ln^2 t$ 驱动。

材料说明：这是项目内部公式的解析/数值审计，无新增外部观测数据，无外部引用事实。由独立 agent 生成和运行，供主研究流程复核；不是外部同行评审。
