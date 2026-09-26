# 保留轻手征模量后的局部 EFT：小交叉项与未受抑制的源驱动

日期：2026-09-26。对应 `light_modulus_v01` 的独立解析审计。所有公式均为本轮对明确作用的推导；没有导入主实验实现，也没有调用观测数据。本文不修改正式稿或历史归档。

**结果：只消去原链的重矢量时，旧的极小端点交叉接触确实能够保留；但剩余轻模量与右端破缺源的直接耦合并不小。原点已经不是驻点，因而原来把轻模量固定后读出的“时钟质量”不能充当完整真空的稳定性判据。**

## Material Passport 与适用范围

- 输入：原 $U(1)^n$ 电荷链、canonical Higgs pairs、乘积超势约束，以及左右两端的非最小 norm-difference currents。
- 近似：树级、刚性超对称、重矢量的领先导数匹配；径向乘积模式处于其约束面。保留全部轻手征场 $C$，不冻结其实部或虚部。
- 源：$W=fX$、$S=|f|^2$；本报告没有为 $f$ 添加宇宙时间依赖。若引用既有超引力标尺 $S=3P m_{3/2}^2$，这里只把它用作一个固定数值输入。
- 未做：完整超引力驻点、一般允许算符、量子修正、全宇宙演化、真实观测拟合，以及黎曼结构到 $1/\ln^2t$ 驱动的推导。
- 检验：[独立代码](light_modulus_response_20260926_checks.py)、[机器结果](light_modulus_response_20260926_checks.json)。110 项检查涵盖精确有理矩阵、局部指标和独立能量积分；不是独立物理证据的计数。

## 1. 轻方向与规范无关的重方向投影

令 $a=0,\ldots,n$ 标记 Higgs pairs，电荷矩阵各行为

\[
Q_0=q e_0^T,\qquad Q_a=e_{a-1}^T-qe_a^T\ (1\leq a<n),
\qquad Q_n=e_{n-1}^T.
\]

注意，左零向量的符号为

\[
w=(1,-q,-q^2,\ldots,-q^n)^T,\qquad
\mathcal W=\sum_{a=0}^nq^{2a},\qquad u=w/\sqrt{\mathcal W}.
\]

它不是交替正负号。$Q^Tu=0$，而消去重矢量产生的投影为

\[
P=Q(Q^TQ)^{-1}Q^T=\mathbf1-uu^T.
\]

写 $\Phi_a=v e^{\xi_a}$、$\widetilde\Phi_a=v e^{-\xi_a}$，并在原点定义规范不变的轻手征坐标

\[
C=\sqrt2v\,u^T\xi,\qquad \sigma=\sqrt2\operatorname{Re}C.
\]

于是 $C$ 在原点 canonical，实部的动能为 $(\partial\sigma)^2/2$。虚部同样必须保留。令 $L_a=\xi_a+\bar\xi_a+2g(QV)_a$，则 $u^TL=\sigma/v$。

端点作用取

\[
K=K_{\rm light}+2v^2\sum_a\cosh L_a
+2v^2\sum_aJ_a\sinh L_a,
\qquad J=c_I k e_0+c_X S_Xe_n,
\]

其中 $k=y^2=-(Z-\bar Z)^2/2$、$S_X=|X|^2$，$c_I,c_X$ 的质量维数为 $-2$。

二次 $L$ 块为 $v^2L^TL+2v^2J^TL$。仅沿重方向极小化，得到

\[
L=u\sigma/v-PJ,
\]
\[
\boxed{K_{\rm eff}=K_{\rm light}+\sigma^2
+2v\sigma\,u^TJ-v^2J^TPJ+\cdots.}
\tag{1}
\]

线性轻流项与重接触项来自同一作用，不能只保留后者。省略号包含 $\sigma$ 的高阶几何项和更高场次项；这些项不改变下面原点处明确给出的指标。

## 2. 旧接触可以匹配，但轻模量隐藏源耦合不小

用

\[
v^2c_I^2P_{00}=\Lambda^{-2},\qquad c_X=q c_I
\]

归一化后，式 (1) 含

\[
\Delta K_{\rm heavy}=-\Lambda^{-2}(k^2+2\tau kS_X+S_X^2),
\qquad
\tau=\frac{q-q^{-1}}{q^n-q^{-n}}.
\tag{2}
\]

同时必须写出

\[
K_{\rm light\,currents}=\beta_I\sigma k+\beta_X\sigma S_X,
\quad
\beta_I=\frac{2u_0}{\Lambda\sqrt{P_{00}}},\quad
\beta_X=\frac{2q u_n}{\Lambda\sqrt{P_{00}}}.
\tag{3}
\]

两条有用的精确恒等式为

\[
\boxed{\beta_I\beta_X=-4\tau/\Lambda^2},\qquad
\boxed{\beta_X^2\Lambda^2=
\frac{4(q^2-1)}{1-q^{-2n}}.}
\tag{4}
\]

因此 $n$ 增大可以使左端流与交叉系数极小，却不能使右端源对轻模量的驱动同样变小。它们是不同的物理量。

采用共同物理矢量尺度 $M=2gv$，并沿用 $M=\Lambda\sqrt{(Q^TQ)^{-1}_{00}}$ 的旧匹配，有

\[
v=\frac{\Lambda\sqrt{P_{00}}}{2gq},\qquad
c_I=\frac{2gq}{\Lambda^2P_{00}},\qquad c_X=q c_I.
\]

在 $q=3,n=127,\Lambda=2\ \mathrm{TeV},g=1$ 上，$v\simeq333.333\ \mathrm{GeV}$，

\[
\begin{aligned}
u_0&=2.39896764\times10^{-61},&u_n&=-0.942809042,\\
\tau&=6.78530514\times10^{-61},&\beta_X\Lambda&=-5.65685425.
\end{aligned}
\]

这项未受链长抑制的耦合是保留轻自由度后的首要结果。

## 3. 原点不是驻点：区分坐标曲率、协变曲率和真空质量

在 $X=0$、局部原点附近，仅凭式 (1) 的这些项就有

\[
V=\frac{S}{1+\beta_X\sigma-2\tau y^2/\Lambda^2}+\cdots,
\qquad
\left.\partial_\sigma V\right|_0=-S\beta_X\ne0.
\tag{5}
\]

所以原点不是完整轻场理论的真空。固定 $C$ 后计算确实给出

\[
\left.\partial_y^2V\right|_0=\frac{4\tau S}{\Lambda^2},
\tag{6}
\]

但这只是非驻点处的普通坐标 Hessian。例如作局部场重定义

\[
\sigma'=\sigma-\frac{2\tau}{\beta_X\Lambda^2}y^2,
\]

式 (5) 的分母在保留阶数变为 $1+\beta_X\sigma'$，因此在固定 $\sigma'$ 下的普通 $yy$ Hessian 为零。这个变换也改变动能，不能据此宣称物理作用消失。

原坐标中，Kähler 动能在原点附近给出

\[
G_{\sigma\sigma}=1,\qquad
G_{yy}=1+\beta_I\sigma+\cdots,\qquad
G_{\sigma y}=\beta_I y+\cdots,
\qquad \Gamma^{\sigma}_{yy}=\beta_I/2.
\]

协变 Hessian 因而是

\[
\left.\nabla_y\nabla_yV\right|_0
=\frac{4\tau S}{\Lambda^2}+\frac{S\beta_I\beta_X}{2}
=\boxed{\frac{2\tau S}{\Lambda^2}}.
\tag{7}
\]

它对上述坐标变换保持不变，也仍然不能直接称为真空质量：背景正在加速离开该点。真实滚动解的线性扰动还需要完整动能、轨迹及背景条件。

一个解释上的交叉检查是：$\chi$ 仍是势的精确平移方向，但在这个非驻点有 $\Gamma^\sigma_{\chi\chi}=-\beta_I/2$，所以 $\nabla_\chi\nabla_\chi V=2\tau S/\Lambda^2$ 也不为零。它没有破坏平移对称性，更不能被解释为该方向获得了真空质量。普通曲率和协变曲率在这里都需要随背景共同解释。

类似地，把 $C$ 从轻谱中删掉会漏掉 $X$ 的度量混合。在原点附近

\[
K_{X\bar X}=1+\beta_X\sigma-4|X|^2/\Lambda^2,
\qquad K_{C\bar X}=\beta_X X/\sqrt2.
\]

inverse metric 中的 Schur complement 使局部 $X$ Hessian 变为

\[
\boxed{m^2_{X,\,\text{local}}=
S\left(4/\Lambda^2+\beta_X^2/2\right).}
\tag{8}
\]

在当前 $q=3,n=127$ 上，它约为 $20S/\Lambda^2$；原来固定 $C$ 的值为 $4S/\Lambda^2$。两者不能混用。式 (8) 仍是非驻点处的局部 Hessian 指标，而非已建立真空的粒子质量。

取数值源 $S=3(5.929225\times10^{54}\ \mathrm{eV}^2)(10^{-15}\ \mathrm{eV})^2=1.7787675\times10^{25}\ \mathrm{eV}^4$，则

| 指标 | 数值 |
|---|---:|
| $\partial_\sigma V$ | $5.03111425\times10^{13}\ \mathrm{eV}^3$ |
| 普通 $yy$ Hessian | $1.20694803\times10^{-59}\ \mathrm{eV}^2$ |
| 协变 $yy$ Hessian | $6.03474013\times10^{-60}\ \mathrm{eV}^2$ |
| 局部 $X$ Hessian | $88.938375\ \mathrm{eV}^2$ |
| 上一行平方根，仅作局部尺度 | $9.43071445\ \mathrm{eV}$ |

当前相互作用仅依赖模平方与 Higgs pair products，保留沿 $u$ 的全局相位变换。因此轻手征场的相位方向在这个经典作用中仍然无势；没有因为实方向被源驱动就同时获得相位质量。

## 4. 固定源平直时空中的局部瞬态：独立能量积分

这是对主实验短程示例的独立核对，范围是**固定 $S$、均匀 Minkowski 背景、沿零源 D-flat 面、从静止出发**。不是宇宙年龄预测，也不是完整超引力寿命。

沿 $d_a=|\Phi_a|^2-|\widetilde\Phi_a|^2=z u_a$，令 $z=-2v^2r$，定义

\[
B(r)=\sum_a\frac{u_a^2}{\sqrt{1+u_a^2r^2}},\qquad
b=-2v^2c_Xu_n>0.
\]

对实正 Higgs vevs，直接从两个 canonical kinetic terms 推出

\[
\mathcal L_{\rm kin}=\frac{v^2}{2}B(r)(\partial r)^2,
\qquad V(r)=\frac{S}{1+br}.
\tag{9}
\]

以 $s=\sqrt S\,t/v$ 为自然单位时间，初始 $r=r'=0$ 的能量守恒为

\[
\frac12B(r)(r')^2+\frac1{1+br}=1.
\]

所以达到指定 $r_*$ 的时间可以不解 ODE，直接由

\[
s_* =\int_0^{r_*}\sqrt{\frac{B(r)(1+br)}{2br}}\,dr
\tag{10}
\]

计算。代码以 $r=t^2$ 去除原点可积奇性，使用 100 位精度。

选定 $br_*=0.1$，表示隐藏度量变为 $1.1$，**不等于观测到10%的物理常数漂移**。本基准结果为

\[
\begin{aligned}
b&=0.9428090415820634,\\
r_*&=0.1060660171779821,\\
s_*&=0.4819126356192804,\\
\hbar v/\sqrt S&=5.202174544832927\times10^{-17}\ \mathrm{s},\\
t_*&=2.506993645851966\times10^{-17}\ \mathrm{s}.
\end{aligned}
\]

这里给出的数字只说明指定局部截断中的源可以很快驱动轻实方向，不能把它直接延伸为遥远宇宙的演化结论。$S/v^4=1.440801675\times10^{-21}$ 是重方向有限源响应的一个小参数；确认整个大场轨迹有效还需要另行匹配，本文没有跨越该范围。

## 5. 对下一步的限制

小的端点交叉并没有保护全部低能自由度。进一步提出稳定机制时，至少需要同时计算轻实方向的驻点、剩余相位、完整轻场动能，以及端点交叉是否保留。不能继续只把式 (2) 当作完整 EFT。

本文没有新增用于达成目标结果的软质量或势函数；若后续引入，必须明确为额外作用假设并重新核查它与原动力学来源的联系。$1/\ln^2t$ 驱动以及黎曼结构到真实相互作用的推导仍然是待完成的问题。

复现：`python reports/light_modulus_response_20260926_checks.py`。初次运行 **110/110** 项通过，未放宽容差或删除失败点。脚本保存版本、输入、精度及自身 SHA-256。另一个同模型家族的内部代理逐项重新推导了实度规、Christoffel 系数、隐藏 Schur complement 和坐标变换；它共享当前上下文，并非盲审。全部核查属于内部计算复核，不属于外部同行评审。
