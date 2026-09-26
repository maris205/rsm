# 原始 Higgs 链保留轻模量时的有限辅助场驻点审计

日期：2026-09-26。独立解析与数值核验；不修改历史实验、正文或模型假设。

**在本报告明确限定的刚性分支中，只要右端超对称破缺源非零、右端耦合非零，原始 rank-$n$ Higgs 链就没有有限的内部驻点。** 把最后一个模量留在低能理论后，重场可以在每个固定模量位置达到唯一的约束极小值，但沿剩余实模量方向的势严格单调。因此，不能把“固定该模量时得到的重场恢复项”直接解释为整个模型已有稳定真空。

这一结果仅针对下面的局域 Kähler 耦合、$X=Z_a=0$、时钟辅助场为零的全局超对称树级分支。它不是对超引力、其他超势、非零 $X,Z_a$ 分支或所有稳定机制的否定。

## 材料与可复现范围

- 输入：[上一阶段 canonical Higgs 模型与有限辅助场约定](higgs_finite_auxiliary_20260926.md)、[电荷链和质量约定](higgs_extra_gauge_20260926.md)。本次恢复原始 $n$ 个规范因子，不加入第 $n+1$ 个规范因子。
- 无新增外部观测数据；全部结论由显式给定的模型推导。没有把黎曼结构、$1/\ln^2t$ 驱动或真实宇宙时间映射作为已经证明的输入。
- 独立脚本：[light_modulus_stationarity_20260926_checks.py](light_modulus_stationarity_20260926_checks.py)。只使用 NumPy、SciPy 和 SymPy，未导入主实验实现。
- 最终机器结果：[checks.json](light_modulus_stationarity_20260926_checks.json)；首次诊断和脚本保留为 `light_modulus_stationarity_20260926_initial_checks.{json,py}`。JSON 记录了输入和脚本的 SHA-256、软件版本以及所有失败与通过项目。
- 本报告由项目内部独立 agent 推导、编写和执行，属于内部可复现审计，不是外部同行评审；检查数不是独立物理证据数。

## 1. 模型、坐标与正定区域

取 $n\ge1$、$q$ 为正整数，原始电荷矩阵为

\[
Q=\begin{pmatrix}
q&0&\cdots&0\\
1&-q&\cdots&0\\
0&1&\ddots&0\\
\vdots&&\ddots&-q\\
0&\cdots&0&1
\end{pmatrix}_{(n+1)\times n}.
\]

它的 primitive left-null vector 是

\[
w=(1,-q,-q^2,\ldots,-q^n)^T,\qquad
W=w^Tw=\sum_{a=0}^{n}q^{2a},\qquad u=w/\sqrt W.
\]

这里从第一个内部行开始都取负号；不能误写成交替符号 $(-q)^a$。有 $Q^Tu=0$、$u_n\ne0$。

定义 $k=y^2=-(T-\bar T)^2/2$，在带电模平方内隐含规范场因子。所检验的模型为

\[
\begin{aligned}
K={}&k+|X|^2+\sum_{a=0}^{n}
(|\Phi_a|^2+|\widetilde\Phi_a|^2+|Z_a|^2)\\
&+c_Ik(|\Phi_0|^2-|\widetilde\Phi_0|^2)
+c_X|X|^2(|\Phi_n|^2-|\widetilde\Phi_n|^2),\\
W={}&fX+\lambda\sum_{a=0}^{n}Z_a(\Phi_a\widetilde\Phi_a-v^2),
\qquad S=|f|^2>0.
\end{aligned}
\]

$c_I,c_X$ 的质量维数均为 $-2$。在 $X=Z_a=0$、$W_{\rm clock}=0$ 的分支上，$W_{\Phi_a}=W_{\widetilde\Phi_a}=0$，而 $W_X=f$ 和 $W_{Z_a}=\lambda(p_a-v^2)$ 均予以保留。$X$ 度量块以及全部 $Z_a$ 度量块分离，故

\[
V=\frac{S}{1+c_Xd_n}
+\lambda^2\sum_a|p_a-v^2|^2
+\frac{g^2}{2}|Q^Tr|^2,
\]

\[
d_a=|\Phi_a|^2-|\widetilde\Phi_a|^2,
\quad s_a=|\Phi_a|^2+|\widetilde\Phi_a|^2,
\quad p_a=\Phi_a\widetilde\Phi_a,
\quad r=d+c_Iks_0e_0.
\]

要求全部 Higgs pair 非零、隐藏度量 $1+c_Xd_n>0$，并处于完整 Kähler 度量正定的开集内。必要条件包括 $|c_Ik|<1$，但它本身不是完整正定性的充分条件。

在固定非零 product $p_a$ 的条件下，$d_a$ 可通过 pair 的实相反缩放局部任意变化。特别地，记 $h_I=c_Ik$，则

\[
r_0=d_0+h_I\sqrt{d_0^2+4|p_0|^2},\qquad
\frac{\partial r_0}{\partial d_0}=1+h_I\frac{d_0}{s_0}>0.
\]

所以在 $|h_I|<1$ 时，$d_0\leftrightarrow r_0$ 是实轴上的可逆映射；反解为

\[
d_0=\frac{r_0-h_I\sqrt{r_0^2+4|p_0|^2(1-h_I^2)}}{1-h_I^2}.
\]

这一步没有取 $F_X\to0$，没有假设径向乘积已经固定，也没有假设 D 项为零。

## 2. 一行驻点障碍：沿剩余模量消去全部 D 力

保持全部 product $p_a$ 不变，作允许的 infinitesimal variation

\[
\delta r=\varepsilon u.
\]

因为 $Q^Tu=0$，D potential 完全不变；乘积 F potential 也完全不变。对于 $n\ge1$，$r_n=d_n$，因此只剩隐藏辅助场的响应：

\[
\boxed{
\left.\frac{dV(r+\varepsilon u,p)}{d\varepsilon}\right|_{\varepsilon=0}
=-\frac{Sc_Xu_n}{(1+c_Xd_n)^2}\ne0.
}
\]

只要 $S>0$、$c_X\ne0$、全部 pair 非零且处于正定区域内部，这个方向总是存在。驻点要求所有允许方向的一阶导数为零，与上式矛盾。此结论对任意有限 $\lambda$、任意有限 $F_X$ 成立；径向位移无法消掉这个 product-preserving force。

在 $k=0$ 也可直接写成

\[
\nabla_dV=-\frac{Sc_X}{(1+c_Xd_n)^2}e_n+g^2QQ^Td
+\text{固定 product 方向上为零的项},
\]

再与 $u$ 点乘即可得到相同障碍。有限 $k$ 只是把有效坐标从 $d_0$ 换成 $r_0$，并不产生一个新驻点。

该证明不排除某个 Higgs 为零的边界分支，不分析非零 $X$ 或 $Z_a$ 的其他分支，也不容许穿过 Kähler 度量退化处把非物理延拓称作稳定解。

## 3. 固定模量后的精确重场谷底

为辨别“重场确实能松弛”与“全系统没有驻点”，取

\[
z=u^Tr,\qquad r=zu+r_\perp,\quad u^Tr_\perp=0,
\qquad A=QQ^T,
\]

并令 $H=Q^TQ$。有

\[
A^+=QH^{-2}Q^T,\quad
AA^+=A^+A=I-uu^T,\quad
h=(A^+)_{nn}=(H^{-2})_{n-1,n-1}>0.
\]

在固定 $z$ 下，product directions 的唯一极小值是 $p_a=v^2$。取 $\lambda\ne0$，无需将 $\lambda$ 送到无穷大：$p_a$ 与 $r_a$ 在所定义的开集内是独立局部坐标，且 product potential 是独立的正二次项。

定义隐藏度量 $D$ 及无量纲有限辅助场响应 $t$：

\[
B=1+c_Xu_nz,\qquad D=B+t>0,
\qquad C=\frac{Sc_X^2h}{g^2}>0.
\]

精确解由

\[
\boxed{tD^2=C,\qquad t>0,\quad D>0}
\]

给出，其重场坐标为

\[
\boxed{r=zu+\frac{Sc_X}{g^2D^2}A^+e_n.}
\]

对任意实 $B$，正分支唯一：等价方程 $D^2(D-B)=C$ 在 $D>\max(B,0)$ 上严格递增，左端由零增长到无穷。固定 $z$ 的势对 $r_\perp$ 严格凸，因为

\[
\nabla_r^2V=g^2A+\frac{2Sc_X^2}{D^3}e_ne_n^T
\]

在 $u^\perp$ 上正定。这里的“唯一”指声明区域内的约束实模量最小值，未包含被保留的全局相位和其他分支。

松弛重场后的势为

\[
\boxed{V_{\rm eff}(z)=\frac{S}{D}
+\frac{g^2t^2}{2c_X^2h}.}
\]

包络求导无需忽略重场响应：

\[
\boxed{\frac{dV_{\rm eff}}{dz}=-\frac{Sc_Xu_n}{D^2}\ne0,}
\qquad
\boxed{\frac{d^2V_{\rm eff}}{dz^2}
=\frac{2S(c_Xu_n)^2}{D^2(D+2t)}>0.}
\]

所以这个势可以单调而凸，正的局部二阶导数不意味着存在稳定极小值。以 $c_X>0$ 为例，$u_n<0$，故 $V_{\rm eff}'(z)>0$，力指向更小的 $z$。

势在上述 $(r,p)$ 坐标中的表达式不显含 $k$；这一事实不等于已经完成 light-field kinetic metric 的匹配，也不等于可以在不存在真空的背景上指定一个不受模量运动影响的物理时钟质量。低能几何和动力学需另行分析。

## 4. 受控域边界比形式上的无穷远更重要

若纯代数地将 $B\to+\infty$，则

\[
t\simeq C/B^2,\qquad
V_{\rm eff}\simeq S/B\to0.
\]

这是一个形式上的大场延拓。它没有保证全程留在正定动能和 EFT 场幅有效区间内。例如在 $k=0,c_I,c_X>0$ 时，沿 $z\to-\infty$ 有 $d_0\simeq u_0z$，因此 clock metric $1+c_Id_0$ 最终变为零。对用于诊断的无量纲参数 $n=4,q=3,S=0.5,c_X=0.3,g=0.7,c_I=0.2$，精确谷底的这个边界位于

\[
z\simeq-429.5637412.
\]

故受控结论是：**当前分支没有有限内部驻点，势将场推向允许域的边界或更大场幅。** 不能无条件声称它描述了一条可无限延伸的真实宇宙轨迹。

作为与上轮模型的交叉理解，额外 $U(1)$ 模型具有 $t(1+t)^2=Sc_X^2/(g^2\eta^2)$。在保持其他参数固定而 $\eta\to0$ 时，$t\sim\eta^{-2/3}$，该模型的静态驻点移向大场，不能得到原始模型中的有限驻点。这个极限同样可能在到达形式极限前先越出有效域。

## 5. 独立检查和实际失败记录

运行：

```bash
python reports/light_modulus_stationarity_20260926_checks.py
```

最终 **1435/1435 项通过**。

- 精确有理矩阵：$n=1,\ldots,5$、$q=2,3,5$，核验左零向量、投影、Moore–Penrose inverse、端点系数和无驻点方向；另作符号坐标反解及包络导数检查。
- **216 个独立非线性约束求解**：直接用 $\Phi_a=e^{x_a}$、$\widetilde\Phi_a=e^{y_a}$ 的完整幅度和 Lagrange multiplier 解 KKT 方程。每次由相同的零 log-amplitude 初值出发，未使用解析谷底作为初值。网格为 $n=1,2,4$、$q=2,3$、$S=0.01,0.5,4$、$c_X=\pm0.3$、$z=-2,0,2$、$k=0,0.4$；取 $v=1,g=0.7,\lambda=2,c_I=0.2$。
- 最大 KKT 残差为 $6.89\times10^{-11}$；最大能量误差 $2.21\times10^{-12}$；最大坐标误差 $1.40\times10^{-11}$。
- 所有被解出的检查点都在正定域内：最小隐藏度量为 $0.43089$；clock/Higgs block 的最小 Schur complement 为 $0.56304$。这验证所采样点的动能正定，不能替代其他标量分支的完整稳定性分析。
- 用独立有限差分检查包络导数，最终最大误差 $6.02\times10^{-12}$。同时通过 KKT multiplier 再次核验一阶导数的符号和数值。
- 另采样正负 $B$ 的 extended algebraic branch，保留完整动能边界诊断，未把超出有效域的点作为物理成功点。

**首次运行是 1423/1435，通过之外有 12 项失败。** 失败均来自二阶中心有限差分的截断误差：最大 $5.25\times10^{-8}$ 超过预设 $2\times10^{-8}$。将差分改为同一步长的四阶中心 stencil 后解决，物理模型、采样点、解析公式和所有容差保持原样。首次脚本和机器结果完整保留。

## 6. 对主研究问题的影响

保留轻模量是正确的低能记账，但它本身没有修好真空。在这一个明确的 microscopic model 中，右端有限 $F_X$ 沿模量产生未被 D 项平衡的力；重场谷底的正常存在不能消除它。

如果下一步添加模量势、改变端点耦合、保留其他辅助场或引入超引力项，应检查其是否能在完整正定域内抵消这里的非零一阶导数，同时重算端点静态交换和低能几何。本报告没有从现有结果挑选一项未计算的额外作用来宣称稳定化成功。
