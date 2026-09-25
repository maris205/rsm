# 几何隔离加入模场稳定之后：完整势、带电分裂与不一致极限

## Material Passport

- 日期：2026-09-25；研究范围：本轮指定的一个 no-scale Kähler 几何与二次模场稳定超势。用户已授权继续深入理论和数值核查。
- 输入：此前指数时钟超势 $w(Z)$、全纯带电质量 $M(Z)$；本轮新增无量纲复模场 $T$，并把 $Z,X,Q_\pm$ 放在同一个对数 Kähler 势中。
- 实际工作：独立推导完整 F 项势、带电二阶作用量、一次 uplift、强稳定渐近极限；用 SymPy 和 80 位 mpmath 作代数与局部驻点检查。没有读取本轮其他代理的结论或导入其他代理的计算代码。
- 输出：[检查脚本](noscale_stabilization_20260925_checks.py)、[机器可读结果](noscale_stabilization_20260925_checks.json)。命令：`timeout 180 python reports/noscale_stabilization_20260925_checks.py`。
- 不在范围：新宇宙背景积分、真实标准模型、全套量子超引力、观测数据拟合、对整个逆对数假设的证实或排除。文中的数值例子只验证局部渐近结构。
- 核查过程记录：首次 64/65 通过，一个辅助量检查误用了矩阵索引方向；将逆度量转置用于反变辅助量后 65/65 通过。随后增加逐项矩阵逆与有限稳定强度检查，最终 70/70 通过。势、带电谱、驻点公式未因这一索引修正改变。来源 HTML 的一次 HTTP 404 已记录，并改读同一工作的 arXiv 原始 PDF。
- 研究性质：模型推导和确定性检查，不涉及 p 值、统计显著性或新观测证据。作者仍需审阅；这里是同一 AI 系统内的独立计算角色，不是外部同行评审。

**核心结论：该几何确实消掉了共同带电软质量中的 $m_G^2$，但加入指定超势稳定后，一般会重新产生时钟跨项与 $m_GM$ 量级的带电双线性分裂。** 因而“软质量被几何隔离”是一个局部成功，不能直接写成“慢时钟已受完整保护”。还有两个需要同时保留的细节：严格强稳定极限的特殊相位仍有局部指数势谷底；有限稳定强度却会移动这个特殊相位，不能用交换极限的方法删除残余时钟力。

## 1. 明确模型与单位

记 $P=M_p^2$，$T$ 无量纲，$Z,X,Q_\pm$ 具有质量量纲，取

$$
X^2=0,\qquad
K=-3P\ln Y,
$$
$$
Y=T+\bar T-
\frac{k(Z,\bar Z)+X\bar X+|Q_+|^2+|Q_-|^2}{3P},
\qquad k=-\frac12(Z-\bar Z)^2,
\tag{S1}
$$
$$
W=w(Z)+W_c+\frac\kappa2(T-\tfrac12)^2+fX+M(Z)Q_+Q_-,
\quad W_c=m_GP e^{i\theta},\quad\kappa>0.
\tag{S2}
$$

$\kappa$ 的单位是 $\mathrm{eV}^3$。$f$ 的单位是 $\mathrm{eV}^2$。物理定义域要求 $Y>0$，并保留 $D_XW=f$ 后才取 nilpotent 最低标量 $x=0$。这里没有假定 $D_TW=0$ 或删掉模场辅助量。

时钟和质量沿用前一轮的局部分支：

$$
Z=\frac F{\sqrt2}(\chi+iy),\quad
w=-\frac{F\sqrt{U_i}}{\sqrt2}
e^{-\sqrt2(Z-Z_i)/F},\quad
M=m_\infty\sqrt{1+\frac{\xi}{(\sqrt2 Z/F)^2}}.
\tag{S3}
$$

该平方根只在此前正实轴附近指定全纯分支。它尚未从黎曼动力学微观推导。几何隔离的历史动机可查 [Berg 等，§2.1–2.2](https://arxiv.org/pdf/1012.1858)：他们明确区分了共同软质量的压低与已有超对称质量项的双线性项。这里的具体系数由 (S1)–(S2) 直接计算，并非从不同模型移植。

## 2. 未冻结模场的精确势

在 $X=Q_\pm=0$，定义

$$
J=W_T=\kappa(T-\tfrac12),\qquad
W_0=w+W_c+\frac\kappa2(T-\tfrac12)^2,
$$
$$
A=w_Z+\frac{k_ZJ}{3P},\qquad
k_Z=-(Z-\bar Z),\quad k_{Z\bar Z}=1.
\tag{S4}
$$

完整势可以压缩为

$$
\boxed{
V(T,Z)=\frac1{Y^2}\left[
|A|^2+|f|^2+\frac{Y|J|^2}{3P}
-\frac{J\bar W_0+\bar J W_0}{P}
\right].}
\tag{S5}
$$

这一等式包括了完整 $T,Z$ 度量和 $F^X$；不是把 $T$ 固定以后再凭 no-scale 名称添加或删除项。一个简明核查是

$$
K_{i\bar j}=\begin{pmatrix}
3P/Y^2&-k_{\bar Z}/Y^2\\
-k_Z/Y^2&1/Y+k_Zk_{\bar Z}/(3PY^2)
\end{pmatrix},
\tag{S6}
$$

并直接代入 $e^{K/P}(K^{i\bar j}D_iW\overline{D_jW}-3|W|^2/P)$。

当 $J=0$ 时，(S5) 退化为 $(|w_Z|^2+|f|^2)/Y^2$，$W_c$ 的显式作用确实消去。但这个事实不意味着 $T=1/2$ 是完整势的驻点。稳定超势的斜率 $J$ 可以在很小的模场位移下保持有限；正是这个量重新进入时钟和带电谱。

## 3. 带电物理谱：共同项消去，双线性项保留

带电度量为 $K_{Q\bar Q}=1/Y$，正则坐标为 $Q_c=Q/\sqrt Y$。直接展开 (S5) 的完整带电推广到二阶，得到

$$
s=m_\psi^2=\frac{|M|^2}{Y},\qquad
\boxed{d=\frac{2V}{3P}},
$$
$$
\boxed{\mathcal B=\frac1Y\left(
\bar A M_Z-\frac{\bar J M}{3P}\right)},\qquad
m_\pm^2=s+d\pm|\mathcal B|.
\tag{S7}
$$

这是任意允许 $T,Z$ 背景上的局部二阶势谱，尚未要求宇宙学驻定。它也不包括滚动背景的速度相关有效质量、非绝热粒子产生或规范圈图。

共同 $m_G^2$ 项在这一特殊对数几何中消失，是可明确指出的改进。但是 $\mathcal B$ 与共同 $d$ 是不同的参数，不能从 $d\simeq0$ 推出费米子与标量重新简并。特别是 $J$ 在稳定后一般接近 $3W_0/Y$；于是轴上有

$$
\mathcal B\simeq \bar w_ZM_Z-\frac{\bar W_0M}{P}
\quad(Y\simeq1).
\tag{S8}
$$

如果 $m_G$ 远大于慢时钟尺度，仍有 $|\mathcal B|\simeq m_G|M|$。这一结果也说明，把重带电质量 $M(Z)$ 作为显式超势项本身就需要额外保护。

## 4. 一次常数 uplift：同时解出参考模场真空

为避免沿时钟轨迹重调常数，先选一个明确的“关闭时钟”参考：$w=w_Z=0$，指定真空能为原来的 $\rho_\Lambda$。这只是校准隐藏/稳定部门的边界条件，不是一次新的宇宙演化解。

设

$$
T_*=\tfrac12+a+ib,\quad
c=\frac{W_c}{\kappa},\quad
\eta=\frac{P\rho_\Lambda}{\kappa^2},\quad
Y_*=1+2a.
$$

要求 $\partial_TV=0$ 和 $V=\rho_\Lambda$，等价于

$$
\boxed{
a-\frac32a^2-\frac12b^2
=3\operatorname{Re}c+6\eta Y_*,\qquad
b(1-a)=3\operatorname{Im}c,}
\tag{S9}
$$
$$
\boxed{
|f|^2=(1-2a)\left[
\frac{\kappa^2(a^2+b^2)}{3P}+\rho_\Lambda Y_*
\right].}
\tag{S10}
$$

取连续连接到 $a=b=0$ 的局部分支，并检查 $Y_*>0$、$|f|^2>0$、$T$ 两实方向 Hessian 为正；不能把所有代数分支都当作可用真空。求出后固定同一个 $f$，恢复 $w(Z)$，重新求 $T$ 响应，沿 $\chi$ 不再调节真空能。原先外加的 $\rho_\Lambda$ 不应再重复加入。

强稳定时，

$$
a+ib\simeq\frac{3W_c}{\kappa},\qquad
|f|^2\simeq \frac{3|W_c|^2}{P}+\rho_\Lambda,
\qquad m_T\simeq\frac\kappa{3P}.
\tag{S11}
$$

这里的物理模场质量由正则化 Hessian 得到；$\kappa$ 本身不是质量。普通强稳定模型中，“很小的 $D_TW$”并不意味着“很小的 $W_T$”，可参见 [Dudas 等，§2.1，式 (9)–(12)](https://arxiv.org/html/1209.0499v2)。该论文使用的 uplift Kähler 势与 (S1) 不同，所以不能把其通用软质量结论直接替换 (S7)。

## 5. 强稳定极限：先消去有限辅助量，再看时钟

当 $\kappa\to\infty$、$W_c,w$ 固定时，模场位移趋零，但 $J=\kappa(T-1/2)$ 不趋零。令

$$
Y_0=1-\frac{\epsilon y^2}{3},\quad
D=1+\frac{\epsilon y^2}{3},\quad
\epsilon=F^2/P.
$$

对 (S5) 的有限 $J$ 配方并取驻点，得到

$$
J_\infty=\frac{3(W_c+w)-k_{\bar Z}w_Z}{D},
\tag{S12}
$$
$$
\boxed{
V_\infty=\frac1{Y_0^2}\left[
|w_Z|^2+|f|^2-
\frac{|3(W_c+w)-k_{\bar Z}w_Z|^2}{3PD}
\right].}
\tag{S13}
$$

这个消去条件对应 $F^T\to0$。在离轴处由于 $T,Z$ 动能混合，不应再把它简化成与混合无关的 $D_TW=0$。

设 $U=U_i e^{-2(\chi-\chi_i)}$、$b_0=1-3\epsilon/2$、$C=\sqrt2Fm_G\sqrt U$，轴上恢复

$$
V_\infty(\chi,0)=b_0U+\rho_\Lambda+3C\cos\theta,
\qquad (V_\infty)_y(\chi,0)=-C\sin\theta.
\tag{S14}
$$

一般相位的跨项依旧存在。离轴的完整强稳定极限为

$$
V_\infty=\frac1{Y_0^2}\left[
|f|^2-\frac{3m_G^2P}{D}
+U\left(1-\frac{\epsilon(9+4y^2)}{6D}\right)
+\frac{C}{D}\{3\cos(y+\theta)+2y\sin(y+\theta)\}
\right].
\tag{S15}
$$

这与上轮相加型 Kähler 模型的离轴势不同。因此也不能照抄上轮重横向场质量或谷底修正。

## 6. 必须保留的特殊相位，以及有限稳定强度的限制

先严格取上述强稳定极限，再取 $m_G\gg\sqrt U/F$ 且 $m_G^2\gg\rho_\Lambda/P$，横向主阶为

$$
V_\infty\simeq b_0U+\rho_\Lambda+3C\cos\theta
-Cy\sin\theta+F^2m_G^2y^2.
$$

所以

$$
y_*\simeq\frac{\sqrt U}{\sqrt2Fm_G}\sin\theta,
\qquad m_y^2\simeq2m_G^2,
$$
$$
V_{\rm valley}\simeq\rho_\Lambda+3C\cos\theta+
\left(b_0-\frac12\sin^2\theta\right)U.
\tag{S16}
$$

恰好正交相位留下 $(b_0-1/2)U=0.49985U$，指数形状保留。这是局部重场谷底的渐近结论，不是完整 FLRW 轨迹、全局吸引子或由给定初态必然到达的证明。

**有限 $\kappa$ 的限制更严格。** 在 $\theta=\pi/2$、$\rho_\Lambda$ 可忽略的参考真空，(S9) 给出

$$
a=\frac92\left|\frac{W_c}{\kappa}\right|^2+\cdots,
\qquad
\operatorname{Re}J_*=\frac92\frac{|W_c|^2}{\kappa}+\cdots.
\tag{S17}
$$

恢复很小的实 $w$ 后，由驻点的包络定理，沿 $y=0$ 有

$$
\left.\frac{\partial V_{\rm eff}}{\partial\operatorname{Re}w}\right|_{w=0}
=-\frac{2\operatorname{Re}J_*}{P Y_*^2}.
\tag{S18}
$$

因此有限 $\kappa$ 的 $W_c$ 正交相位仍有 $m_G^2/m_T$ 量级的线性时钟驱动。模场相对位移非常小，不代表它对极慢时钟的相对力足够小。以“$m_T\gg m_G$”替代实际慢时钟预算，会错过这种尺度差异。

可以选择不同的、依赖稳定参数的相位使 $\operatorname{Re}J_*=0$。例如 $\rho_\Lambda=0$ 且 $a=0$ 时，精确参考条件是 $\operatorname{Re}c=-b^2/6$、$\operatorname{Im}c=b/3$。这是另一个需给出保护来源的关系，尚未自动由 (S1)–(S2) 保证；它也不能顺带消掉 $|\mathcal B|$。在有限 $\kappa$ 数值中必须明确区分“$W_c$ 相位为 $\pi/2$”与“稳定后 $\operatorname{Re}J_*=0$”。

## 7. 共同软项被消去之后，双线性圈图仍然在

在 $Y\simeq1$、$d$ 远小于隐藏诱导分裂、$m_G\ll|M|$ 的局部极限，$|\mathcal B|^2\simeq m_G^2s$。沿用此前明确的有限局部一圈核：

$$
V_1=\frac{f_{\rm CW}(s+d+|\mathcal B|)
+f_{\rm CW}(s+d-|\mathcal B|)-2f_{\rm CW}(s)}{32\pi^2},
$$
$$
f_{\rm CW}(x)=x^2[\ln(x/\mu^2)-3/2].
$$

仅保留双线性主项便有

$$
V_1\simeq\frac{m_G^2s}{16\pi^2}\ln\frac s{\mu^2},
\qquad
\partial_\chi V_1\simeq\frac{m_G^2s'}{16\pi^2}
\left(\ln\frac s{\mu^2}+1\right).
\tag{S19}
$$

在匹配点 $\mu^2=s$，这一主阶势值为零，**力仍不为零**。所以只展示势曲线重合会漏掉时钟问题。式 (S19) 是固定有限减法条件下的局部重带电贡献，未证明它在任意 UV 完成或任意可调场依赖反项下都是不变观测量。全套数值预算仍应使用 (S7) 并包含 $T(\chi),y(\chi)$、$M(Z)$ 和分裂项的全部导数。

## 8. 有效理论域包括稳定模场

此几何有

$$
K_{X\bar X}=1/Y,\quad |F^X|=|f|/\sqrt Y,
\quad |F_{X,\mathrm{can}}|=|f|/Y.
$$

故一个坐标正则化后的 nilpotent 尺度标志是

$$
\Lambda_{\rm VA}\sim\sqrt{|f|/Y}.
\tag{S20}
$$

这只是一阶幂次计数标志，不是已经构造出 UV 完成的证明。如果把重带电粒子和稳定模场都作为当前 EFT 内的显式自由度，至少应同时检查 $|M|/\sqrt Y$ 和 $m_T$ 是否低于已说明的 cutoff。不能无限增大 $\kappa$ 压低 (S17)，同时完全不检查 $m_T$ 是否越界。

若较高能理论先把超出 cutoff 的自由度匹配掉，则其 Wilson 系数需要另行给出；低能表达式本身不能替代这项工作。$T$、隐藏 goldstino、规范异常传递和更高维算符的量子反馈也尚未全部加入。

## 9. 独立局部数值核查

检查脚本采用 $P=1,|W_c|=1,\rho_\Lambda=10^{-12}$ 的无量纲例子，选 $\kappa=10^3,10^5,10^7$ 和三个相位。它先解 (S9)–(S10)，验证势值和两方向驻定，再在固定 $f$ 下加入很小的时钟项并重新解两模场方向。三个 $\kappa$ 不是实际 eV 预测，尤其不能拿这些单位直接评价 nilpotent EFT 的物理 cutoff。

| $\kappa$ | 正交相位 $\operatorname{Re}J_*$ | 相对 $4.5/\kappa$ | 最小物理模质量 / $(\kappa/3)$ |
|---:|---:|---:|---:|
| $10^3$ | $4.5000708769\times10^{-3}$ | $1.0000157504$ | $0.9984943580$ |
| $10^5$ | $4.5000000071\times10^{-5}$ | $1.0000000016$ | $0.9999849994$ |
| $10^7$ | $4.5000000000\times10^{-7}$ | $1.0000000000$ | $0.9999998500$ |

九个参考驻点的模场 Hessian 均为正。它们验证 (S11)、(S17) 的收敛，却也显示有限 $\kappa$ 泄漏并未严格为零。这里没有把 $T$ 正曲率误当作包括 $Z$ 在内的整个多场系统稳定，也没有把局部势驻点误当作宇宙吸引解。

## 10. 对下一步模型的具体要求

本轮允许提出的较强结论是：**对数几何对共同软质量有确切作用，但指定的超势稳定与显式 $M(Z)Q_+Q_-$ 质量项留下了未受保护的通道。** 要真正改善慢时钟的量子稳定性，下一种明确作用量至少应同时回答：

1. 稳定完成后，为什么影响时钟的 $\operatorname{Re}J$ 或对应跨项足够小，并且这一关系受结构保护？
2. 带电双线性 $\bar A M_Z-\bar JM/(3P)$ 如何受保护，同时仍有非平凡质量/规范阈值变化？
3. 隐藏破缺尺度、带电阈值和稳定模场质量是否都位于已经说明的有效理论域？
4. 残余量子异常、更高阶 Kähler 修正和匹配反项如何限定，而非仅在一个背景点取消？

这些是本模型导出的条件，尚未证明一般不可能满足。下一步可以检验不同的稳定机制或质量生成机制，但必须把新的 $K,W$ 和匹配边界写出来；换一个几何名称或重新选一个漂亮相位，不足以代替计算。
