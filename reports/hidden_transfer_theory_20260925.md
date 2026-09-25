# 隐藏破缺向慢时钟传递：指定 nilpotent 超引力模型的理论核查

## Material Passport

- 日期：2026-09-25；类型：指定模型的推导、符号复核与一手来源核查。
- 输入：上一轮移位型 Kähler 势、同一指数超势与带电质量律；本轮新增一个 $X^2=0$ 的隐藏超场以及常数超势 $W_c$。
- 权限与边界：用户已批准继续理论实验；只核查本轮明确模型，不改观测、不做新数据拟合、不推断整个超引力理论类是否可行。
- 核验命令：`timeout 120 python reports/hidden_transfer_theory_20260925_checks.py`。本次 29/29 个符号恒等式通过，退出码 0；无执行异常。它们是内部代数检查，不是独立同行评审或物理真实性证明。
- 输出：[符号脚本](hidden_transfer_theory_20260925_checks.py)、[检查结果](hidden_transfer_theory_20260925_checks.json)、[来源与读取范围](hidden_transfer_theory_20260925_sources.json)。数值预算与背景由本轮实验目录另外记录。

**结论先行：这个最小隐藏区并不会自动与慢时钟隔离。** 一般相位产生与 $m_G\sqrt U$ 成正比的树级交叉项；带电场还收到 $d\simeq m_G^2$ 和 $|\mathcal B|\simeq m_G M$。但不能把这一结果写成“任何相位都不能保持树级慢变化”：恰好正交相位的重横向场可以移动到很近的谷底，留下有限的指数势修正。该特殊树级情况仍需面对带电圈图、相位保护和有效理论截止尺度。

## 1. 明确作用量与新增假设

用约化普朗克质量 $M_p$，记 $F=F_\chi$、$\epsilon=F^2/M_p^2$，取

$$
Z=\frac{F}{\sqrt2}(\chi+iy),\qquad X^2=0,
$$
$$
K=-\frac12(Z-\bar Z)^2+X\bar X+|Q_+|^2+|Q_-|^2,
$$
$$
W=W_\chi(Z)+W_c+fX+M(Z)Q_+Q_-,\qquad
W_c=m_G M_p^2e^{i\theta},\quad m_G\geq0,
$$
$$
W_\chi=-\frac{F\sqrt{U_i}}{\sqrt2}
\exp\left[-\frac{\sqrt2(Z-Z_i)}F\right],\qquad
M(Z)=m_\infty\sqrt{1+\frac{\xi}{(\sqrt2Z/F)^2}}.
\tag{H1}
$$

质量的全纯分支限定在此前正实轴附近；这里没有引入 $XQ_+Q_-$、$X\bar X Q\bar Q$ 或别的可调直接耦合。规范全纯高能边界维持上一轮选择。$K$ 的实方向移位性质不是 $W$ 或整个作用量的精确对称性，不能由此宣称所有隐藏区作用都受保护。

规定一次常数真空能条件

$$
|f|^2-3m_G^2M_p^2=\rho_\Lambda.
\tag{H2}
$$

每个 $m_G$ 对应由 (H2) 确定的 $f$；这是一族明确边界条件。只抵消常数项，不抵消任何依赖 $\chi,y$ 的项，不在演化途中重调。原 FLRW 计算外加的同一 $\rho_\Lambda$ 应由本势替换，不能再次加进 Friedmann 方程。辐射、尘埃仍是现象学流体，未得到完整局部超对称物质作用量。

nilpotent 约束不能在求导前把 $X$ 的全部超场贡献删掉。其纯玻色势先用常规超引力 F 项，保留 $|D_XW|^2$，再取最低标量 $x=0$。约束解需要 $F^X\neq0$；这一规则及其适用边界可查 [Dall’Agata 与 Zwirner，§2.2，式 (16)—(17)](https://arxiv.org/html/1411.2605v3)。本轮只约束 $X^2=0$，没有额外约束移除 $Z$ 的费米子或横向标量。

## 2. 完整玻色势与轴上受力

由标准 F 项势

$$
V=e^{K/M_p^2}\left(|D_ZW|^2+|D_XW|^2
-\frac{3|W|^2}{M_p^2}\right)_{X=Q_\pm=0}
\tag{H3}
$$

出发。其普通超引力定义与 Kähler 协变导数见 [Martin，§7.6](https://arxiv.org/html/hep-ph/9709356v7)。设

$$
U=U_i e^{-2(\chi-\chi_i)},\quad b=1-\frac32\epsilon,
\quad C=\sqrt2 Fm_G\sqrt U,\quad T=\theta+y.
$$

直接计算得到

$$
\boxed{
V(\chi,y)=e^{\epsilon y^2}\left[
U(b+\epsilon^2y^2)+\rho_\Lambda+2F^2m_G^2y^2
+C\{(3-2\epsilon y^2)\cos T+2y\sin T\}\right].}
\tag{H4}
$$

特别是 $\rho_\Lambda$ 在离轴处乘有 $e^{\epsilon y^2}$，所以它不是可以在整个二维场空间独立抽出的裸宇宙学常数。即使 $m_G=0$，这一隐藏来源也给横向方向新增曲率；在 $y=0$ 的背景则与旧外加 $\rho_\Lambda$ 相同。

在 $y=0$：

$$
V=bU+\rho_\Lambda+3C\cos\theta,
$$
$$
V_\chi=-2bU-3C\cos\theta,\qquad
V_y=-C\sin\theta,
\tag{H5}
$$
$$
V_{\chi\chi}=4bU+3C\cos\theta,\qquad
V_{\chi y}=C\sin\theta,
$$
$$
V_{yy}=2\epsilon(1-\epsilon/2)U+2\epsilon\rho_\Lambda
+4F^2m_G^2+(1+2\epsilon)C\cos\theta.
\tag{H6}
$$

正则场为 $F\chi,Fy$；正则质量平方是相应 Hessian 除以 $F^2$。

| 相位 | 新增轴向梯度 | 横向梯度 | 轴是否为不变轨道 |
|---|---:|---:|---|
| $0$ | $-3C$ | $0$ | 是 |
| $\pi/2$ | $0$ | $-C$ | 否 |
| $\pi$ | $+3C$ | $0$ | 是 |

所以只取 $\theta=\pi/2$ 并强制 $y=0$ 会漏掉真实运动方程。初始交叉项的二维梯度模为

$$
\sqrt{(\Delta V_\chi)^2+V_y^2}
=C\sqrt{1+8\cos^2\theta}\ge C.
\tag{H7}
$$

它给出对原轴轨迹的受力预算，**不是**沿所有可能新轨迹的不可行性证明。尤其对重场，初始梯度大与远离原点的位移大并不等价。

## 3. 正交相位的重场谷底：必须保留的例外

当 $m_G\gg\sqrt U/F$，并且 $m_G^2\gg\rho_\Lambda/M_p^2$ 时，在很小的 $y$ 附近保留主项：

$$
V\simeq bU+\rho_\Lambda+3C\cos\theta
-Cy\sin\theta+2F^2m_G^2y^2.
$$

局部谷底在

$$
y_*\simeq\frac{\sqrt{2U}}{4Fm_G}\sin\theta,
\qquad m_y^2\simeq4m_G^2.
\tag{H8}
$$

代回给出渐近的绝热势

$$
V_{\rm eff}\simeq\rho_\Lambda+3C\cos\theta+
\left(b-\frac14\sin^2\theta\right)U.
\tag{H9}
$$

因此恰好 $\theta=\pi/2$ 时，树级指数形状仍可保留，系数从 $b$ 变为 $b-1/4$。这既不是原模型完全不变，也不是树级必然失败。式 (H9) 是局部重场消去后的势；从指定初态出发的短时振荡、能量转移和何时可取绝热极限还需分别计算。未展示的快速振荡不能被当作已经积分过。

偏离正交相位时仍有 $3C\cos\theta$。要使其轴向力相对 $2bU$ 小于给定预算 $\delta$，需

$$
|\cos\theta|\lesssim
\frac{2\delta b\sqrt U}{3\sqrt2 Fm_G}.
\tag{H10}
$$

这个精确相位选择是否由对称性保护尚未说明；也没有计算其量子稳定性。更一般的 no-scale、隔离几何或相关超势结构会改变上述跨项，不能把指定相加型模型的结果外推成一般排除。

## 4. 完整重带电谱：隐藏破缺也进入 $d$ 和 $\mathcal B$

在任意背景上记 $w=W_\chi+W_c$、$A=D_Zw$。把 (H3) 直接展开到 $Q_\pm$ 二阶，定义

$$
s=e^{K_0/M_p^2}|M|^2,\qquad
m_{3/2}^2=e^{K_0/M_p^2}\frac{|w|^2}{M_p^4},\qquad
d=m_{3/2}^2+\frac{V}{M_p^2},
$$
$$
\mathcal B=e^{K_0/M_p^2}\left[
A^*\left(M_Z+\frac{K_ZM}{M_p^2}\right)
-\frac{w^*M}{M_p^2}\right],
\qquad m_\pm^2=s+d\pm|\mathcal B|.
\tag{H11}
$$

$m_\psi^2=s$。$f$ 为常数且 $M$ 不依赖 $X$，故没有额外 $F_X M_X$ 双线性项，但 $|f|^2$ 通过 $V/M_p^2$ 进入共同标量质量。常数真空能调节不能删除 $d$ 中的 $m_{3/2}^2$。这些结构与普通隐藏区传递的一般公式一致，见 [Brignole、Ibáñez 与 Muñoz，式 (11)、(13)、(14)](https://arxiv.org/html/hep-ph/9707209v3)；此处系数由本模型直接展开得到，不把该文的驻定真空条件套到滚动宇宙背景。

实轴上 $r=\partial_\chi\ln M=-\xi/[\chi(\chi^2+\xi)]$，于是

$$
s=M^2,
$$
$$
\boxed{d=m_G^2+\frac{(1-\epsilon)U+\rho_\Lambda+2C\cos\theta}{M_p^2}},
$$
$$
\boxed{\mathcal B=M\left[
\frac{\sqrt{2U}}F\left(r+\frac\epsilon2\right)
-m_Ge^{-i\theta}\right]}.
\tag{H12}
$$

其中

$$
m_{3/2}=\left|m_Ge^{i\theta}
-\frac{F\sqrt U}{\sqrt2M_p^2}\right|
\tag{H13}
$$

才是轴上的局部引力微子质量；$m_G$ 只是常数超势的参数，二者在慢尺度附近不能混用。物理重场质量平方需要正性检查；不能只凭共同软质量为正就省略 $\mathcal B$。

在 $m_G\gg\sqrt U/F$ 且 $m_G\ll M$ 的重带电弱分裂极限，有 $d\simeq m_G^2$、$|\mathcal B|^2\simeq m_G^2s$。指定有限平直局部一圈核的主阶为

$$
V_1\simeq\frac{m_G^2s}{16\pi^2}(3L-2),\qquad
L=\ln(s/\mu^2),
$$
$$
\partial_\chi V_1\simeq
\frac{m_G^2\partial_\chi s}{16\pi^2}(3L+1).
\tag{H14}
$$

在 $\mu^2=s_i$、初始 $L=0$ 处也有非零力。若只加入共同软质量而遗漏 $\mathcal B$，可能错误得到初始力消失。原点的一个常数减法不能改变任何力。完整预算仍应使用 (H12) 的全部导数，不把渐近式当成全区间等式。

即使正交相位的树级谷底保留指数形状，这个 $m_G^2s$ 圈图尺度仍在。它是指定模型与有限减法边界下的量级诊断；有场依赖的重整化反项可改变有限势，因此不能宣称它在所有 UV 完成里都是不可取消的可观测量。当前没有引入抵消势斜率、曲率的自由调节。还未包括全部引力圈图、隐藏区 UV 阈值、规范异常破缺传递和真实标准模型。固定高能规范函数只令直接树级规范质量消失，并不消除量子异常传递；相关普通例子见 [Martin，§7.8](https://arxiv.org/html/hep-ph/9709356v7)。

## 5. nilpotent 约束和 100 GeV 阈值的使用域

本模型

$$
|F^X|=e^{\epsilon y^2/2}|f|,
\qquad |f|^2=3m_G^2M_p^2+\rho_\Lambda>0.
\tag{H15}
$$

因此有限背景上约束的分母非零；$m_G=0$ 也不使 $X$ 方向恢复超对称，因为 $\rho_\Lambda>0$。但非零 F 项只是约束解存在的必要条件之一，不证明有效理论可以用到任意高能。

一个必要的参数尺度标志是

$$
\Lambda_{\rm VA}\sim\sqrt{|F^X|},\qquad
\Lambda_{\rm VA}|_{y=0}\sim
(3m_G^2M_p^2+\rho_\Lambda)^{1/4}.
\tag{H16}
$$

这里的等号只能当作选定的量级标志；精确强耦合系数、sgoldstino 质量与真正 UV 阈值未给定。nilpotent 非线性理论的微扰幺正范围由破缺 F 项限制，背景能量密度的四次方根也不能与传播模的实际能量混为一谈，参见 [Dall’Agata 与 Zwirner，§5](https://arxiv.org/html/1411.2605v3)。[Argurio 等，§2，式 (3)—(6)](https://arxiv.org/html/1705.06788v2)进一步区分了 sgoldstino 质量、生成其质量的 EFT 尺度和真正低能解耦域；本轮没有提供这些更高能自由度。

若要**由当前 nilpotent EFT 本身**计算 $m_i=100\,\mathrm{GeV}$ 粒子的匹配，需要该质量及有关外部尺度位于其受控域内。忽略极小 $\rho_\Lambda$ 时，$\Lambda_{\rm VA}\gtrsim m_i$ 给出仅作为必要量级标志的

$$
m_G\gtrsim\frac{m_i^2}{\sqrt3M_p}
\simeq2.37\times10^{-6}\,\mathrm{eV}.
\tag{H17}
$$

要求 $\Lambda_{\rm VA}\ge10m_i$ 会把此质量参数标志提高 100 倍。**通过这个标志仍不足以证明存在受控 UV 完成。** 在低 $m_G$ 分支，宇宙背景的低频玻色演化仍可作为形式低能测试，但 100 GeV 的显式带电圈图超出本 nilpotent EFT 已说明的适用域；它只能标为形式预算，不能一边选极小隐藏破缺保护时钟，一边把该重场阈值宣称为已受控的完整计算。UV 理论可预先匹配成低能 Wilson 系数，但本轮没有指定它。

## 6. 可用于论文的结论边界

本轮把“隐藏破缺会不会影响慢时钟”变成了可计算问题：相加型最小隐藏区既给树级交叉项，也给带电软分裂；常数真空能调整不会同时消掉这些反馈。正交相位存在特殊树级谷底，必须明确保留，不能为了得到整齐结论忽略它。更强的下一步条件是同一个明确作用量需要同时说明相位/跨项保护、重带电圈图保护与 nilpotent 适用域。

这不是对逆对数思想、所有超引力完成或所有隔离机制的排除。作为例子，Brignole 等在同一来源的 no-scale 部分给出不同的软项结构；那是需要另行计算的新 $K,W$，不能未经推导借来抵消本模型中的项。本轮也未从黎曼动力学导出 (H1)，未接入真实电磁物质，未证明新的观测常数变化。
