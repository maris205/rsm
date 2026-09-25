# 最小重带电对的规范匹配：全纯边界、物理边界与 Kähler 抵消

2026-09-25。针对 `sugra_shift_v01` 的有界理论核验；22项符号检查通过。没有观测拟合，也没有加入电子或完整标准模型。

**结论：在本轮明确指定的最小 U(1) 模型中，如果固定的是 UV 全纯规范函数，重带电对完全退耦后的低能规范动能系数不保留显式的 $K/M_{\rm Pl}^2$ 项。此前仅用物理质量算出的标准型候选符号翻转，属于固定 UV 物理耦合的另一种边界。两者不能混用。** 当前超引力质量分裂留下极小的非全纯阈值修正；这一结论不等于现实电磁耦合已经完成匹配。

## 1. 范围、约定和原始公式

采用 Einstein 标架，规范项为 $-F_{\mu\nu}F^{\mu\nu}/(4g^2)$，$\alpha=g^2/(4\pi)$。背景中 $Q_\pm=0$，电荷分别为 $\pm1$，超势含 $M(Z)Q_+Q_-$。用 $\mathcal Z_\pm$ 表示带电场动能系数，避免与中性时钟超场 $Z$ 混淆。定义

$$
k=K_0/M_{\rm Pl}^2,\qquad \mathcal Z_p=\mathcal Z_+\mathcal Z_-,
\qquad s=m_\psi^2=\frac{e^k|M|^2}{\mathcal Z_p}.
\tag{G1}
$$

生产模型取 $\mathcal Z_\pm=1$；暂留它们以检查 Konishi 项。所有差分 $\Delta X=X(Z)-X(Z_0)$ 以各轨迹今天为参考，比较的 UV 和 IR 动量均在**固定物理单位**下取值。UV 满足 $m_\pm,m_\psi\ll\mu_{\rm UV}\ll M_{\rm Pl}$，IR 满足 $p\ll m_\pm,m_\psi$。不采用随场变化的测量单位或截止。

来源是 [Kaplunovsky–Louis, *Field Dependent Gauge Couplings in Locally Supersymmetric Effective Quantum Field Theories*, hep-th/9402005v2](https://arxiv.org/pdf/hep-th/9402005v2)。其式 (2.21)、(3.3) 定义

$$
b=\sum_r n_rT(r)-3T(G),\qquad c=\sum_r n_rT(r)-T(G).
\tag{G2}
$$

由式 (3.4)、(3.7) 移项，在一圈阶得到

$$
g^{-2}(\mu)=\operatorname{Re}f_{\rm UV}
+\frac{b}{16\pi^2}\ln\frac{\Lambda^2}{\mu^2}
+\frac{c}{16\pi^2}k
+\frac{T(G)}{8\pi^2}\ln g^{-2}
-\sum_r\frac{T(r)}{8\pi^2}\ln\det\mathcal Z_r .
\tag{G3}
$$

在此 U(1) 中 $T(G)=0$，每个单位电荷手征超场给 $T=1$，所以 $b_{\rm UV}=c_{\rm UV}=2$。低能不剩任何带电场，$b_{\rm IR}=c_{\rm IR}=0$。重对是非 Higgs 质量，规范对称性没有破缺。于是固定尺度间的差分为

$$
\Delta g_{\rm UV}^{-2}
=\Delta\operatorname{Re}f_{\rm UV}
+\frac{\Delta k}{8\pi^2}
-\frac{\Delta\ln\mathcal Z_p}{8\pi^2}.
\tag{G4}
$$

这里 $+ck/(16\pi^2)$ 已包含相应超 Weyl/Kähler 匹配；不能另加一份相同异常，也不能把它与阈值贡献视为不同圈阶。$\mathcal Z$ 项为负号。式 (G3) 的超对称原式并不是“滚动、破缺背景中完整有效作用量精确等于此式”的声明：本应用取重粒子局部、绝热展开的首阶规范动能系数，并在下一节明确加入已给定的小质量分裂。完整引力圈图、导数算符、隐藏区和实际轻带电粒子仍未处理。

## 2. 用物理质量退耦后，显式 Kähler 项如何抵消

单位电荷复标量与 Dirac 费米子的普通一圈系数分别为 $1/3$、$4/3$。因此重对的两个复标量和一个 Dirac 粒子满足总系数 $2$，与式 (G2) 一致。设物理谱为

$$
m_\pm^2=s+d\pm h,\qquad h=|\mathcal B|,\qquad
R=(1+d/s)^2-(h/s)^2>0.
\tag{G5}
$$

固定有限减法方案、远离各阈值时，局部质量对数贡献为

$$
\Delta g_{\rm IR}^{-2}-\Delta g_{\rm UV}^{-2}
=-\frac1{16\pi^2}\left[
\frac13\Delta\ln m_+^2+\frac13\Delta\ln m_-^2
+\frac43\Delta\ln s\right]
=-\frac{\Delta\ln s}{8\pi^2}-\frac{\Delta\ln R}{48\pi^2}.
\tag{G6}
$$

质量对数的负号也可由 $d(g^{-2})/d\ln\mu=-b/(8\pi^2)$ 检查：在固定高能耦合下提高退耦质量，低能逆耦合会减少。代入式 (G1)、(G4) 得

$$
\boxed{\displaystyle
\Delta g_{\rm IR}^{-2}
=\Delta\operatorname{Re}f_{\rm UV}
-\frac{\Delta\ln|M|^2}{8\pi^2}
-\frac{\Delta\ln R}{48\pi^2}.}
\tag{G7}
$$

$\Delta k$ 和 $\Delta\ln\mathcal Z_p$ 均抵消。**这没有撤销“阈值必须使用物理质量”的规则；正是使用物理质量并加入同圈阶异常项后，才得到此抵消。** 在严格超对称极限 $d=h=0$，也可不依靠逐粒子阈值，直接用 KL 式 (3.22)—(3.29) 得到

$$
f_{\rm IR}=f_{\rm UV}-\frac1{4\pi^2}\ln\frac{M}{M_*}
+\hbox{场无关常数}.
\tag{G8}
$$

$M_*$ 是固定参考尺度；式 (3.29) 的 $c_{\rm UV}-c_{\rm IR}=2$ 正好消去重对物理质量中的 Kähler 因子。这是此特定场谱的结果；若有轻带电场、非阿贝尔矢量或别的场依赖 UV 项，式 (G7) 需要相应改写。

令 $r=d/s$、$v=h/s$，小分裂给

$$
\ln R=2r-r^2-v^2+\frac23r^3+2rv^2+O((|r|+|v|)^4).
\tag{G9}
$$

首个共同质量分裂项是 $O(d/s)$，对称的 $\pm h$ 分裂首项为 $O(h^2/s^2)$。先前谱中 $d/s\sim10^{-90}$、$h/s\sim10^{-46}$，修正极小但不能通过双精度中直接形成 `1+r+v` 来可靠计算；应使用稳定展开或高精度。若对时间求导，$s,d,h$ 全部随场变化，尤其

$$
(\ln R)'=\frac{2(1+r)r'-2vv'}{R}.
\tag{G10}
$$

有限分裂的式 (G6) 是指定物理谱的局部一圈分量阈值；不能把超对称的全纯非重整化结论推广为任意大 soft 破缺下的完整理论。展示极小修正的数值精度也不意味着未知高圈、曲率或隐藏区修正具有同等小的误差。

## 3. 两种边界代表不同的理论输入

**A：固定 UV 全纯规范函数。** 在已声明的 Kähler 表示中令 $f_{\rm UV}=f_0$，则式 (G7) 剩下全纯质量和分裂项。标准型候选中的 $e^k$ 不会独自产生先前局部阈值的反向趋势。不同候选仍具有不同的 $\chi(t)$、$d$ 和 $h$，所以曲线不必相同。

**B：固定 UV 物理耦合。** 在同一个固定物理尺度 $\mu_{\rm UV}$，要求沿实轴轨迹 $g_{\rm UV}$ 不随 $\chi$ 变化，则必须指定

$$
\Delta\operatorname{Re}f_{\rm UV}
=-\frac{\Delta k}{8\pi^2}
+\frac{\Delta\ln\mathcal Z_p}{8\pi^2}.
\tag{G11}
$$

此时式 (G7) 重新变为仅物理质量阈值的式 (G6)。这是一种可定义的条件边界，不能同时再要求同一 Kähler 表示中的 $f_{\rm UV}$ 为常数。

式 (G11) 不是整个复场空间中任意非全纯函数都可被全纯 $f$ 抵消的主张。例如标准型 $k=|Z|^2/M_{\rm Pl}^2$、$\mathcal Z_p=1$ 时，沿 $Z=\bar Z$ 可以选

$$
f_{\rm UV}(Z)=f_0-\frac{Z^2-Z_0^2}{8\pi^2M_{\rm Pl}^2}.
\tag{G12}
$$

它满足实轴边界，但离轴时不满足同样的常数条件，因为 $\operatorname{Re}Z^2\ne|Z|^2$；一般 $K$ 是非调和函数。该 $Z$ 依赖的规范函数还会带来其它局部超对称相互作用。本轮只作规范动能匹配，不把它们的完整反馈宣称为已求解。移位型候选在所考察实轴上 $k=0$，故 A、B 的这一差别消失；这不是二者在整个复场空间上都相同。

今天的单点 $\alpha_0$ 校准不能选择 A 或 B：它只固定一个常数，两种边界的差异在于整个场依赖函数。所谓“常数全纯 $f$”本身也依赖所指定 Kähler 表示，转换表示时必须连同下一节的异常变换一起携带，不能重新把变换后的 $f$ 人为设为常数。

## 4. Kähler 和 Konishi 变换检查

写 $j(Z)=J(Z)/M_{\rm Pl}^2$。KL 式 (3.9)、(3.10)、(3.13) 要求

$$
k\longmapsto k+j+\bar j,\qquad
M\longmapsto e^{-j}M,\qquad
f_{\rm UV}\longmapsto f_{\rm UV}-\frac{c_{\rm UV}}{8\pi^2}j.
\tag{G13}
$$

本例 $c_{\rm UV}=2$，式 (G8) 中两项的变换正好相消；低能 $c_{\rm IR}=0$ 与此一致。物理质量 $s$ 本身也不变。从标准型 $K=|Z|^2$ 变到移位型 $K=-(Z-\bar Z)^2/2$ 的真正表示变换取 $J=-Z^2/2$，因此还需

$$
W\longmapsto e^{Z^2/(2M_{\rm Pl}^2)}W,\qquad
f_{\rm UV}\longmapsto f_{\rm UV}+\frac{Z^2}{8\pi^2M_{\rm Pl}^2}.
\tag{G14}
$$

先前保持相同 $W$ 而更改 $K$ 的两个候选仍是不同物理模型；本轮显式 $k$ 抵消并未使它们成为等价模型。

若仅重定义带电场 $Q'_\pm=\Upsilon_\pm(Z)Q_\pm$，则

$$
\mathcal Z'_\pm=\frac{\mathcal Z_\pm}{|\Upsilon_\pm|^2},\qquad
M'=\frac{M}{\Upsilon_+\Upsilon_-},\qquad
f'_{\rm UV}=f_{\rm UV}-\frac{\ln(\Upsilon_+\Upsilon_-)}{4\pi^2}.
\tag{G15}
$$

最后一项是 KL 式 (3.19) 的 Konishi 系数和符号。式 (G1)、(G4)、(G7) 在同步变换下均不变。省略它才会产生对带电场归一化的虚假依赖。

## 5. 对逆对数平方构想的实际含义

忽略极小分裂并取边界 A，如果所假设的质量律为 $|M|^2\propto1+\xi/\chi^2$，那么

$$
\Delta g_{\rm IR}^{-2}
=-\frac1{8\pi^2}\ln\frac{1+\xi/\chi^2}{1+\xi/\chi_0^2}.
\tag{G16}
$$

低能读数首先是**逆平方质量响应的对数**。应区分两种展开。若把 $|\xi/\chi^2|\ll1$ 本身作为小参数，主阶是 $\xi(\chi^{-2}-\chi_0^{-2})$；当前输入的 $\xi/\chi^2$ 接近1，不能使用这一小参数假设。但这并不是获得逆平方小变化响应的必要条件：相对今天定义

$$
X=\frac{\xi(\chi^{-2}-\chi_0^{-2})}{1+\xi/\chi_0^2},\qquad
\ln\frac{1+\xi/\chi^2}{1+\xi/\chi_0^2}=\ln(1+X).
$$

当 $|X|\ll1$ 时仍可展开为 $X+O(X^2)$，即带固定参考系数 $\xi/(1+\xi/\chi_0^2)$ 的逆平方差，而无需 $\xi/\chi^2$ 很小。本轮轨迹 $\max|X|$ 约为0.0155—0.0201，因此存在可量化的小变化近似，不能把它说成严格相等；[完整轨迹报告](../experiments/gauge_matching_v01/reports/results_cn.md)分别计算质量对数近似与时间映射误差。如果动力学另外给出近似 $\chi\sim\ln t$，式 (G16) 才形成与 $1/\ln^2t$ 有关的条件函数。这仍没有从黎曼结构推出质量律或真实相互作用。

以 $\alpha_0$ 作今天校准，可写

$$
\frac{\alpha(t)}{\alpha_0}-1
=\left[1+4\pi\alpha_0\Delta g_{\rm IR}^{-2}(t)\right]^{-1}-1.
\tag{G17}
$$

此处分式是把一圈逆耦合结果转换为读数的代数关系，并不意味着已经计算全部高圈。由于 IR 场谱为空带电部门，这里的 $\alpha$ 只是固定电荷归一化下的**条件 U(1) 规范动能读数**。没有电子、标准模型阈值、实验测量算符及其物质耦合，就不能把它直接命名为已预测的实验精细结构常数，也不能直接将 ppm 与天文或原子钟数据相减作拟合。

## 6. 读取和核验记录

- KL v2 的来源身份由 arXiv 摘要页与 PDF 核对；§2.3 式 (2.20)—(2.22) 在印刷页19；§3.1 式 (3.3)—(3.8) 在印刷页24—26；§3.2 式 (3.13)、(3.19) 在印刷页28、30；§3.3 式 (3.21)—(3.29) 在印刷页32—36。读的是这些明确段落，未声称阅读全文。
- `pdftotext -layout` 提取原文；对式 (3.27)、(3.28) 另以本地 PDF 第35页图像核对符号、括号和系数。预检工具返回 `UNAVAILABLE`，原因是环境缺少 `pypdf`，故不把自动结构预检说成通过；本地渲染与印刷页码另行核对。网页截图一次超时，随即使用已缓存同哈希 PDF 渲染，没有替换来源。
- [符号脚本](gauge_matching_theory_20260925_checks.py) 独立实现上述代数，不导入主匹配程序；[22/22结果](gauge_matching_theory_20260925_checks.json) 包括 Kähler/Konishi 不变性、K 抵消、两种边界、非全纯限制与分裂展开。命令：`python reports/gauge_matching_theory_20260925_checks.py`，退出码0。
- [来源与本地材料哈希](gauge_matching_theory_20260925_sources.json) 记录版本、读取范围和异常。上述均是内部理论复核，不是外部同行评审或自然界验证。
