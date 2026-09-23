# 最小作用量补全 v0.2：独立推导与适用边界

2026-09-22。范围：为 [v0.1 最小框架](../draft_v01/minimal_framework_v01.md)提供一种明确的驱动场与电磁响应实现，核对变分方程、能量交换、对数解及其非唯一性。本报告没有重新拟合观测，没有声称完成原子理论、黎曼谱到物理场的推导或完整宇宙学。

**结论：** 指数势的正则标量场可以在一个常状态方程的 scaling 阶段生成 $\chi=\ln(t/t_*)$，再由规定的电磁响应产生 $1/\ln^2(t/t_*)$。这是具有能量收支的经典有效实现。对数背景属于已知指数势解；平方指数仍是响应假设，不能说已经从黎曼零点或该作用量的其他项唯一导出。

## 1. 定义、作用量与已有研究的位置

采用 $c=\hbar=1$、度规号差 $(-+++)$，$M_{\rm Pl}=(8\pi G)^{-1/2}$。$\chi$ 无量纲，$F_\chi>0$ 为固定质量尺度，$\phi=F_\chi\chi$ 是规范归一化标量。$F_{\mu\nu}=\partial_\mu A_\nu-\partial_\nu A_\mu$ 是电磁场强，勿与 $F_\chi$ 混淆。

$$
S=\int d^4x\sqrt{-g}\left[
\frac{M_{\rm Pl}^2}{2}R
-\frac{F_\chi^2}{2}g^{\mu\nu}\partial_\mu\chi\partial_\nu\chi
-V_s e^{-2\chi}
-\frac{B(\chi)}4 F_{\mu\nu}F^{\mu\nu}
\right]+S_m[g,A,\Psi],
\tag{1}
$$

$$
\mathcal A(\chi)=1+\beta(\chi^{-2}-\chi_0^{-2}),
\qquad B(\chi)=\mathcal A(\chi)^{-1}.
\tag{2}
$$

$V_s>0$ 的量纲为质量四次方；$\beta,\chi_0,\mathcal A,B$ 无量纲。有效域取 $\chi>0$、$\mathcal A>0$，实际计算采用远离 $\chi=0$ 和 $\mathcal A=0$ 的闭区间。$\chi_0$ 是同一背景解在参照时刻的场值，$B(\chi_0)=1$ 是电磁归一化条件。这里 $S_m$ 的微观参数先不显式依赖 $\chi$；物质中的电磁束缚能及其标量响应并不会因此消失。

指数势与流体的 scaling 解已有系统分析；使用标量依赖的光子动能系数研究 $\alpha$ 变化也已有先例。本候选的研究对象是具体响应函数、它与既有数值动机的关系以及可检验限制，不能把标量—电磁作用量本身称作首创。参见 [Copeland、Liddle 与 Wands，1998，方程 (1)—(10) 与表 I](https://arxiv.org/pdf/gr-qc/9711068)，以及 [Copeland、Nunes 与 Pospelov，2004，第 2 节](https://arxiv.org/pdf/hep-ph/0307299)。下列计算按本报告的符号重新推导。

## 2. 变分方程与符号

记 $\mathcal I=F_{\mu\nu}F^{\mu\nu}$，并定义电流约定 $j^\nu=-(1/\sqrt{-g})\delta S_m/\delta A_\nu$。电磁与标量方程分别为

$$
\nabla_\mu(BF^{\mu\nu})=j^\nu,
\qquad \nabla_{[\rho}F_{\mu\nu]}=0,
\tag{3}
$$

$$
F_\chi^2\Box\chi-V_{,\chi}-\frac14B_{,\chi}\mathcal I=0,
\quad V_{,\chi}=-2V,
\quad B_{,\chi}=2\beta\chi^{-3}B^2.
\tag{4}
$$

平直 FLRW 中 $ds^2=-dt^2+a(t)^2d\mathbf x^2$，均匀场有 $\Box\chi=-\ddot\chi-3H\dot\chi$，所以

$$
F_\chi^2(\ddot\chi+3H\dot\chi)-2V
=-\frac14B_{,\chi}\langle\mathcal I\rangle.
\tag{5}
$$

右侧的负号随这里的作用量与度规约定固定。若采用含有效质量 $m_a(\chi)$ 的粗粒化物质作用量，定义 $J_\chi=(1/\sqrt{-g})\delta S_m/\delta\chi$，则式 (4) 左侧加 $J_\chi$、式 (5) 右侧加 $J_\chi$；必须避免把已纳入有效质量的电磁束缚贡献再计入一次。

Einstein 方程为 $M_{\rm Pl}^2G_{\mu\nu}=T^{\chi}_{\mu\nu}+T^{\rm em}_{\mu\nu}+T^m_{\mu\nu}$，其中

$$
T^{\chi}_{\mu\nu}=F_\chi^2\partial_\mu\chi\partial_\nu\chi
-g_{\mu\nu}\left[\frac{F_\chi^2}{2}(\partial\chi)^2+V\right],
\quad
T^{\rm em}_{\mu\nu}=B\left(F_{\mu\rho}F_\nu{}^\rho-\frac14g_{\mu\nu}\mathcal I\right).
\tag{6}
$$

## 3. 驱动能量来自哪里

均匀标量具有

$$
\rho_\chi=\frac{F_\chi^2}{2}\dot\chi^2+V,
\qquad p_\chi=\frac{F_\chi^2}{2}\dot\chi^2-V.
\tag{7}
$$

式 (5) 乘以 $\dot\chi$ 给出精确的标量能量方程

$$
\dot\rho_\chi+3H(\rho_\chi+p_\chi)
=-\frac14\dot B\langle\mathcal I\rangle\equiv Q_\chi.
\tag{8}
$$

把电磁场与带电物质合为 rest 子系统、对相同背景作一致平均，则

$$
\dot\rho_{\rm rest}+3H(\rho_{\rm rest}+p_{\rm rest})=-Q_\chi.
\tag{9}
$$

这同时包含场间交换与宇宙膨胀项，没有无来源的能量增加。存在显式 $J_\chi$ 时，$Q_\chi$ 另加 $J_\chi\dot\chi$，rest 端加相反项。协变形式也可核对：

$$
\nabla_\mu T_\chi^{\mu\nu}=
\left(\frac14B_{,\chi}\mathcal I-J_\chi\right)\nabla^\nu\chi,
\quad
\nabla_\mu(T_{\rm em}^{\mu\nu}+T_m^{\mu\nu})=
-\left(\frac14B_{,\chi}\mathcal I-J_\chi\right)\nabla^\nu\chi.
\tag{10}
$$

在局部正交标架中，$\mathcal I=2(|\mathbf B_{\rm mag}|^2-|\mathbf E|^2)$。自由传播的理想平面波是零不变量场，$\mathcal I=0$；几何光学极限下各向同性自由辐射的平均源也可忽略。**有电磁能量不等于有非零 $F^2$ 源，反过来辐射源为零也不等于所有物质源为零。** 静电束缚、磁场、介质和带电粒子的有效质量贡献需要单独计算；慢变 $B$ 条件及粗粒化误差也要检查。[标量—电磁宇宙学中的源项处理](https://arxiv.org/pdf/hep-ph/0307299)。

作用量没有人为写入的显式宇宙时间函数。指定一个背景解后，探针所见的 $B[\chi(t)]$ 才成为非自治驱动。这给出了 v0.1 外部参数历史的一种物理来源；尚未证明原有离散更新方程就是该场论的离散化。

## 4. 对数背景的直接代入验证

在电磁及有效物质标量源可忽略时，先试取

$$
\chi(t)=\ln(t/t_*),\quad \dot\chi=t^{-1},\quad
\ddot\chi=-t^{-2},\quad H=h/t,
\tag{11}
$$

其中 $h$ 为常数、$t>t_*$。代入式 (5) 得

$$
\frac{F_\chi^2(3h-1)-2V_s t_*^2}{t^2}=0,
\qquad
\boxed{V_s t_*^2=\frac{F_\chi^2}{2}(3h-1)}.
\tag{12}
$$

正势要求 $h>1/3$。$h=1/3$ 要求零势，是另一个退化情况；$h<1/3$ 不属于这里的正势分支。于是

$$
\rho_\chi=\frac{3hF_\chi^2}{2t^2},\quad
p_\chi=\frac{(2-3h)F_\chi^2}{2t^2},\quad
w_\chi=\frac{2}{3h}-1.
\tag{13}
$$

现在加入独立守恒的单一流体 $p_f=w\rho_f$，要求两部分均非零且都按 $t^{-2}$ 缩放。$\rho_f\propto a^{-3(1+w)}$ 给出

$$
h=\frac{2}{3(1+w)},\quad -1<w<1,
\qquad w_\chi=w.
\tag{14}
$$

Friedmann 约束进一步给出

$$
\boxed{\Omega_\chi=\frac{\rho_\chi}{3M_{\rm Pl}^2H^2}
=\frac{F_\chi^2}{2hM_{\rm Pl}^2}
=\frac{3(1+w)F_\chi^2}{4M_{\rm Pl}^2}},
\tag{15}
$$

$$
\rho_f=\frac{3h^2M_{\rm Pl}^2-3hF_\chi^2/2}{t^2}.
\tag{16}
$$

因此 $0<F_\chi^2<2hM_{\rm Pl}^2$ 保证流体密度为正。式 (12)—(16) 同时满足标量方程、流体守恒、Friedmann 与 Raychaudhuri 方程，并非只在指定 $H$ 上代入的一条轨迹。

## 5. 吸引子条件、健康条件与全宇宙历史限制

规范场下势为 $V=V_s\exp(-\lambda\phi/M_{\rm Pl})$，其中 $\lambda=2M_{\rm Pl}/F_\chi$。标准的膨胀平直 FLRW、正势、常 $-1<w<1$、无额外交换项系统中，非零流体 scaling 分支在

$$
\lambda^2>3(1+w)
\quad\Longleftrightarrow\quad
F_\chi^2<\frac{4M_{\rm Pl}^2}{3(1+w)}
\tag{17}
$$

时是晚期吸引子，$\Omega_\chi=3(1+w)/\lambda^2$。等号为分支合并边界，不应直接套用严格吸引结论；$w=-1$ 和 $w=1$ 也需要单独讨论。这是已有相平面结果在本符号下的重写，未在本轮重证明全部非线性稳定性。[原始分析及线性稳定性附录](https://arxiv.org/pdf/gr-qc/9711068)。

必要的经典健康条件包括 $F_\chi^2>0$、$B>0$ 和限定有效域；规范标量与 Maxwell 主部没有由这些系数造成的负动能。该最小二阶模型的局部传播主部仍使用同一度规光锥。$V_{,\phi\phi}=4V/F_\chi^2>0$ 也可直接核对，但以上均不替代带物质源的稳定性、辐射修正、第五力、等效原理和实际宇宙学限制。

**不能把一个 scaling 解当作从辐射期到今天的精确时间律。** 辐射 $h=1/2$ 要求 $V_s t_*^2=F_\chi^2/4$；物质 $h=2/3$ 要求 $V_s t_*^2=F_\chi^2/2$。同一组固定 $V_s,F_\chi,t_*$ 不能同时满足两者，更不能把 $\Lambda$ 加速期当作常 $h$。在固定作用量中，必须跨时代积分 $\chi(t)$，允许过渡偏离；不同渐近阶段的形式常数 $t_*$ 可以不同，但那不是重新自由拟合每个时代响应的许可。

tracking 分支满足 $w_\chi=w$，所以在尘埃或辐射背景中自身不提供晚期加速，也没有由此解决 Hubble 张力。要保留已观测的加速背景，需要另有合适组分或修改模型，并重新求解。旧项目采用物质加 $\Lambda$ 的 $t(z)$ 再人为代入精确对数律，属于另一种规定背景的现象学假设，不是本作用量已经保证的预测。

## 6. 电磁响应与可检验关系

在统一低能定义和相同重整化约定下，局部缓变场可将光子动能规范化，电磁耦合为 $e_{\rm eff}=e_*/\sqrt B$，故该有效模型定义

$$
\frac{\alpha(\chi)}{\alpha_0}=\frac1{B(\chi)}
=1+\beta(\chi^{-2}-\chi_0^{-2}).
\tag{18}
$$

这里不把探测能标变化造成的 QED 运行当作宇宙漂移；忽略场在原子尺度的梯度也属于近似条件。对于完整背景解，直接预测应使用

$$
\delta_\alpha(t)=\beta[\chi(t)^{-2}-\chi_0^{-2}],
\qquad
D(t)=\frac{\dot\alpha}{\alpha}
=-\frac{2\beta\dot\chi}{\chi^3\mathcal A(\chi)}.
\tag{19}
$$

只有当同一观测区间符合式 (11) 时，才恢复 v0.1 的精确式，且 $\Gamma=\beta$：

$$
\delta_\alpha=\Gamma\left[\ln^{-2}(t/t_*)-\ln^{-2}(t_0/t_*)\right],
\quad D_0=-\frac{2\Gamma}{t_0\ln^3(t_0/t_*)}.
\tag{20}
$$

跨时代时仍能消去 $\beta$，条件是背景解的形状已固定且 $\dot\chi_0\ne0$：

$$
\delta_\alpha(t)=D_0\left\{-\frac{\chi_0^3}{2\dot\chi_0}
[\chi(t)^{-2}-\chi_0^{-2}]\right\}.
\tag{21}
$$

若 $\dot\chi_0=0$，这个消幅度式不适用；当前瞬时漂移可以为零而历史变化非零。若电磁源反作用使 $\chi(t)$ 本身依赖 $\beta$，式 (21) 仍是给定解上的恒等式，但不能再把花括号当成与待拟合幅度无关的固定传递函数。

规范场的无量纲电磁响应为

$$
q_\alpha(\phi)\equiv M_{\rm Pl}\frac{\partial\ln\alpha}{\partial\phi}
=-\frac{2\beta M_{\rm Pl}}{F_\chi\chi^3\mathcal A},
\quad
\frac D H=q_\alpha\frac{\dot\phi}{M_{\rm Pl}H}.
\tag{22}
$$

对正则均匀标量，$\dot\phi^2=3M_{\rm Pl}^2H^2\Omega_\chi(1+w_\chi)$；因此在向正方向滚动的分支

$$
\frac D H=q_\alpha\sqrt{3\Omega_\chi(1+w_\chi)}.
\tag{23}
$$

这是响应强度、背景能量与时间漂移之间的同一套参数关系，但 $q_\alpha$ 中仍有独立 $\beta$，不能声称仅由 $\Omega_\chi$ 就预测了 $D$。物质的电磁质量占比可使同一响应进入组成依赖的作用力；地面与天文系统共享均匀背景也必须受约束。[标量—电磁耦合的物理限制与自然性讨论](https://arxiv.org/pdf/hep-ph/0307299)。

$\beta=0$ 恢复恒定 $\alpha$ 并切断这项显式电磁交换，但 $\rho_\chi$ 一般仍非零。因此“常数零假设”在电磁观测中的嵌套成立，不等于完整宇宙背景自动退回无新增组分的标准模型；后者还需要控制标量能量。

## 7. 为什么平方仍不是唯一结论

保持同一个势和背景，将式 (2) 换成

$$
B_s(\chi)=\left[1+\beta_s(\chi^{-s}-\chi_0^{-s})\right]^{-1},
\qquad s>0,
\tag{24}
$$

在正值有效域同样得到 $\delta_\alpha=\beta_s[L^{-s}-L_0^{-s}]$。背景的指数势本身不选择 $s=2$。这个反例直接说明：把指定时间函数补为某个作用量，并不等于时间函数已获得唯一微观推导。相关的标量—电磁研究也明确指出，自由选择响应函数会削弱 $\alpha(t)$ 的预测性。[Copeland、Nunes 与 Pospelov，第 2—3 节](https://arxiv.org/pdf/hep-ph/0307299)。

本框架保留 $s=2$ 的合理表述是“由已发表的算术数值工作激发、预先指定的响应候选”。若要让“黎曼”成为不可替代的物理内容，仍须从已定义的算术状态或算符推出该响应，或者提出一般平滑响应模型做不到的独立关系。作用量 (1) 本身没有零点、素数序列或相应谱算符。

在给定 scaling 阶段，式 (12) 将 $t_*$ 与势的归一化联系起来：$t_*^2=F_\chi^2(3h-1)/(2V_s)$。由于 $V_s$ 尚无独立测定，这并不自动确定 $t_*$，更不证明它等于普朗克时间。固定响应函数后，$\chi=0$ 相对势的位置具有模型含义；不能一边保持式 (2) 不变，一边把任意场原点平移说成纯单位换算。

## 8. 本轮可以固定下来的主张

| 主张 | 状态与条件 |
|---|---|
| 存在具有正则动能与规范电磁耦合的经典有效实现 | 已给出，限定 $\chi>0,B>0,F_\chi>0$；不是 UV 完备性结论 |
| 子系统的非自治驱动可以来自动力学场 | 是，指定背景解后得到；完整系统保留反作用与能量交换 |
| 精确对数背景存在 | 是，常 $w$、忽略平均标量源、满足 scaling 参数条件 |
| $1/\ln^2$ 已从黎曼性质推导出来 | 否，平方仍在电磁响应中作为假设输入 |
| 能从同一个模型连接历史 $\alpha$ 与现在钟漂移 | 是，使用实际背景解与同一响应；式 (20) 的简式限于对数段 |
| 现有数据已发现该场或该响应 | 否，本轮没有新观测拟合，既有约束仍允许零电磁幅度 |
| 对数解能精确贯穿真实宇宙全部阶段 | 否，需跨辐射、物质、加速阶段积分并重做约束 |

下一步最有区分度的计算，是固定一组正则参数后求同一作用量的多组分宇宙背景，量化式 (20) 与式 (19) 的差异，再评估是否能使用旧钟—高红移联合似然。不能先用旧曲线完成拟合，再把“存在一个 scaling 解”当作完整动力学验证。
