# 四次 Kähler 稳定化的独立推导核查

日期：2026-09-25。状态：候选有效作用量的局部解析核查；不是新的宇宙学拟合、标准模型构造或完整超引力量子有效作用量。

这个候选与超势稳定化不同：保留完整的模量多重态，以 Kähler 四次项给模量质量，并让中心的低阶几何导数保留 no-scale 抵消。结果是**局部树级可以同时稳定模量、消除常数超势导致的主要时钟交叉项，并抑制带电软分裂**。但新重模量自身的一圈反馈仍需计算；保护没有因此自动扩展到全部量子修正。

## 1. 指定模型与精确局部势

取无量纲复模量 $T$，有量纲时钟 $Z=F_\chi(\chi+iy)/\sqrt2$，以及带电场 $Q_\pm$：

$$
K=-3M_P^2\ln\Omega,\qquad
\Omega=f(T,\bar T)-\frac{k(Z,\bar Z)+|Q_+|^2+|Q_-|^2}{3M_P^2},
$$
$$
k=-\frac12(Z-\bar Z)^2,\qquad
W=w(Z)+W_c+M(Z)Q_+Q_-,\quad W_T=0,
$$
$$
t=T+\bar T,\quad v=2\operatorname{Im}T,\quad \delta=t-1,
\qquad f=t+a\delta^4+bv^4.
$$

这里 $a,b>0$ 是局部稳定所需的符号。令

$$
p=f_T=1+4a\delta^3-4ibv^3,\quad
q=f_{T\bar T}=12a\delta^2+12bv^2,
$$
$$
E=|p|^2-q\left[\Omega+\frac{|k_Z|^2}{3M_P^2}\right].
$$

从完整 $(T,Z)$ Kähler 度规求逆、代入超引力 F 势得到，在 $Q_\pm=0$ 时

$$
\boxed{V_F=\frac{|w_Z|^2}{\Omega^2}
+\frac{q\,|3W-k_{\bar Z}w_Z|^2}{3M_P^2\Omega^2E}.}
$$

这包含复 $T$、复 $Z$ 的局部几何，不是仅在轴上猜出的势。健康局部域需要 $\Omega>0,E>0$。实 $Z$ 轴上 $k=k_Z=0$，定义 $D=|p|^2-fq$，则

$$
\boxed{V_F=\frac{U}{f^2}+\frac{3q|W|^2}{M_P^2f^2D}},
\qquad U=|w_Z|^2.
$$

中心 $t=1,v=0$ 的 $f=p=D=1,q=0$，所以 $V_F=U$，任意常数 $W_c$ 均消去。四次形变的前一至三阶导数在中心不改变原几何；四阶导数给出模量质量。这里仍沿用原有指数 $w$ 与逆对数质量函数，**没有从黎曼结构推导真实相互作用或证明 $1/\ln^2t$ 必然出现**。

## 2. 质量、补偿子与真实带电谱

先取 clock-off：$w=w_Z=0$，且 $m_G=|W_c|/M_P^2$。模量中心是局部 Minkowski 极小值，规范化实自由度为 $\sigma_t=\sqrt{3/2}M_P\delta$、$\sigma_v=\sqrt{3/2}M_Pv$，质量

$$
m_t^2=48a m_G^2,\qquad m_v^2=48b m_G^2.
$$

该质量系数经过独立符号求导与规范化核查，也与 2019 年原始论文的独立实、虚模量质量吻合。$F^T\ne0$，模量费米子提供被 gravitino 吸收的 goldstino；不能再把它作为一个独立同质量费米子去抵消两个模量标量的 Coleman–Weinberg 项。

采用 $F^I=-e^{K/(2M_P^2)}K^{I\bar J}D_{\bar J}\bar W$ 的符号，在实 $Z$ 轴定义复量 $\mathfrak m=\bar W/(M_P^2 f^{3/2})$，则

$$
\frac{F^C}{C}=\mathfrak m+\frac{K_IF^I}{3M_P^2}
=-\mathfrak m\frac{fq}{D}.
$$

中心 $q=0$ 时此补偿子组合严格为零；物理 $m_{3/2}=|\mathfrak m|$ 则可以非零。该组合的抵消本身不能替代完整、规范一致的 anomaly-mediated gaugino 计算。

展开完整带电势，并除以带电动能系数 $Y=1/f=e^{K_0/(3M_P^2)}$，得到

$$
m_\psi^2=\frac{|M|^2}{f},\qquad
m_\pm^2=m_\psi^2+d\pm|\mathcal B|,
$$
$$
\boxed{d=\frac{2V_F}{3M_P^2}},\qquad
\boxed{\mathcal B=\frac{\bar w_Z M_Z}{f}
+\frac{qM\bar W}{M_P^2fD}}.
$$

其中 $d$ 的恒等式也可由 $\partial_I\partial_{\bar J}\ln Y=K_{I\bar J}/(3M_P^2)$ 与 F 势恒等式得到。clock-off 中心 $d=\mathcal B=0$；它没有上一轮加和式隐藏部门的 $d\simeq m_G^2$、$|\mathcal B|\simeq m_G|M|$。这里的 $V_F$ 只包含指定作用量的势；若人为另加常数，不能把它直接塞进这个超引力软质量恒等式。

## 3. 小时钟驱动引起的模量移动

固定一瞬间 $\chi$，且 $U\ll m_G^2M_P^2$，中心仍有 $V_t=-2U$。真实局部谷底因此略微偏移：

$$
\delta_*\simeq\frac{UM_P^2}{36a|W|^2},\quad v_*=0,
\qquad q_*\simeq\frac{U^2M_P^4}{108a|W|^4}.
$$

势的局部修正为

$$
V_{\rm eff}\simeq U-\frac{U^2M_P^2}{36a|W|^2}.
$$

故带电 $\mathcal B$ 的额外 $W_c$ 项量级为 $M U^2/(a m_G^3M_P^4)$，而不是重新回到 $M m_G$。这些式子要求谷底位移小、模量质量大于变化率，并处于正度规域；不能用于 $m_G\to0$ 或近 $W_c+w=0$ 的区域。

clock-on 的几何中心并非完整驻点，因此该处普通坐标二阶导数的 $O(U/M_P^2)$ 修正不应直接叫作物理质量。更精确的质量需在谷底结合协变 Hessian、动能度规及轻重混合求本征值；本报告环压力测试只使用 $m_G\gg H$ 时的重质量主项。

时钟虚部并未因此变重。无二次形变时，中心

$$
m_y^2=\frac{4U}{3M_P^2},
$$

只有宇宙演化量级。后续量子受力必须检查 $(\chi,y)$ 两个方向；只检查实轴斜率会漏掉特殊相位的横向反馈。

## 4. 小二次 Kähler 形变容纳树级正真空能

可进一步明确更改作用量

$$
f=t+\eta(t-1)^2+a(t-1)^4+bv^4,
\quad p=1+2\eta\delta+4a\delta^3-4ibv^3,
\quad q=2\eta+12a\delta^2+12bv^2.
$$

以上势、谱、补偿子公式继续成立。中心 clock-off 势为

$$
V_0=\frac{6\eta|W_c|^2}{M_P^2(1-2\eta)},
\quad
\eta\simeq\frac{\rho_\Lambda}{6m_G^2M_P^2}>0.
$$

这只是一次指定的真空能系数调节，并未自然解释宇宙学常数。真实极小值需要一起求解，leading 为

$$
\delta_*\simeq\frac{\eta}{6a}+\frac{UM_P^2}{36a|W|^2},
\qquad
\mathcal B_\eta\simeq\frac{2\eta M\bar W}{M_P^2}.
$$

最小值相对中心真空能差是 $O(\eta^2m_G^2M_P^2/a)$，应在精密数值中调节谷底值，而不把中心近似当作精确匹配。对于 $m_G\gg H$，由暗能量产生的额外 B 项仅为 $O[M\rho_\Lambda/(m_GM_P^2)]$。

横向时钟也可以显式检查。令 $w_r=-F_\chi\sqrt U/\sqrt2$ 为实轴值，$W_c=m_GM_P^2e^{i\theta}$、$\epsilon=F_\chi^2/M_P^2$，则中心

$$
3W-k_{\bar Z}w_Z=3W_c+w_r(3+2iy)e^{-iy}.
$$

以 $q=2\eta,D=1-q$ 记，

$$
V_y=-\frac{2q m_Gw_r\sin\theta}{D},
$$
$$
V_{yy}=\frac{4\epsilon U}{3}+\frac{q}{3M_P^2D}
\left[8w_r^2+6\operatorname{Re}(\bar W_cw_r)
+6\epsilon\left(2+\frac qD\right)|W|^2\right].
$$

当 $\eta$ 按上述正真空能选择、$m_G\gg H$ 时，主项给 $m_y^2\simeq4(U+\rho_\Lambda)/(3M_P^2)>0$，横向偏移受 $H/m_G$ 抑制。这是局部层级展开，不是任意参数的全局稳定定理。

## 5. 新模量自身的量子反馈与真空能减除

实轴中心的两个重模量质量由 $|W|^2$ 给出，但**横向求导必须先回到完整势**。定义

$$
X(\chi,y)=\frac{|S(\chi,y)|^2}{9M_P^4},\quad
S=3W_c+w_r(3+2iy)e^{-iy}.
$$

实轴上 $X=|W_c+w_r|^2/M_P^4$。忽略小谷底位移与 $U$ 级轻重混合时，两个实模量的重质量主项为 $c_tX,c_vX$，$c_t=48a,c_v=48b$。完整势的三阶导数核对给出

$$
\partial_y m_t^2\big|_0=-\frac{32a m_Gw_r\sin\theta}{M_P^2},\qquad
\partial_y m_v^2\big|_0=-\frac{32b m_Gw_r\sin\theta}{M_P^2}.
$$

它们是错误地直接对 $48a|W|^2/M_P^4$ 横向求导的三分之一。$X$ 的使用在这里是局部重质量展开，不宣称任意横向位移处全谱都可由一个函数描述。仅这两个标量的平直时空有限 Coleman–Weinberg 诊断是

$$
V_{1,T}=\sum_{c=c_t,c_v}\frac{c^2X^2}{64\pi^2}
\left[\ln\frac{cX}{\mu^2}-\frac32\right],
$$
$$
\partial_A V_{1,T}=\sum_c\frac{c^2X\,\partial_AX}{32\pi^2}
\left[\ln\frac{cX}{\mu^2}-1\right],\quad A=\chi,y.
$$

一般相位的纵向主项为 $O[(a^2+b^2)m_G^3w_\chi/M_P^2]$。$\theta=\pi/2$ 消去纵向线性交叉，但横向交叉最大：

$$
X_\chi\big|_0=\frac{2\operatorname{Re}(\bar Ww_\chi)}{M_P^4},\qquad
X_y\big|_0=-\frac{2m_Gw_r\sin\theta}{3M_P^2}.
$$

所以不能凭特殊相位宣布量子稳定；应比较二维受力，再必要时重新求横向谷底与宇宙背景。

上述表达式只是一组已指定自由度的局部压力测试。它不包含完整 gravitino、引力规范固定、鬼场及其混合，也没有计算完整重整化群改正或证明所有修正相加必为这个值。反过来，也不能在没有完整计算时假定 gravitino 会精确取消它。

**真空能调节与减除方案必须一致。** 如果只是减掉 $V_{1,T}(\chi_i)$ 常数，所得是指定有限势差的诊断；若要通过同一个 Kähler 参数 $\eta$ 在作用量内取消 clock-off 环真空能，则必须相应重调

$$
\eta_{\rm rad}\simeq-\frac{V_{1,T}(X_0)}{6M_P^2X_0}.
$$

它同时改变 $q$、带电 B 与树级斜率。leading 总环斜率相应变为

$$
\left[\frac{dV_{1,T}}{dX}-\frac{V_{1,T}(X_0)}{X_0}\right]\partial_AX.
$$

在 $X=X_0$，每个标量的系数是 $c^2X_0[\ln(cX_0/\mu^2)-1/2]/(64\pi^2)$，仍然有 $m_G^3w/M_P^2$ 量级。这个区别可以改变预算的数值系数，故不能把树级暗能量匹配与简单环常数减除合并宣称为“完成了一圈 de Sitter 真空”。

另一个配套影响是 $\eta_{\rm rad}$ 使 $V_{\rm tree}$ 本身含补偿环真空能的项，故实际树级 $d=2V_{\rm tree}/(3M_P^2)$ 可达 $m_G^4/M_P^2$ 量级；不能用已经互相抵消后的总势代换该恒等式。将这个 $d$ 再送入带电阈值是有用的配套压力测试，但在环计数上已涉及更高阶。完整模量一圈有效作用量对 $Q^2$ 的直接修正、以及同阶其余贡献未在本解析核查求出，因此配套预算也不能被表述为完整两圈约束。

## 6. 有效理论范围与文献位置

该四次稳定化是已有的有效作用量构造。2013 年原始论文 Sec. 4、PDF Eq. (63) 给出正号四次项，并将其作为可行性构造而非已导出的微观机制；2019 年论文 Eq. (2)、Eq. (8) 后质量段给出独立实、虚模量质量及 $\Lambda_T>m_T$ 的有效理论一致性条件。[2013 原始论文](https://arxiv.org/pdf/1307.3537)，[2019 原始论文](https://arxiv.org/html/1903.05267v2)。注意 2013 HTML 会将同一式重编号为 Eq. (67)，因此上述定位使用 PDF 编号。

若四次项由尺度 $\Lambda_T$ 的局部算符产生，数量级可记 $a,b\sim M_P^2/\Lambda_T^2$。需有 $m_t,m_v\ll\Lambda_T$，且若要计算 100 GeV 阈值，还需 $100\,\mathrm{GeV}\ll\Lambda_T$。例如 $a,b=0.1,1,10$ 的几何截断标度接近 $M_P$，在此仅给参数计数；它不是发现了真实 UV 完成。系数小于 1 时不能把估算截止推到量子引力尺度以上。大系数提升模量质量，同时降低几何截止，不能无限放大。

局部 $D>0$ 域也有限；该多项式不能被直接当作任意大模量位移的健康理论。保留完整线性模量多重态后，上一轮 nilpotent 有效理论的低截止判断不能照搬，但仍需说明新几何 EFT 自己的适用范围。没有在本核查中加入额外 nilpotent 场或外部上抬常数。

## 材料与复现记录

- 数学输入：本报告给出的 K、W；原项目的指数时钟超势形式；没有外部观测数据或统计拟合。
- 原始来源：上述两篇论文用于四次稳定化构造和范围核对。没有把其宇宙模型、退相干结论或完整 UV 完成当作本模型已证事实。1984 年历史前驱只由 2013 文献追溯，原文正文未独立核读，故不据其引用具体式子。
- [独立检查脚本](noscale_kahler_20260925_checks.py) 与 [机器可读结果](noscale_kahler_20260925_checks.json)：35 个符号恒等式；脚本不导入主实验代码。涵盖复几何势、模量质量、带电展开、谷底 leading、正真空能形变、横向时钟导数、完整势重模量横向质量导数及 CW/真空能反项链式求导。
- 范围：上述局部解析检查；主实验的高精度谷底与预算扫描由另外的实现完成。本报告不重复使用其结果作为独立核验。
- 执行记录：单次检查进程均设 180 秒超时；阶段性 27/27、31/31、32/32 通过，最终补充三个重模量横向导数核验。内部审查发现早期手推曾将实轴 $m_T^2\propto|W|^2$ 直接用于横向导数，遗漏 $S(y)$ 几何因子；在正式扫描前已修正为三分之一系数并记录相应独立检查。没有符号检查失败或科学结果丢弃。
- 透明度：AI 内部分工推导和检查，同模型家族的独立实现；不等于外部同行评审。用户未作原文阅读声明。
