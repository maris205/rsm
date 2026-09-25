# 明确 no-scale 基准：确实取消隐藏常数传递，但还没有稳定宇宙时钟

本报告只研究下述未经稳定的完整超引力作用量。结论分为两层：no-scale 几何确实从势和局部带电谱中消除了上一轮危险的 $W_c$ 依赖；同一作用量却留下自由模量，不能把它固定为常数后便宣布获得了完整的慢时钟宇宙解。这里没有新增观测拟合，也没有完成模量稳定或全超引力量子计算。

## Material Passport

- 日期：2026-09-25。任务：独立推导这一明确基准，先完成解析与核验，不阅读其他本轮分支的结论。
- 复用材料：`sugra_shift_v01` 的参数与 `shift_axis` 的 1025 个 $\chi$ 坐标。旧轨迹只是固定截面的采样点，不作为新模型的解。
- 新材料：[可执行核验](noscale_exact_20260925_checks.py)、[结果与输入哈希](noscale_exact_20260925_checks.json)。本轮主报告只改这些指定文件。
- 主来源阅读范围：Ellis 等 2015 年论文的摘要／元数据，以及 Ellis 等 2017 年论文第 2.2 节式 (4)–(8) 与邻接段落。本文的 shift-clock 推导与数值是本轮计算，不能归给来源作者，也不声称已有物理实验证据。
- 内部 AI 数学与程序审计，不等同外部同行评议；不存在统计抽样或 $p$ 值。

## 1. 作用、量纲与完整抵消

记 $P=M_p^2$，$T$ 无量纲，$Z,Q_\pm$ 具有质量量纲，$K$ 具有质量平方量纲：

$$
K=-3P\ln Y,\qquad
Y=T+\bar T-\frac{k(Z,\bar Z)+|Q_+|^2+|Q_-|^2}{3P},
\qquad k=-\frac12(Z-\bar Z)^2,
$$
$$
W=W_c+w(Z)+M(Z)Q_+Q_-,\quad
W_c=m_G P e^{i\theta},\quad W_T=0.
$$

要求 $Y>0$。定义 $\varphi^a=(Z,Q_+,Q_-)$ 及 $k_{\mathrm{all}}=k+|Q_+|^2+|Q_-|^2$，有 $k_{a\bar b}=\delta_{a\bar b}$。完整正定线元可写为

$$
ds_K^2=\frac{3P}{Y^2}\left|dT-\frac{k_a d\varphi^a}{3P}\right|^2
+\frac{1}{Y}\sum_a |d\varphi^a|^2.
$$

必须使用包括 $T$ 的整个逆度规。由直接矩阵求逆可得

$$
K^{I\bar J}D_IW D_{\bar J}\bar W
=\frac{Y^2}{3P}|D_TW|^2+
Y\sum_a\left|D_aW+\frac{k_a}{3P}D_TW\right|^2
=\frac{3|W|^2}{P}+Y\sum_a|W_a|^2.
$$

其中 $D_TW=-3W/Y$。因此整个 F 势精确成为

$$
\boxed{V_F=\frac{|w_Z+M_ZQ_+Q_-|^2+
|M|^2(|Q_+|^2+|Q_-|^2)}{Y^2}.}
$$

这不是在实轴展开后恰巧抵消：对于允许域内的轴外 $Z$ 与非零 $Q_\pm$，常数 $W_c$ 都不出现在此 F 势中。若规范动能为常数，未破缺 U(1) 的 D 势从 $Q^4$ 开始，不改本报告的二次谱。该限定不替代完整规范／引力圈图。

标准 no-scale 文献已有对数 Kähler 几何和某些零 soft 边界的研究；文献也明确区分了模量取值假设和动力学固定问题。因此本结果是把已知结构应用并核验于指定 clock 模型，不是新发现 no-scale 抵消。[No-scale SU$5$ super-GUTs，第 2.2 节](https://link.springer.com/article/10.1140/epjc/s10052-017-4805-x)；[Phenomenological Aspects of No-Scale Inflation Models，摘要](https://arxiv.org/abs/1503.08867)。

在本报告采用的辅助场符号约定中，$F^I=-e^{K/(2P)}K^{I\bar J}D_{\bar J}\bar W$。于是

$$
F^a=-Y^{-1/2}\overline{W_a},\qquad
F^T=\frac{Y^{-1/2}}{P}\left(\bar W-\frac13\sum_a k_a\overline{W_a}\right).
$$

实轴 $Q=0,y=0$ 上 $F^T=\bar W/(P\sqrt t)$，其中 $t=T+\bar T$。它一般不为零，删除它正好会删除实现抵消的部门。物理引力微子质量为

$$
m_{3/2}=\frac{|W_c+w|}{P t^{3/2}},
$$

不应把参数 $m_G$ 无条件当作 $m_{3/2}$。本模型没有 nilpotent 约束；上一轮 nilpotent 场的 VA 截止标记不能直接搬来作为这一线性模量作用量的截止定理。这里也未提供 UV 完成。

## 2. 轴外势与真正的带电谱

沿用

$$
Z=\frac{F_\chi}{\sqrt2}(\chi+iy),\qquad
w=-\frac{F_\chi\sqrt{U_i}}{\sqrt2}
e^{-\sqrt2(Z-Z_i)/F_\chi},\qquad
U=U_i e^{-2(\chi-\chi_i)},\quad
\varepsilon=F_\chi^2/P.
$$

在 $Q=0$ 时，$Y=t-\varepsilon y^2/3$，有精确结果

$$
V=\frac{U}{(t-\varepsilon y^2/3)^2},\qquad
V_y|_{y=0}=0.
$$

这一基准本身没有上一轮的 $\rho_\Lambda$ 常数部门，也没有 $b=1-3\varepsilon/2$ 因子。添加宇宙常数或 uplift 是另一项作用假设，必须重新分析其模量依赖，不能默认已经包含。

在任意固定背景点的局部无导数二次谱中，带电场的规范化为 $q_\pm=Q_\pm/\sqrt Y$，其费米子质量为 $M_{\mathrm{phys}}=M/\sqrt Y$。实轴上记

$$
s=\frac{|M|^2}{t},\qquad
d=\frac{2U}{3Pt^2}=\frac{2V}{3P},\qquad
\mathcal B=\frac{\overline{w_Z}M_Z}{t},
$$
$$
\boxed{m_\pm^2=s+d\pm|\mathcal B|,\qquad m_\psi^2=s.}
$$

抵消掉的是 $m_G^2$ 型共同分裂和 $m_GM$ 型 B 项；**不是全部分裂都为零**。剩余 $d$ 来自 $Y^{-2}$ 的展开，剩余 B 项来自 clock 自身的 $w_ZM_Z$。只要 clock 有非零势，笼统使用“所有 soft 参数为零”就会漏项。

沿实轴 $M_Z=\sqrt2 M_\chi/F_\chi$。若 $M^2=m_\infty^2(1+\xi/\chi^2)$，

$$
r_M\equiv \frac{M_\chi}{M}=-\frac{\xi}{\chi(\chi^2+\xi)},\qquad
|\mathcal B|^2=\frac{2U|M|^2r_M^2}{F_\chi^2t^2}.
$$

这些是局部势的谱，不能混同于任意弯曲、滚动宇宙背景中的完整传播频率；背景导数和曲率修正仍需另外处理。

## 3. 模量是动力学变量：静态 runaway 与动态例外

在 $y=0$、$\operatorname{Im}T$ 固定的子空间，动能为

$$
\mathcal L_{\mathrm{kin}}=-\frac{3P}{4t^2}(\partial t)^2
-\frac{F_\chi^2}{2t}(\partial\chi)^2.
$$

令

$$
\sigma=\sqrt{\frac32}M_p\ln t,\qquad
a=\sqrt{\frac23}\frac1{M_p},
$$

得到

$$
\mathcal L_{\mathrm{kin}}=-\frac12(\partial\sigma)^2
-\frac12F_\chi^2 e^{-a\sigma}(\partial\chi)^2,
\qquad V=U_i e^{-2(\chi-\chi_i)-2a\sigma}.
$$

因此 clock 的局部规范化也随模量改变。不能在任意 $t(t_{\mathrm{cos}})$ 下继续把 $F_\chi\chi$ 当成独立的全局正则单场。有限 $U>0$ 时

$$
V_\sigma=-\sqrt{\frac83}\frac V{M_p}\ne0,\quad
V_{\sigma\sigma}=\frac{8V}{3P}.
$$

势沿 $t\to\infty$ 下降，正则距离无限；该方向没有静态极小值。正的二阶导数并不是“已经稳定的模量质量”。若 clock 势主导 Friedmann 方程，则此曲率尺度仅 $8H^2$，也不构成可被无条件积分掉的高质量部门。

实轴的静态横向势曲率换算为 $m_{y,\mathrm{pot}}^2=4V/(3P)>0$。这同样只是势部分；滚动多场扰动还包含场空间几何与速度项。Im $T$ 没有势，未得到稳定。

一个必须保留的例外是：没有静态极小并不等于所有固定 $t$ 的运动均不可能。齐次 FLRW 模量方程是

$$
\ddot\sigma+3H\dot\sigma+
\frac a2 G_{\chi\chi}\dot\chi^2-2aV=0,
\qquad G_{\chi\chi}=F_\chi^2/t.
$$

故固定 $t$ 的必要条件为 $K_\chi=G_{\chi\chi}\dot\chi^2/2=2V$。结合 clock 方程，可得一支辐射跟踪解：

$$
\dot\chi=2H,\quad \dot H=-2H^2,\quad
V=F_\chi^2H^2/t,\quad \Omega_\chi=\varepsilon/t,\quad w_\chi=1/3.
$$

若 $t\ge\varepsilon$，余下能量可由辐射提供。这是特殊动态平衡，不是模量静态稳定，也未在此证明吸引域或完整扰动稳定。它不能说明此前晚期 dust+$\Lambda$ 轨迹同样允许 $t=1$。本轮应说“旧固定模量晚期背景尚未得到作用量支持”，而非宣称一切 no-scale 宇宙演化均失败。

## 4. 理想基准的带电一圈检查

固定局部背景和物理匹配标尺 $\mu=100\,\mathrm{GeV}$，有限 CW 核为

$$
V_1=\frac{f(s+d+h)+f(s+d-h)-2f(s)}{32\pi^2},\qquad
f(x)=x^2\left[\ln(x/\mu^2)-\frac32\right],\quad h=|\mathcal B|.
$$

在 $|d|/s,h/s\ll1$ 时，

$$
V_1=\frac{4sd(L-1)+2h^2L}{32\pi^2}+\cdots,
\qquad L=\ln(s/\mu^2).
$$

这里保留了共同 $d$ 和 B 分裂、及其完整 $\chi$ 导数；没有靠删掉斜率或改有限反项压低反馈。所有 $m_G$ 依赖先在作用量级消失，因而此局部带电核也不含 $m_G$。残余项受 $M^2/P$ 和 $M_\chi^2/F_\chi^2$ 等因子控制。

用旧参数、固定 $t=1$ 的 1025 点截面得到

| 局部诊断 | 数值 |
|---|---:|
| $\max \lvert V_{1,\chi}\rvert/(2U)$ | $1.4271135835\times10^{-35}$ |
| $\max d/s$ | $6.9281795078\times10^{-91}$ |
| $\max \lvert\mathcal B\rvert/s$ | $5.2370687672\times10^{-46}$ |

这表明**这一明确理想 no-scale 基准，确实消除了上一轮的那一个带电阈值障碍**。但它尚未证明稳定完成中的同样抵消、全量子有效作用的稳定性或观测可行性。计算未包含模量稳定部门、引力微子／引力圈、规范 anomaly、高圈效应及完整弯曲时空有效作用。在宇宙背景下，曲率项可能改变这一极小局部核的系数；本数值不应作为完整量子误差上界。

没有据此制作新的 $\alpha(t)$ 曲线：因 $t$ 尚未固定，物理阈值 $M/\sqrt t$ 和 clock 轨迹均可能改变。拟议的 $1/\ln^2 t_{\mathrm{cos}}$ 慢变化仍是需要通过一致背景实现的假设。

## 5. 补充：理想中心的 clock 费米子是否有 $m_G$ 阶质量

在 $w\to0$、$Q=y=0$ 的 Minkowski no-scale 点，$F^T\ne0,F^Z=0$。直接求导得到

$$
K_{ZZ}=-1/t,\quad
\Gamma^T_{ZZ}=1/(3P),\quad
D_TW=-3W/t.
$$

故

$$
D_ZD_ZW=K_{ZZ}W/P-\Gamma^T_{ZZ}D_TW=0.
$$

由于 $D_ZW=0$，去掉 goldstino 方向的费米子质量修正不改变此式：幸存的 clock Weyl 费米子质量为零，不能无依据添加一个 $m_G$ 阶 clock 费米子抵消某个模量标量圈。$T$ 费米子是被引力微子吸收的 goldstino，不再是独立传播的自旋 $1/2$ 粒子。此点 $m_{3/2}=|W_c|/(Pt^{3/2})$。

若某项 Kähler 修正在中心保留至三阶导数，它会继承此中心费米子结论；其四阶导数可以改变标量谱。这只是中心的局部说明，不证明修正后的势、移位谷底或量子超迹如何。在 $w_Z\ne0$ 时 goldstino 方向发生混合，也不能照搬这个理想极限。

## 6. 核验、问题日志与结论范围

最终 **52/52** 个检查通过，包括：四场度规的精确逆矩阵；任意常数超势的完整抵消；轴外复数直接 K/W 势；二次带电谱；正则模量和辐射跟踪恒等式；理想点费米子投影；3 个 $\chi$ 点与 $t=0.5,1,2$ 的 190 位直接 supertrace 及其导数。直接 CW 与稳定展开的比较门限为相对 $10^{-70}$。这些是数学／实现核验，而不是模型获证。

首轮为 45/46：一般 $F^T$ 期望式误把 $k_a$ 写成其共轭，直接矩阵求逆指出并纠正了指标。错误不涉及轴上 $k_Z=0$ 的结果、标量势、带电谱或 CW。新增费米子检查后得到上述 52/52。工具日志还记录首次 `rg` 正则转义错误，以及 arXiv 旧文 HTML 返回 406；改用有限文件读取和官方摘要／期刊正文，没有隐藏失败或伪造读完全文。

来源响应 SHA-256 与读取范围保存在结果 JSON。期刊下载经过 cookie 查询重定向，哈希仅标识当次返回 HTML；公式另由 web 正文读取核对。未声称对整篇论文逐页精读。

这一步支持继续研究保护机制的理由是明确的：有一个完整作用使危险项解析取消，且局部数值复算符合解析。但完成物理模型的条件也同样明确：必须稳定或控制模量，保留相应辅助场，并重新计算新增部门的势、谱和阈值。不能把这些条件压缩成一句“固定 $T=1$”略过。
