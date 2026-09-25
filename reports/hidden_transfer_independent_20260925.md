# 隐藏部门传递：独立推导与适用范围核查

## Material Passport

- 日期：2026-09-25。
- 任务：独立核查指定 nilpotent 隐藏部门的两实场势、轴向力、带电谱及尺度预算。
- 材料：父实验 `experiments/sugra_shift_v01/results/inputs.json` 与 `trajectories.csv`；任务给出的新 $K,W$。未读取本轮另一理论审查者文件，也未导入主实验实现。
- 方法：从超引力势直接展开；以实部、虚部独立重建势，并从带电场的二阶展开取质量；24 项 SymPy 恒等式检查。补充读取一手理论文献，独立资料检索者未访问本项目。
- 产物：[检查脚本](hidden_transfer_independent_20260925_checks.py)、[24/24 检查和父材料哈希](hidden_transfer_independent_20260925_checks.json)。命令：`timeout 120 python reports/hidden_transfer_independent_20260925_checks.py`。
- 边界：这是同一模型家族代理的独立实现核查，非外部同行评审。未重算宇宙轨迹，未拟合观测；以下数值预算冻结在旧移位型轨迹上。

## 1. 约定与完整势

令 $F=F_\chi$、$P=M_{\rm Pl}$、

$$
 Z=\frac{F}{\sqrt2}(\chi+iy),\quad
 \epsilon=F^2/P^2,\quad
 U=U_i e^{-2(\chi-\chi_i)},\quad
 W_\chi=-\frac{F\sqrt U}{\sqrt2}e^{-iy}.
$$

指定候选为

$$
K=-\tfrac12(Z-\bar Z)^2+X\bar X+|Q_+|^2+|Q_-|^2,
\qquad X^2=0,
$$

$$
W=W_\chi+W_c+fX+M(Z)Q_+Q_-,\qquad
W_c=m_G P^2e^{i\theta},\qquad
|f|^2=3m_G^2P^2+\rho_\Lambda.
$$

$m_G\ge0$ 是常数超势定义的质量参数。完整模型的实际引力微子质量为

$$
m_{3/2}=e^{K/(2P^2)}|W_\chi+W_c|/P^2,
$$

因此 $m_G=0$ 并不表示实际 $m_{3/2}=0$。$\rho_\Lambda$ 已由隐藏部门生成，Friedmann 方程不得另加同一份外部常数密度。

在 $X=Q_\pm=0$ 上，标准 F 项势给出

$$
\boxed{V=e^{\epsilon y^2}\left\{
\rho_\Lambda+U\left(1-\frac32\epsilon+\epsilon^2y^2\right)
+2F^2m_G^2y^2+C\left[(3-2\epsilon y^2)\cos(y+\theta)+2y\sin(y+\theta)\right]
\right\}},
$$

其中 $C=\sqrt2 Fm_G\sqrt U$。这由以下表达式直接展开，而非假设一个额外软质量：

$$
D_ZW=\sqrt U(1+i\epsilon y)e^{-iy}-i\sqrt2Fym_Ge^{i\theta},
\quad D_XW=f,
$$

$$
V=e^{K/P^2}\left(|D_ZW|^2+|f|^2-3|W_\chi+W_c|^2/P^2\right).
$$

Nilpotent 的基本势与非零辅助场条件可分别对照 Bergshoeff 等的式 (1.3)、(3.1)。这里添加滚动 $Z$ 后的交叉项是本报告推导，不能归属于该文献的既有宇宙模型。[Pure de Sitter Supergravity](https://arxiv.org/pdf/1507.08264)

## 2. 三个相位及重横向场

记 $c=1-3\epsilon/2$。轴上导数为

$$
V=\rho_\Lambda+cU+3C\cos\theta,\quad
V_{,\chi}=-2cU-3C\cos\theta,\quad
V_{,y}=-C\sin\theta,
$$

$$
V_{,\chi\chi}=4cU+3C\cos\theta,\quad
V_{,\chi y}=C\sin\theta,
$$

$$
V_{,yy}=2\epsilon\rho_\Lambda+2\epsilon(1-\epsilon/2)U
+4F^2m_G^2+(1+2\epsilon)C\cos\theta.
$$

物理质量矩阵为以上 Hessian 除以 $F^2$。

| 相位 | 轴上新增 $\chi$ 力 | 轴上 $y$ 梯度 | 解释 |
|---|---:|---:|---|
| $0$ | $-3C$ | $0$ | 推动原方向滚动；大质量时不再是原慢轨迹的小修正 |
| $\pi/2$ | $0$ | $-C$ | 实轴不是静止谷底；不能只检查 $\chi$ 方程 |
| $\pi$ | $+3C$ | $0$ | 可反转滚动；负的交叉势也可能使原初始数据无正的 Friedmann 解 |

即便 $m_G=0$，新模型也有 $e^{\epsilon y^2}\rho_\Lambda$。它在轴上不改旧背景，却提供额外的横向质量

$$
\Delta m_y^2=2\rho_\Lambda/P^2=6\Omega_\Lambda H_{\rm ref}^2.
$$

因此不能直接沿用上一轮“横向始终很轻”的判断。

对 $\theta=\pi/2$ 和 $m_G/H\to\infty$，先取

$$
 y_*\simeq\frac{\sqrt U}{2\sqrt2Fm_G}
$$

再消去重场，得到

$$
\boxed{V_{\rm eff}\simeq\rho_\Lambda+(c-1/4)U}.
$$

这一结果也通过在无量纲势中固定 $v=(m_G/H_{\rm ref})y$、取大质量极限再求极小值验证。横向位移虽随 $1/m_G$ 变小，残留势修正仍为 $-U/4$。所以“$y$ 很小”不等于隐藏部门解耦；但也不能因实轴上有大的 $y$ 力，就排除重场谷底的慢运动。若仍把初态放在 $y=0$，它与谷底之间有约 $U/4$ 的势能差，需要演化或平均快振荡后才能说明宇宙历史。本报告未做该演化。

## 3. 带电粒子的质量及隐藏部门传递

轴上 $K=K_Z=0$，取正的实 $M(\chi)$，定义 $r=M_{,\chi}/M$。带电 Dirac 费米子与两复标量满足

$$
s=m_\psi^2=M^2,\qquad m_\pm^2=s+d\pm|B|,
$$

$$
\boxed{d=m_G^2+\frac{\rho_\Lambda}{P^2}
+(1-\epsilon)\frac{U}{P^2}
+\frac{2C\cos\theta}{P^2}},
$$

$$
\boxed{B=\frac{\sqrt{2U}M}{F}\left(r+\frac\epsilon2\right)
-m_GM e^{-i\theta}}.
$$

一般复背景可写为

$$
s=e^{K/P^2}|M|^2,\quad
d=m_{3/2}^2+V/P^2,
$$

$$
B=e^{K/P^2}\left[(D_ZW)^*\left(M_Z+K_ZM/P^2\right)-MW^*/P^2\right].
$$

该式已由完整势对带电场的二阶展开验证：实相位与虚相位两个展开分别核对 $\mathrm{Re}\,B$、$\mathrm{Im}\,B$。在 $m_G\gg H$ 的相关范围，即便交叉项的相位把轴向树级力消掉，仍有 $d\simeq m_G^2$、

$$
|B|\simeq m_G M.
$$

这里没有额外的接触项隔离、物质 Kähler 度量调节或实际标准模型嵌入。带电谱是指定最小作用量的结果。

## 4. 固定旧轨迹上的两个预算

第一项仅衡量旧实轴轨迹承受的额外树级力：

$$
R_{\rm tree}=\frac{C\sqrt{9\cos^2\theta+\sin^2\theta}}{2cU}.
$$

用父实验完整 1025 点轨迹，要求全区间 $R_{\rm tree}\le0.1$ 得

| 相位 | $m_G$ 上界 |
|---|---:|
| $0,\pi$ | $5.58605\times10^{-35}\ \mathrm{eV}$ |
| $\pi/2$ | $1.67582\times10^{-34}\ \mathrm{eV}$ |

第二行只是实轴偏离预算，不能用来否定上一节的重场谷底。

第二项是指定有限局部 Coleman–Weinberg 项的形式预算。令

$$
f_0(x)=x^2[\ln(x/\mu^2)-3/2],\quad
U_1=\frac{f_0(s+d+|B|)+f_0(s+d-|B|)-2f_0(s)}{32\pi^2},
$$

并固定 $\mu=m_i=100\ \mathrm{GeV}$。在 $H\ll m_G\ll M$ 时，记 $L=\ln(s/\mu^2)$，主项为

$$
U_1\simeq\frac{m_G^2s(6L-4)}{32\pi^2},\qquad
\boxed{U_{1,\chi}\simeq\frac{m_G^2s_{,\chi}(3L+1)}{16\pi^2}}.
$$

这里必须同时保留公共质量 $d$ 和双线性分裂 $B$。只给两标量加入公共软质量，会得到不同系数。固定旧轨迹上的主项估计为

$$
\max\frac{|U_{1,\chi}|}{|V_{\chi,\rm old}|}
\simeq 2.51517\times10^{32}\,(m_G/\mathrm{eV})^2.
$$

其 10% 上界约 $1.99396\times10^{-17}\ \mathrm{eV}$。这是指定有限项和固定重整化尺度的预算，不是与 UV 完成、反项或重整化条件无关的绝对自然性定理。尤其是低破缺尺度下，下面的 EFT 检查会限制此形式计算的物理解释。

## 5. Nilpotent 描述的能标限制

Goldstino 的导数相互作用给出数量级

$$
\Lambda_G\sim\sqrt{|F_X|}\sim\sqrt{|f|}
=(3m_G^2P^2+\rho_\Lambda)^{1/4}
$$

的强耦合标记；不是已知精确系数的硬截止。Komargodski–Seiberg 式 (3.6) 的四 Goldstino 算符直接显示其幂次计数；Dall'Agata–Zwirner 第 5 节讨论了 nilpotent 引力理论的相应尺度、背景依赖以及完整超多重态可改变高能范围。[From Linear SUSY to Constrained Superfields](https://arxiv.org/pdf/0907.2441)，[On sgoldstino-less supergravity models of inflation](https://arxiv.org/pdf/1411.2605)

本候选以小的晚期密度固定 $f$，故冻结参数给出

$$
\rho_\Lambda=2.51814\times10^{-11}\ \mathrm{eV}^4,\qquad
\Lambda_G(m_G=0)\sim2.24011\ \mathrm{meV}.
$$

这里的 meV 结论依赖 $|f|^2=3m_G^2P^2+\rho_\Lambda$ 的约束；仅说 $m_{3/2}\sim H_0$，而不固定真空能或其他部门，并不能推出它。

若把 $\Lambda_G\gtrsim100\ \mathrm{GeV}$ 暂作数量级一致性标记，需

$$
m_G\gtrsim2.37105\times10^{-6}\ \mathrm{eV}.
$$

这一门槛比实轴树级 10% 上界高约 28 个数量级；即使考虑 $\pi/2$ 的经典重场谷底，形式局部环修正在此门槛约为旧树级力的 $1.414\times10^{21}$ 倍。由此可见，当前最小候选中“保留原慢时钟”“直接解析 100 GeV 阈值”和“该选定有限环修正小”没有共同的已展示窗口。

这不是禁止高能粒子存在的结论。低能有效理论可以包含来自更高能粒子的 Wilson 系数；只是这些系数的匹配需要更高能完成，不能由低截止的 nilpotent 作用量自行计算。本候选也没有排除带有额外场、隔离机制、不同物质度量或不同反项条件的其他完成。

## 6. 可支持的结论

本轮的独立代数核查支持三个明确判断：隐藏部门的常数超势会通过交叉项改动时钟势；特殊相位可以形成重横向场谷底，但不消除所有反馈；带电部门仍接收 $m_G^2$ 和 $m_GM$ 量级的破缺，且重阈值的计算必须检查 EFT 能标。

因此，下一步如继续建设物理模型，应具体给出抑制传递的作用量或更高能完成，并重新匹配，而不能把相位选择、初始常数扣除或形式上的低软质量预算当作已经解决现实部门的问题。上述判断均限于本次明确指定的候选，不是对逆对数尺度设想的一般否定。
