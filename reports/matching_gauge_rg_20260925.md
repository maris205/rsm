# 带电超场阈值与选定 CW 部门的 RG 闭合

2026-09-25。独立接口核查，承接 [no-scale 量子接口](noscale_quantum_sources_20260925.md) 和 [KL 匹配](gauge_matching_theory_20260925.md)。本报告给出两项有界结果：**重带电对退耦后，补偿场项和物质度规线性辅助项可以抵消，但场依赖超势质量仍留下规范阈值；选定 CW 行列式的显式尺度依赖可以由局部反项精确抵消，但有限匹配系数不能因此确定。** 这不是完整的 $O(Q^2)$ 超引力一环有效作用量。

## 1. 约定与适用范围

$P=M_P^2$，$k=K_0/P$，$F^i=-e^{k/2}K^{i\bar j}D_{\bar j}\bar W$，并定义本场基底中的 $F_{\rm comp}=e^{k/2}\bar W/P+K_iF^i/(3P)$。只考虑电荷 $\pm1$ 的一对 $Q_\pm$；超势质量为 $M(Z)Q_+Q_-$，无 Kähler 双线性质量来源，无其他轻带电粒子。

本文取 $\beta_g=b g^3/(16\pi^2)$，所以 $b_H=2,b_L=0$。规范作用量与费米子相位约定为 $\frac14\int d^2\theta f W^\alpha W_\alpha+\mathrm{h.c.}$、$W_\alpha|=-i\lambda_\alpha$、$\mathcal L\supset-\frac12M_\lambda\lambda\lambda+\mathrm{h.c.}$；直接规范函数项因而是 $M_\lambda^{\rm direct}=g^2F^i\partial_i f/2$。标量势中的双线性定义为 $V\supset b_Q q_+q_-+\mathrm{h.c.}$。不同文献的规范超场和 gaugino 相位记号可能改变中间式的整体符号，以下相对符号由常质量退耦和全纯规范函数两条路径共同固定。

这是固定常背景、近超对称重阈值、导数展开首阶的匹配。要求 $|d|/s\ll1,|b_Q|/s\ll1$，$s=|\widehat M|^2$；高阶分裂、时空导数及完整引力环图没有由此算完。

## 2. 先把质量超场规范归一化

把一般物质度规写成

\[
\mathcal Z_r=e^{k/3}e^{\ell_r},\qquad
\widehat M=e^{k/2}\frac{M}{\sqrt{\mathcal Z_+\mathcal Z_-}}
=e^{k/6}M e^{-(\ell_++\ell_-)/2}.
\]

在选定背景附近，令局部手征归一化因子 $N_r$ 满足
$N_r|=e^{\ell_r|/2}$，$F^{N_r}/N_r=F^i\partial_i\ell_r$。这是去除动能中线性辅助分量的局部构造，并非把任意实度规写成全局手征函数。$Q_r^c=C N_r Q_r$ 后，超势中的阈值超场是

\[
\boldsymbol{\mathcal M}=\frac{\boldsymbol C\boldsymbol M}{\boldsymbol N_+\boldsymbol N_-},\qquad
\frac{F^{\mathcal M}}{\mathcal M}
=F_{\rm comp}+A_M-A_\Sigma\equiv S_M,
\]
\[
A_M=F^i\partial_i\ln M,\qquad
A_\Sigma=F^i\partial_i(\ell_++\ell_-),\qquad
\frac{b_Q}{\widehat M}=-S_M.
\]

共同软质量 $d$ 来自额外的 $\theta^2\bar\theta^2$ 信息，不能由 $S_M$ 单独恢复。上述手征重标度及对应异常项的原始框架见 [D’Eramo–Thaler–Thomas，§2、§3，尤其式 (23)–(32)](https://arxiv.org/html/1202.1280v3)。在本项目中，最后一行也可直接由既有规范归一化软项公式导出，因而是对特定质量来源的限制，不是所有双线性模型的恒等式。

## 3. $b_L=0$ 不等于所有低能阈值为零

将完整的 super-Weyl、Kähler 和 Konishi 组合写在阈值以上，在本谱和上述记号下为

\[
M_{\lambda,H}=M_\lambda^{\rm UV}
+\frac{2g^2}{16\pi^2}(F_{\rm comp}-A_\Sigma).
\]

带电对的质量超场匹配给

\[
\Delta M_{\lambda,Q}=-\frac{2g^2}{16\pi^2}S_M
=\frac{2g^2}{16\pi^2}\frac{b_Q}{\widehat M},
\]
\[
\boxed{M_{\lambda,L}=M_\lambda^{\rm UV}
-\frac{2g^2}{16\pi^2}F^i\partial_i\ln M.}
\]

这里各侧耦合以匹配点的一环精度取同一 $g$，或更严格地比较 $M_\lambda/g^2$。$M_\lambda^{\rm UV}$ 表示未计算的独立规范函数边界；若声称它为零，需要明确采用固定全纯 UV 规范函数这一条件。普通近超对称阈值的辅助分量必须一起匹配，原始论证见 [Giudice–Luty–Murayama–Rattazzi，§3.1，式 (10)–(15)](https://arxiv.org/html/hep-ph/9810442)。本报告已换用第一节的 $b$ 和费米子相位约定，未逐字搬用其超场分量符号。

有两个有辨识力的极限：$M$ 与中性场无关时，$A_M=0$，阈值会取消重对带来的补偿量和线性物质度规项；$M(Z)$ 有辅助分量时，$A_M\ne0$，即使 $b_L=0$ 仍有非零局部 gaugino 质量系数。后者是匹配遗留项，不能用低能异常介导的空谱公式抹去。$d$ 以及更高阶分裂会给出这里未保留的非全纯修正。

第二条独立推导是

\[
f_L=f_H-\frac{2}{8\pi^2}\ln\frac{M}{M_*},
\qquad
\frac{g^2}{2}F^i\partial_i f_L
=\frac{g^2}{2}F^i\partial_i f_H
-\frac{2g^2}{16\pi^2}A_M.
\]

$M_*$ 为固定参考尺度。质量对数的符号、系数及 Kähler 项抵消来自 [Kaplunovsky–Louis，式 (3.27)–(3.29)](https://arxiv.org/pdf/hep-th/9402005v2)，本次重新检查了缓存 PDF 的相应页。真正进行 Kähler 或物质场全纯重定义时，$f_H$ 的异常变换必须一并携带；$f_H=\mathrm{const.}$ 只是一个指定表示中的边界。脚本分别核查了 Kähler、Konishi 的这一不变性。

在四次几何的 $g\simeq1,y=0$ 初始中心，

\[
A_M\simeq\frac{\sqrt{2U}}{F_\chi}
\frac{\xi}{\chi(\chi^2+\xi)}>0.
\]

使用已有输入，条件式给出 $M_{\lambda,L}\simeq-6.08\times10^{-38}\,\mathrm{eV}$，仅用于确认“非零残留”的量级。这里的负号是第一节相位约定下的参数符号，不是负物理质量。它没有包含未知 UV 规范算符，也不是现实电磁理论、真实 gaugino 或观测预言。

## 4. 实际质量层次决定匹配顺序

$m_G=1\,\mathrm{eV},a=1$ 样本有 $m_T\simeq\sqrt{48}\,m_G=6.928\,\mathrm{eV}$，而 $|\widehat M_Q|\simeq10^{11}\,\mathrm{eV}$，两者之比为 $1.4434\times10^{10}$。所以适用的 Wilson 匹配顺序是：在 $T$ 仍活跃的理论中处理约 $100\,\mathrm{GeV}$ 的带电对；在中间能区保留 $T$ 及由重对生成的局部算符；随后在约 $6.9\,\mathrm{eV}$ 处理重模量标量。完整 $T$ 多重态的Goldstino（超对称破缺费米子）与引力的混合另须处理，不能仅由两条标量质量宣称全多重态退耦已经完成。

也可以直接算统一常背景 1PI 行列式，但那是另一种计算组织，必须保留相应轻重混合贡献。先积分掉轻 $T$ 后得到的低能局部展开，不能直接拿去描述 $100\,\mathrm{GeV}$ 带电粒子的在壳过程。若改变 $m_G$ 使层次交换，匹配顺序也需要改变。

## 5. 带电 CW 的尺度依赖可以精确闭合

令 $h=|b_Q|$，$m_\pm^2=s+d\pm h>0$。对两个复标量和一个 Dirac 费米子，选定有限减法的行列式为

\[
V_{1,Q}=\frac{F(s+d+h)+F(s+d-h)-2F(s)}{32\pi^2},
\quad F(x)=x^2[\ln(x/\mu_Q^2)-3/2].
\]

直接求导，无需小分裂展开，即得

\[
\boxed{\left.\frac{\partial V_{1,Q}}{\partial\ln\mu_Q}\right|_{s,d,h}
=-J_Q,\qquad J_Q=\frac{2sd+d^2+h^2}{8\pi^2}.}
\]

这也能由独立的软矩阵公式核查：取 $m_{\rm soft}^2=dI_2$、超势质量矩阵的两个非对角元为 $M$、双线性矩阵对应元为 $b_Q$，一般真空能 beta 函数给 $[2d^2+4sd+2|b_Q|^2]/(16\pi^2)$，恰好等于 $J_Q$。一环势、减法方案与真空能 beta 函数的原始研究来源为 [Martin，hep-ph/0111209v2，式 (3.2)–(3.13)、(7.1)–(7.11)](https://arxiv.org/html/hep-ph/0111209v2)。该文处理一般可重整化理论；这里仅用这些行列式和局部 RG 恒等式，不把其结果冒充整个非线性超引力模型的 beta 函数。

引入投影到本部门的局部势系数 $C_Q(\chi,y;\mu)$，在固定本阶谱的意义下

\[
C_Q(\chi,y;\mu)=C_Q(\chi,y;\mu_*)+J_Q(\chi,y)\ln\frac\mu{\mu_*}
\]

即可使 $V_{1,Q}+C_Q$ 及它对 $\chi,y$ 的导数都不依赖显式 $\mu$。谱参数的一环运行插入 $V_1$ 属下一环；完整 RG 还含参数运行和背景场异常维数，见来源式 (7.1)–(7.3)。这里的 $C_Q$ 是局部作用量在背景上的投影，不能误称为只需调整一个场无关宇宙学常数。

若只保留两个实模量标量，则同样有

\[
J_T=\frac{m_t^4+m_v^4}{32\pi^2},\qquad
\partial_{\ln\mu_T}V_{1,T}=-J_T.
\]

这仅闭合选定标量行列式的显式对数。既没有据此计算完整超引力的局部反项，也没有确定任何反项的有限值。

改变减法常数 $3/2\to3/2+\Delta c$ 会令
$\Delta V_{1,Q}=-(\Delta c/2)J_Q$；同一理论的局部系数必须作相反平移。因此只改变 CW 中的 $\mu$ 或减法常数、却固定或遗漏匹配系数，是改变减法边界的诊断，不是物理量的合法尺度扫描。之前那些随尺度变化的力预算根，不能解释为新的粒子质量限值。

## 6. 对主实验选定系数基底的独立化简

主实验取 $d=-cX$、$h^2=\beta^2sX$、$p=sX$。保留 $sd,h^2$ 而明确略去 $d^2$ 及更高阶分裂，有

\[
V_{1,Q}^{\rm sel}=\frac{p}{8\pi^2}[A L+c],\quad
A=-c+\frac{\beta^2}{2},\quad L=\ln\frac{s}{\mu^2}.
\]

用局部项 $p\nu/(8\pi^2)$ 补齐后，

\[
\nu(\mu)=\nu_0+2A\ln\frac{\mu}{M_i}
\]

恰好取消 $L$ 的尺度变化。这与第五节的 $J_Q$ 完全一致；被略去的 $d^2$ 对应另一个 $X^2$ 局部结构，不能在不说明阶数的情况下丢掉。

$\partial_cV\equiv K_c$、$\partial_{\beta^2}V\equiv K_{\beta^2}$、$\partial_\nu V\equiv K_\nu$ 满足

\[
K_c+2K_{\beta^2}=K_\nu.
\]

因此这里三种系数在该阶只有两种独立势函数；$c\to c+\delta,\beta^2\to\beta^2+2\delta,\nu_0\to\nu_0-\delta$ 不改变这个选定势。谱和其他算符一般仍会改变，所以不能把这称为整个理论的物理等价。固定 $\nu_0=0$ 后的预算范围是明确边界下的条件结果；改变 $\nu_0$ 是改变有限匹配输入，与只换 $\mu$ 的操作不同。

## 7. 没有被这些恒等式完成的工作

本报告完成规范阈值的首个辅助分量和 $Q=0$ 的选定真空行列式对数核对。$O(Q^2)$ 物质度规、Kähler 双线性、混合重轻图、辅助场一致消除、完整引力/ghost 以及有限 UV 规范函数仍需匹配。模量一环生成的 $\delta d$ 或 $\eta_{\rm rad}$ 对带电一环的线性插入，以及一环 $\delta b_Q$ 与树级 $b_{Q,0}$ 的干涉，可以属于二环子集；若仅取一环 $\delta b_Q$ 的绝对值平方，则通常已到原理论三环阶。层次增强不会自动使同阶其余项消失。

一个参考真空能条件只能固定一个组合，不能同时确定所有局部场依赖系数及两个时钟方向的导数。RG 闭合解决的是任意减法尺度不应成为物理结论的问题，并不证明有限 Wilson 系数小、保护具有技术自然性或量子稳定已经完成。

## 8. 验证与材料记录

[独立脚本](matching_gauge_rg_20260925_checks.py) 不导入主实验计算函数。`python reports/matching_gauge_rg_20260925_checks.py` 得到 [34/34 项通过](matching_gauge_rg_20260925_checks.json)：精确 CW 尺度导数、独立矩阵收缩、谱导数闭合、有限减法平移、三阶分裂展开、阈值两路径、常质量及遗漏阈值负对照、场基底变换、主选定基底退化和九个 100 位精度直接行列式尺度比较。直接数值比较最大相对残差约 $6.81\times10^{-26}$；检查项互有关联，不是 34 条独立物理证据。

[来源清单](matching_gauge_rg_20260925_sources.json) 保存 URL、实际读取范围、缓存哈希和检索异常；文献全文缓存在忽略的 `build/` 下。初始文本搜索有一次正则转义错误，改用多个字面模式后成功；KL HTML 返回内部错误后使用先前已核验 PDF，这些不会影响推导结果。无新观测，无统计显著性检验，无外部审稿；属于 AI 辅助的项目内部独立推导与复核。

## 9. 对 `eft_matching_v01` 主实现的有界对抗性复核

本次另逐段读取了 `experiments/eft_matching_v01/protocol.md`、`code/run_matching.py` 与已生成的核、预算和 RG 摘要。主代码把有限边界 $\nu_0$ 的改变与仅改变重整化尺度分开，正确标明单系数预算不构成联合界；$c,\beta^2,\nu$ 的响应核、初始 $L=0$ 的稳定 $c$ 核及二维 $X$ 一阶导数与本报告独立推导一致。尺度检查比较相对于 $2U$ 的绝对误差，没有把极小零附近的分量相对误差伪装成精确。

针对容易混淆的基线 $B$ 干涉项，[另一份独立符号脚本](matching_gauge_rg_20260925_adversarial.py) 从 $W(\chi+iy)$、$M(\chi+iy)$ 的局部展开及 $S=3W-\bar k_Zw_Z$ 重新求导，不导入生产函数；[10/10 恒等式通过](matching_gauge_rg_20260925_adversarial.json)。特别是基线 $B_0$ 的 $\eta$ 部分含 $M-\bar k_ZM_Z$，所以其横向导数带 $-rF^T$；额外项 $\Delta B=\beta M F^T$ 则带 $+rF^T$。主代码中这两个不同符号是正确的，不能为了外观一致而改成同号。它们共享 $\partial_yF^T=-i m_G r_w/3$，与完整辅助组合给出的三分之一因子相容。

**有界结论：未发现上述选定部门的公式或环阶声明存在新的阻断问题。** 该结论仅覆盖局部中心的一阶 jet、指定带电行列式近似和明确有限边界；不等于完成全部 $O(Q^2)$ 匹配，也不替代主实验独立高精度核验的最终通过状态。审核时获知独立数值代理正在保留原始失败记录后提高横向差分精度；这一生产归档门槛应以其最终记录为准。本节自身的 10 项检查采用精确符号约简，没有修改数值容差。
