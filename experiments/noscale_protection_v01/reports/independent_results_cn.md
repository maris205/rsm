# 四次几何保护：独立高精度复核

本复核使用冻结协议和旧输入，未读取或导入主程序 run_protection.py。先直接从完整 Kähler 度规与超势建立高精度势、模量谷和局部谱，再读取主输出比较。独立结果支持主扫描所采用的局部、绝热近似及其条件预算；它没有把两个模量实标量的有限圈图扩展成完整超引力量子有效作用。

## Material Passport：材料与复算入口

- 输入：[协议](../protocol.md)、sugra_shift_v01 的 inputs.json、旧轨迹的 1025 个 shift_axis 坐标。没有重跑旧 ODE。
- 实现：[independent_check.py](../code/independent_check.py)。
- 精确点：[81 点 tree_eta，220 位](../results/independent_tree_eta_220.csv)；[27 点 scalar_retuned_eta，220 位](../results/independent_scalar_retuned_eta_220.csv)。相同点另有 160 位结果。
- [原始 CW 超迹结果](../results/independent_direct_CW.csv)、[精确点与主结果比较](../results/independent_exact_main_comparison.csv)、[最终摘要](../results/independent_summary.json)。
- 日期：2026-09-25；内部 AI 计算审计，非外部同行评审，无新增观测、统计显著性或自由拟合。

单次命令均以 timeout 180 执行。精确点可通过参数 --mode tree_eta --dps 220 或 --mode scalar_retuned_eta --dps 220 重建；--audit 执行精度、原始 CW、主摘要和预算根的复核。较低精度的三几何批次用 --a 0.1、1、10 分开保存。

## 完整 K/W 与容易漏掉的横向导数

令 $P=M_p^2$，$f=f(\delta,v)$，$p=f_T$，$q=f_{T\bar T}$，

$$
\Omega=f-\frac{F_\chi^2y^2}{3P},\qquad
k_Z=-i\sqrt2 F_\chi y,\qquad
\widetilde D=|p|^2-q\left(\Omega+\frac{|k_Z|^2}{3P}\right).
$$

直接求逆完整 $(T,Z)$ Kähler 矩阵得到

$$
V=\frac{1}{\Omega^2}\left[
|w_Z|^2+\frac{q}{3P\widetilde D}|S|^2
\right],\qquad
S=3W-\bar k_Z w_Z.
$$

此式同时用独立的矩阵逆、$D_IW$ 收缩与 $-3|W|^2/P$ 原表达复算，而非只在实轴检验。每个精确点还检验 $y=0.07$ 的轴外势。沿指定指数超势，

$$
S(y)=3W_c+w_{\rm real}(3+2iy)e^{-iy},\qquad S_y(0)=-iw_{\rm real}.
$$

若 $r_w=w_{\rm real}/(m_G P)<0$、$W_c=m_GP e^{+i\theta}$，则

$$
\left.\frac{\partial}{\partial y}\frac{|S|^2}{9P^2}\right|_0
=-\frac23m_G^2 r_w\sin\theta.
$$

它是同一约定下 $\partial_y|W|^2/P^2$ 的三分之一。若改用正数 $-r_w$，公式符号相应改变；主输出的实际符号与完整计算一致。不能把轴上的 $m_T^2\simeq48a|W|^2/P^2$ 直接当作轴外质量函数来求 $y$ 导数。

完整带电二次展开给出

$$
B=\frac1\Omega\left[
\overline{w_Z}M_Z+
\frac{q\bar S}{3P\widetilde D}
(M-\bar k_ZM_Z)
\right].
$$

独立计算直接差分整个轴外表达，包含 $M_y=iM_\chi$、$(M_Z)_y=i\sqrt2M_{\chi\chi}/F_\chi$、$S_y$ 及谷底位移。$B_y$ 不能由 $W$ 相位一项代替。

## 一次真空调节与模量谷

每个 $(a,m_G,\eta\text{方案})$ 先在 clock-off 完整势上同时求真空高度与驻点，之后不随 $\chi$ 再调 $\eta$：

$$
V_{\rm off}(\delta_0,\eta)=V_{\rm target},\qquad
\partial_\delta V_{\rm off}=0.
$$

tree_eta 取 $V_{\rm target}=\rho_\Lambda$；retuned 取协议指定的 $\rho_\Lambda-V_{1,T}^{\rm fixed}(w=0)$。这里没有加入完整圈图的 $Q^2$ 匹配或宣称完成所有反项重整化。

开启 clock 后重新解 $\partial_\delta V=0$，用 $\delta$ 的预测量级缩放变量，避免双精度的 1+delta 抹掉位移。隐函数导数是

$$
\delta_\chi=-V_{\delta\chi}/V_{\delta\delta},\qquad
\delta_y=-V_{\delta y}/V_{\delta\delta}.
$$

模量曲率、带电谱和圈图梯度均包含这项谷底导数。两个实模量的检查量是规范化后的 **T 方向局部势曲率**；其正性不等于全宇宙多场扰动谱、全局吸引域或稳定宇宙演化的证明。

81 个 tree_eta 点与 27 个 retuned 点均在 160 位、220 位重复。所有点满足正定度规、正模量局部曲率和正带电质量平方。一阶位移近似的最大相对差约 $3.08\times10^{-41}$；一次 $\eta$ 近似的最大相对差约 $5.64\times10^{-41}$。各点及最坏残差全部保留。

正交相位的极小 $\chi$ 模量质量导数含 $U/P$ 项，与只差分 $48a|W|^2/P^2$ 的次阶结果可以有较大的相对差；主导的 $y$ 梯度则由完整 K/W 重现。因此比较使用协议规定的非零二维向量／同组尺度，而不把趋零的次阶分量当作主预测。差异保留在逐点 CSV。

## 原始 CW、全扫描与预算

对 $a=b=1$ 的 27 个 tree_eta 点，以及额外 27 个 retuned 点，直接计算

$$
V_{1,Q}=\frac{f(s+d+h)+f(s+d-h)-2f(s)}{32\pi^2},
\quad f(x)=x^2[\ln(x/\mu_Q^2)-3/2],
$$

并用三个质量平方的完整导数求 $\chi,y$ 梯度。在 180 位与 240 位重复，消除重简并谱相减误差。稳定展开保留共同分裂至 $d^3$ 和双线性分裂至 $|B|^4$，与原始超迹比较目标为相对 $10^{-40}$。

全部 873 行基准扫描核为 5238 行方案扫描的对应子集。独立 NumPy 实现重建 **5238×1025=5,368,950** 个采样位置的全部摘要，并另求 54 个 10% 分量和预算根。这一步重建遵循协议保留的绝热阶次；108 个完整势精确点另外检验近似适用性，而非将同一近似的两次实现误称为精确理论验证。

最终汇总审计 **527/527** 通过；220 位精确点本身为 **540/540**，160 位重复也全部通过。最大误差如下，完整最坏位置、预算根与输入哈希见[独立摘要](../results/independent_summary.json)：

- 精确点 160→220 位相对变化：$4.82\times10^{-83}$。
- 180→240 位原始 CW 变化：$8.03\times10^{-94}$。
- 原始 CW 与稳定展开的最大相对差：$1.22\times10^{-105}$。
- 精确点与主输出的同组二维归一误差：$1.33\times10^{-14}$。
- 全扫描各摘要字段的最大相对差：$5.64\times10^{-15}$。
- 54 个独立预算根与主结果的最大相对差：$2.82\times10^{-12}$。

复核支持的是给定作用、参数、有限圈图分量和固定匹配方案下的条件预算。它不支持“完整量子保护已经证明”“现实精细结构常数变化已预言”或“黎曼零点已经导出真实相互作用”等更强结论。

## 实际发现与问题日志

1. 初版把已舍入的 $\xi,\chi_i$ 转成高精度十进制后，沿用半数质量系数，使初始质量偏离固定 $100\,\mathrm{GeV}$ 约 $7\times10^{-18}$。在 $L=0$ 处这能制造假导数。最终使用 $M^2=m_i^2(1+\xi/\chi^2)/(1+\xi/\chi_i^2)$，落实初始质量并重跑全部精确点。原尝试日志保留。
2. 初版仅线性分裂展开，在三个额外 retuned、$m_G=10^7\,\mathrm{eV}$、初始匹配点与直接核相差最高 $4.83\times10^{-9}$。这是 $d^2$ 等项在最低阶项接近零时可分辨的结果。保留 [523/526 原审计](../results/independent_initial_leading_expansion_audit.json)，随后补全稳定展开；没有放宽 $10^{-40}$ 门槛。
3. 添加元数据后，一次审计启动出现括号错误，在科学计算开始前纠正；[日志](../results/independent_attempt_audit_syntax.log)保留。报告写入还曾遇到一次工具侧 JavaScript 字符串分隔错误，重写后成功，没有改变科学结果。
4. 归档前仅删除函数声明的一处行尾空格，并重跑六个精确批次及顺序汇总审计以更新执行哈希。八个科学 CSV 逐字节相同、所有科学统计不变，七个最终 JSON 均匹配现有源码哈希；527/527 保持通过。[归档核验记录](../results/independent_archive_format_verification.json)及带 _archive_final 后缀的新日志保留，旧日志没有覆盖。

这些修正涉及数值归一、展开及记录机制，没有重新选择物理数据、隐藏不利结果或按结果调节 $\eta(\chi)$。
