# 独立复算的方法和边界

2026-09-25。内部数值复核，由独立代理实现；不等于外部同行评审，也不产生新的观测证据。

## 材料与问题

输入是本目录冻结的 `protocol.md`、`results/inputs.json`、主求解器保存的背景和线性模 CSV。九组初态和四个波数均照用预定值，没有因结果而改变参数。旧实验数据和正式论文均未修改。生产求解器与本复核程序共享的是物理问题、参数、比较网格和结果格式，不共享 ODE 实现。

## 不同积分变量的背景核查

主求解器使用 $N=\ln a$。本程序以 $\tau=H_{\rm ref}t$ 为自变量，独立积分 $(N,u,y,w_\chi,w_y)$，其中 $u=\chi-\chi_i$、$w_\chi=\dot\chi/H_{\rm ref}$ 和 $w_y=\dot y/H_{\rm ref}$。方程为

\[
\frac{dN}{d\tau}=E,\quad
\frac{du}{d\tau}=w_\chi,\quad
\frac{dy}{d\tau}=w_y,\qquad
\frac{dw_A}{d\tau}=-3Ew_A-V_{,A},
\]

\[
E^2=\Omega_r e^{-4N}+\Omega_m e^{-3N}+\Omega_\Lambda+
\frac{\epsilon}{3}\left[\frac{w_\chi^2+w_y^2}{2}+V\right].
\]

$V$ 的单位为 $F_\chi^2H_{\rm ref}^2$。所有模型匹配相同的物理初速；不把相同 $d\chi/dN$ 当成相同物理初速。到达 $N=0$ 时停止，再用单调的 $N(\tau)$ 反解主程序的网格。另积分损失量 $L'=3E(w_\chi^2+w_y^2)$，核查 $\rho_\chi+L$ 守恒。这是连续性方程的数值恒等式，不意味着标量部门在膨胀背景中保持独立不变的能量密度。

## 横向线性模

使用本程序自己的轴向背景，独立积分

\[
\frac{d^2\delta y}{d\tau^2}+3E\frac{d\delta y}{d\tau}+
\left[\kappa^2 e^{-2N}+V_{,yy}\right]\delta y=0,
\quad\kappa=k_{\rm com}/H_{\rm ref}.
\]

主程序约定的两个基解初态 $(\delta y,\delta y_N)=(1,0),(0,1)$ 在这里分别转换成 $(\delta y,\delta y_\tau)=(1,0),(0,E_i)$。比较时除回本程序的 $E$ 得到 $\delta y_N$。因此误把不同时间变量的初速视为相同这一实现错误能被此复核发现。

## 不依赖简化实势公式的核查

复场取 $z=Z/F_\chi=(\chi+iy)/\sqrt2$，$k=K/M_{\rm Pl}^2$、$w=W/(F_\chi\sqrt{U_i})$。以 70 位精度直接计算

\[
V/U_i=e^k\left[\frac{\epsilon}{k_{z\bar z}}
|w_z+k_zw|^2-3\epsilon |w|^2\right].
\]

计算 $k_z$、$k_{z\bar z}$ 时把 $z$ 与 $\bar z$ 当成独立变量，然后代入共轭关系。对这一个完整复表达式进行高精度实变量求导，获得势、梯度和 Hessian。全局对照单独采用 $V/U_i=|w_z|^2$。预定采样是 $u=0,0.75,2$ 与 $y=0,0.001,0.1$ 的笛卡尔积，每个模型九点。分别对照本程序的符号求导和主程序的解析函数。主程序的势函数只在这一步用于比较，没有进入独立 ODE。

## 数值精度和解释范围

proper-time 积分使用 DOP853，`rtol=2e-12`、`atol=2e-14`、最大步长 $5\times10^{-4}$。背景和线性模每一分量的验收误差是

\[
\frac{\max|x_{\rm independent}-x_{\rm primary}|}
{\max(1,\max|x_{\rm primary}|)}<10^{-6}.
\]

所有比较量均已使用前述无量纲单位；这个定义允许跨零点，不把相对误差的除零发散误判成不一致。势核查容差为 $2\times10^{-11}$；独立连续性积分核查容差为 $10^{-8}$。程序记录每个检验、通过数、输入和代码的 SHA-256，以及完整独立轨迹。结果见 `results/independent_checks.json`；执行结果和解释另记于 `independent_results_cn.md`。

这些检验只能确认给定树级玻色部门的代数与数值一致性，不验证完整超引力宇宙、可见物质耦合、量子保护的所有高阶贡献、一般初态稳定域，或实验允许性。正横向曲率也不自动意味着重场、快速衰减或所有扰动部门均稳定。
