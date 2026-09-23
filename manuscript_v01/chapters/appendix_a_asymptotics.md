# 附录A 有限格点纯包络模型的晚时行为

## A.1 假设与前向全局存在

本附录只考虑$\varepsilon_{\rm ar}=0$的式(3.1)—(3.4)。假设$N$为固定正整数，$L>0$、$a=L/N>0$，$I>0$、$g\ge0$、$v^2\ge0$，并要求

$$
R_i>R_*>0,\qquad u_i=\dot R(0)>0,\qquad
0<b<m_0^2\chi_i^2,\qquad
M_\infty^2=m_0^2-\frac{b}{\chi_i^2}>0.
\tag{A.1}
$$

初始场坐标和动量有限，总能量$H(0)$有限。以下常数可以依赖固定参数及初值，不主张在$u_i\to0$、$M_\infty^2\to0$或$N\to\infty$时一致有界。

在$R>R_*$上，有限维右端光滑，局部解唯一。纯包络反馈非负；在局部解的整个存在区间内有

$$
\ddot R=\frac{bQ_2}{IR\chi^3}\ge0,\qquad
\dot R\ge u_i,\qquad R\ge R_i+u_i\tau,
\qquad M_{\rm env}^2=M_\infty^2+\frac{b}{\chi^2}\ge M_\infty^2.
\tag{A.2}
$$

故轨迹始终远离对数奇点，全部Hamilton能量项非负。记$E_0=H_{\rm field}(0)$；式(3.5)表明场能单调下降，因而

$$
\begin{aligned}
0&\le H_{\rm field}(\tau)\le E_0,\\
Q_2(\tau)&\le Q_{\max}:=\frac{2E_0}{M_\infty^2},\\
u_i&\le\dot R(\tau)\le U:=\sqrt{u_i^2+\frac{2E_0}{I}}.
\end{aligned}
\tag{A.3}
$$

对于固定$a>0$，场能还分别控制每个$q_\ell$和$p_\ell$。在任意有限时间区间，$R\le R_i+U\tau$；状态因此位于不接触奇点的有界区域。常微分方程的延拓条件排除有限时间终止，证明前向全局存在。该论证并不涵盖反向时间，也没有建立连续PDE的统一正则性估计。

## A.2 时钟速度的正极限

因为$H_{\rm field}$单调有下界，存在$E_\infty\in[0,E_0]$。能量守恒与正速度分支给出

$$
\dot R(\tau)\longrightarrow
u_\infty=\sqrt{u_i^2+\frac{2(E_0-E_\infty)}{I}}
\ge u_i>0,
\qquad
\frac{R(\tau)}{\tau}\longrightarrow u_\infty.
\tag{A.4}
$$

第二个极限来自$R(\tau)=R_i+\int_0^\tau\dot R(s)\,\mathrm ds$及收敛函数的时间平均。它不要求各场模式趋向静止；场可在有限剩余能量下继续振荡。时钟获得的总能量至多为$E_0$。

## A.3 尾部误差与逆对数平方响应

由于$\dot R>0$且$R\to\infty$，可沿轨迹把速度$u=\dot R$视为$R$的函数。链式求导和式(A.3)给出

$$
\frac{\mathrm d(u^2)}{\mathrm dR}
=\frac{2bQ_2}{IR\chi^3},\qquad
0\le u_\infty^2-u^2(R)
\le\frac{2bQ_{\max}}{I}
\int_R^\infty\frac{\mathrm d\widetilde R}
{\widetilde R\ln^3(\widetilde R/R_*)}
=\frac{bQ_{\max}}{I\chi^2}.
\tag{A.5}
$$

因此$0\le u_\infty-u\le bQ_{\max}/(2Iu_i\chi^2)$。定义

$$
\tau_{\rm eff}=\frac{R_*}{u_\infty},\qquad
\mathcal L(\tau)=\ln(\tau/\tau_{\rm eff}).
\tag{A.6}
$$

只在充分晚时使用$\mathcal L$。由式(A.2)—(A.3)的上下线性界，$\chi=\mathcal L+O(1)$，故式(A.5)首先推出$u_\infty-u=O(\mathcal L^{-2})$。又有

$$
\int_{\tau_0}^{\tau}\mathcal L(s)^{-2}\,\mathrm ds
=O\!\left(\tau\mathcal L(\tau)^{-2}\right),\qquad
\frac{\mathrm d}{\mathrm d\tau}
\left[\frac{\tau}{\mathcal L(\tau)^2}\right]
=\mathcal L^{-2}-2\mathcal L^{-3},
\tag{A.7}
$$

其中$\tau_0$取到足够大，使对数为正。式(A.7)可由其右侧导数与被积函数的比趋于1得到。积分速度差，并将有限初值项包含在余项中，遂有

$$
R(\tau)=u_\infty\tau+O\!\left(\tau\mathcal L^{-2}\right),\qquad
\chi(\tau)=\mathcal L+
\ln\!\left[\frac{R(\tau)}{u_\infty\tau}\right]
=\mathcal L+O(\mathcal L^{-2}).
\tag{A.8}
$$

最后对$z^{-2}$使用中值定理：其在$z\sim\mathcal L$处的导数为$O(\mathcal L^{-3})$，而式(A.8)的自变量误差为$O(\mathcal L^{-2})$。因此

$$
\boxed{
M_{\rm env}^2(R(\tau))-M_\infty^2
=\frac{b}{\mathcal L(\tau)^2}+O(\mathcal L(\tau)^{-5})
}.
\tag{A.9}
$$

这证明命题1，并给出强于渐近等价的余项。式(A.8)只控制$R-u_\infty\tau=O(\tau/\ln^2\tau)$，并不保证这一累积位置差有界。固定时间平移对晚时对数的影响为$O(\tau^{-1})$，不改变主项，但不能由此用一次时间平移吸收全部反馈历史。$\tau_{\rm eff}$也只是由完整解决定的渐近尺度，不是全时间外部轨迹的自由拟合参数。

## A.4 指数选择与紧支撑残差的边界

上述论证说明已选响应能够在反馈下保留，而不证明平方指数唯一。若另行选择$s>0$及正下界，令$M_s^2=M_{s,\infty}^2+b\chi^{-s}$，同样的积分可得

$$
0\le u_\infty^2-u^2\le\frac{bQ_{\max}}{I\chi^s},\qquad
M_s^2(R(\tau))-M_{s,\infty}^2
=b\mathcal L^{-s}+O(\mathcal L^{-(2s+1)}).
\tag{A.10}
$$

这里$Q_{\max}$、$u_\infty$及$\mathcal L$必须按新模型重新定义。因而指数的物理选择仍需额外机制或外部检验。

对于紧支撑残差，势能中新增$\varepsilon_{\rm ar}wS$可能改变反馈符号，并可能破坏所需正能条件；不能沿用式(A.2)证明它必然穿过支撑。只有独立建立以下条件后，尾部论证才可恢复：解在有限$\tau_e$到达残差支撑上端之外，状态有限，$\dot R(\tau_e)>0$，且该区域具有同一个正的纯包络下界。此后$[wS]'=0$、$\ddot R\ge0$，轨迹不会返回支撑，以$\tau_e$为起点重复证明即可。有限起点平移不改变式(A.9)的晚时阶，但这些出支撑条件不能由残差幅度小或有限窗口数值稳定自动推出。
