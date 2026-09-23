# 附录C 不同驱动候选与条件观测接口

## C.1 对数坐标不等于恒定动能系数的标量

本稿的时钟为自由内部坐标$R$，其动能为$I\dot R^2/2$。令$R=R_*e^\chi$，正则一形式要求$P_R\,dR=P_\chi\,d\chi$，故

$$
P_\chi=RP_R,\qquad
E_R=H_{\rm clock}=\frac{P_\chi^2}{2IR_*^2e^{2\chi}},\qquad
\mathscr L_R=\frac12IR_*^2e^{2\chi}\dot\chi^2.
\tag{C.1}
$$

因此，对$\chi$的变量替换保留了原模型，却给出随坐标变化的动能系数。它不能直接换成常数$F_\chi^2$，否则就是另一个候选。纯包络的时钟方程相应变为

$$
IR^2(\ddot\chi+\dot\chi^2)=\frac{bQ_2}{\chi^3}.
\tag{C.2}
$$

作为对比，在均匀膨胀背景上另行选择具有恒定正动能系数的标量，势取$V_s e^{-2\chi}$。忽略物质源后，方程为

$$
F_\chi^2\left[\frac{d^2\chi}{dt^2}
+3H_{\rm cos}(t)\frac{d\chi}{dt}\right]
-2V_s e^{-2\chi}=0.
\tag{C.3}
$$

$H_{\rm cos}$在此表示宇宙膨胀率，区别于本稿总Hamilton量$H$。在给定$H_{\rm cos}=h_{\rm cos}/t$的背景上，直接代入可知$\chi=\ln(t/t_*)$要求

$$
V_s t_*^2=\frac{F_\chi^2}{2}(3h_{\rm cos}-1).
\tag{C.4}
$$

当$V_s>0$时，需要$h_{\rm cos}>1/3$。若联立引力与物质方程，还须满足背景能量份额与流体标度条件；指数势的相应标度解已有标准研究。[7](../references.md#r7) 此处只用式(C.3)—(C.4)说明候选之间的区别，不把这一背景解当作式(C.1)的另一写法，也不在当前论文中合并两者的数值约束。

## C.2 当前漂移与过去变化的消参关系

第6章的响应假设以$t_0$为参考时刻，$\alpha_0=\alpha(t_0)$。令$\mathcal L_\alpha(t)=\ln(t/t_*)$，在固定历史、$s>0$及$1+\delta_\alpha>0$的范围内，由式(6.1)求导得到

$$
\frac{d\ln\alpha}{dt}
=-\frac{s\Gamma}{t\mathcal L_\alpha(t)^{s+1}[1+\delta_\alpha(t)]},
\qquad
D_0=-\frac{s\Gamma}{t_0\mathcal L_{\alpha0}^{s+1}}.
\tag{C.5}
$$

用第二式消去模型幅度$\Gamma$，即得到式(6.2)。这是模型参数之间的代数关系，不是用一个实测漂移均值去除另一个观测量。$D_0=0$仍是允许的模型值，应与所有数据一起拟合。

若改用某条已标定的实际背景$\chi_{\rm phys}(t)$，令$\chi_0=\chi_{\rm phys}(t_0)$，并要求整个使用区间内$\chi_{\rm phys}(t)>0$及$1+\delta_\alpha(t)>0$。另外假定$\delta_\alpha=\Gamma[\chi_{\rm phys}^{-s}-\chi_0^{-s}]$后，在$\chi'_0=(d\chi_{\rm phys}/dt)_{t_0}\ne0$时有

$$
D_0=-\frac{s\Gamma\chi'_0}{\chi_0^{s+1}},\qquad
\delta_\alpha(t)=-D_0\frac{\chi_0^{s+1}}{s\chi'_0}
\left[\chi_{\rm phys}(t)^{-s}-\chi_0^{-s}\right].
\tag{C.6}
$$

当当前背景恰好停住时，不能用$D_0$唯一参数化过去变化，应保留$\Gamma$；若反馈又使背景依赖$\Gamma$，式(C.6)虽然沿每条解成立，传递关系却不再是幅度无关的固定函数。本稿仍需物理机制才能把内部时钟和这种背景联系起来。

对实际跃迁，定义$K_a=\partial\ln\nu_a/\partial\ln\alpha$，其他相关参数固定，则链式求导给出$d\ln(\nu_a/\nu_b)/dt=(K_a-K_b)d\ln\alpha/dt$。数值应用必须由独立原子计算或测量确定$K_a,K_b$，并核验慢变近似与环境条件，不能把敏感度作为自由幅度补偿任意时间曲线。
