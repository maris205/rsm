# 局部非线性正则格点到连续场：独立推导与限制

2026-09-22。本报告是 `micro_macro_v01` 的数学说明，方程由下述 Hamilton 量直接推导。它不报告实验运行结果，也不把形式连续极限当作完整收敛定理。本轮研究的是一个固定系数的已知型非线性场基准，用于检验后续微观—宏观研究的方法；尚未将黎曼零点、非自治驱动或真实粒子嵌入模型。

研究设定：周期区间 $x\in[0,L)$，$L=2\pi$，$v=\omega_0=g=1$，初值 $\phi(x,0)=A\cos x$、$\pi(x,0)=0$，$A=0.5$，终点 $T=20$。这里 $g$ 是四次势系数，与其他项目中同名的结构参数没有预设关系。

## 1. 格距与正则动量的归一化

令 $a=L/N$、$x_j=ja$，周期下标 $q_{j+N}=q_j$。取 $q_j=\phi(x_j)$，而正则动量为

$$
p_j=a\pi_j,\qquad \pi_j\simeq\partial_t\phi(x_j),
\qquad \{q_j,p_k\}=\delta_{jk}.
\tag{1}
$$

这是离散 Lagrange 量

$$
\mathcal L_a=a\sum_j\left[
\frac{\dot q_j^2}{2}
-\frac{v^2}{2a^2}(q_{j+1}-q_j)^2
-\frac{\omega_0^2q_j^2}{2}-\frac{gq_j^4}{4}
\right]
\tag{2}
$$

通过 $p_j=\partial\mathcal L_a/\partial\dot q_j$ 给出的结果。对应的 Hamilton 量为

$$
H_a(q,p)=\sum_j\left[
\frac{p_j^2}{2a}
 +\frac{v^2}{2a}(q_{j+1}-q_j)^2
 +a\left(\frac{\omega_0^2q_j^2}{2}+\frac{gq_j^4}{4}\right)
\right].
\tag{3}
$$

因此

$$
\dot q_j=\frac{p_j}{a},\qquad
\dot p_j=\frac{v^2}{a}(q_{j+1}-2q_j+q_{j-1})
-a(\omega_0^2q_j+gq_j^3),
\tag{4}
$$

$$
\ddot q_j=v^2\frac{q_{j+1}-2q_j+q_{j-1}}{a^2}
-\omega_0^2q_j-gq_j^3.
\tag{5}
$$

若程序储存 $(q_j,\pi_j)$，可以直接使用 $\dot q_j=\pi_j$、$\dot\pi_j=$ 式 (5) 右侧；但此时 $\{q_j,\pi_k\}=\delta_{jk}/a$。连续与离散辛形式的对应是

$$
\sum_jdq_j\wedge dp_j
=a\sum_jdq_j\wedge d\pi_j
\longrightarrow\int_0^L d\phi(x)\wedge d\pi(x)\,dx.
\tag{6}
$$

对一个固定均匀网格，辛形式前的整体常数 $a$ 不影响“是否保持辛形式”的判断；但它影响 Hamilton 量、括号及不同网格间物理量的归一化，不能任意丢弃。

## 2. 连续极限与首个空间修正

对光滑场进行 Taylor 展开得到

$$
\frac{\phi(x+a)-2\phi(x)+\phi(x-a)}{a^2}
=\phi_{xx}+\frac{a^2}{12}\phi_{xxxx}
+\frac{a^4}{360}\phi_{xxxxxx}+O(a^6).
\tag{7}
$$

固定 $L$ 并令 $a\to0$，主阶连续方程是

$$
\boxed{\phi_{tt}=v^2\phi_{xx}-\omega_0^2\phi-g\phi^3},
\tag{8}
$$

其 Hamilton 泛函为

$$
H_0[\phi,\pi]=\int_0^L\left[
\frac{\pi^2}{2}+\frac{v^2\phi_x^2}{2}
+\frac{\omega_0^2\phi^2}{2}+\frac{g\phi^4}{4}
\right]dx .
\tag{9}
$$

保留首个格距修正的方程为

$$
\boxed{\phi_{tt}=v^2\phi_{xx}-\omega_0^2\phi-g\phi^3
+\frac{a^2v^2}{12}\phi_{xxxx}},
\tag{10}
$$

注意四阶导数前为**正号**。相应的形式修正能量为

$$
H_{(2)}=H_0-\frac{a^2v^2}{24}\int_0^L\phi_{xx}^2\,dx.
\tag{11}
$$

式 (7) 是一致性展开；把残差阶提升为解误差阶还需要有限时间内的光滑性、稳定性和高频控制。本基准用加密实验验证所选初值下的误差阶，不声称仅凭 Taylor 展开就完成了任意初值的数学证明。

## 3. 空间线性色散及四阶截断的高频问题

在零场附近线性化，令 $q_j\propto\exp[i(kx_j-\Omega_a t)]$。这里 $k=2\pi n/L=n$，格点模式处于第一 Brillouin 区内。直接代入式 (5) 得

$$
\boxed{\Omega_a(k)^2=\omega_0^2+
\frac{4v^2}{a^2}\sin^2\!\left(\frac{ka}{2}\right)}.
\tag{12}
$$

连续与修正模型分别给出

$$
\Omega_0(k)^2=\omega_0^2+v^2k^2,
\qquad
\Omega_{(2)}(k)^2=\omega_0^2+v^2k^2-\frac{a^2v^2k^4}{12}.
\tag{13}
$$

展开精确格点式可核对

$$
\Omega_a^2=\Omega_0^2-\frac{a^2v^2k^4}{12}
+\frac{a^4v^2k^6}{360}+O(a^6k^8).
\tag{14}
$$

**式 (10) 不能作为全部波数有效的连续理论。** 任意固定 $a>0$ 下，$k\to\infty$ 时 $\Omega_{(2)}^2\to-\infty$，导致指数增长，式 (11) 的空间能量也无下界。对 $v>0$，其线性不稳定阈值为

$$
k^2>\frac{6}{a^2}\left[1+\sqrt{1+\frac{a^2\omega_0^2}{3v^2}}\right].
\tag{15}
$$

精确格点式 (12) 在本参数下始终为正；因此这个高频病态属于截断展开，并非原格点不稳定。阈值约 $\sqrt{12}/a$ 虽高于一维格点 Nyquist 波数 $\pi/a$，也不能据此把靠近 Nyquist 的模式称作长波，因为那里 $|ka|$ 并不小。

展开的受控精度要求实际占有的物理波段满足$|ka|\ll1$；数值截止$K_{\rm num}$则用于容纳尾部模式并排除不稳定的高频，两者不应混为一谈。本轮取$K_{\rm num}=24$，最粗网格的$K_{\rm num}a=2.356$并不远小于1，因此不能声称整条保留频带都有长波精度。应以实际尾谱、截止收敛和格距收敛评估这个初值的有效性。若做 Fourier–Galerkin 比较，令下式$K=K_{\rm num}$并明确采用投影 $P_K$：

$$
\phi_{tt}=v^2\phi_{xx}-\omega_0^2\phi
-gP_K(\phi^3)+\frac{a^2v^2}{12}\phi_{xxxx},
\qquad \phi=P_K\phi.
\tag{16}
$$

非线性会产生更高谐波，因此不能只检查初始模式低频。需要监测尾部谱能量，说明投影误差；用更高分辨率检验参考解，并避免把 Fourier 混叠误差误认为格点物理。对保留模式的投影与对原格点状态的观测滤波也是两种不同操作。

## 4. Störmer–Verlet、辛性与精确线性相位

定义 $U_a(q)=H_a(q,0)$。固定时间步长 $\Delta t$ 的 kick–drift–kick 更新为

$$
p^{n+1/2}=p^n-\frac{\Delta t}{2}\nabla U_a(q^n),
\quad q^{n+1}=q^n+\frac{\Delta t}{a}p^{n+1/2},
\quad p^{n+1}=p^{n+1/2}-\frac{\Delta t}{2}\nabla U_a(q^{n+1}).
\tag{17}
$$

每个 kick 和 drift 都是子 Hamilton 流，所以其组合保持式 (6) 的离散辛结构。该方法二阶、时间可逆。相关结构与修正方程的标准背景见 [Hairer、Lubich 与 Wanner，Störmer–Verlet 综述](https://www.unige.ch/~hairer/preprints/gniverlet.html)；下面的单模公式由式 (17) 直接推导。

对于 $H=p^2/(2a)+a\Omega_a^2q^2/2$，令 $z=\Delta t\Omega_a$，则

$$
\binom{q^{n+1}}{p^{n+1}}
=\begin{pmatrix}
1-z^2/2&\Delta t/a\\
-a\Delta t\Omega_a^2(1-z^2/4)&1-z^2/2
\end{pmatrix}
\binom{q^n}{p^n}.
\tag{18}
$$

此矩阵行列式为一，其特征相位满足 $\cos\theta=1-z^2/2$，因此在稳定振荡分支

$$
\boxed{\Omega_{\Delta t}(k)=\frac{2}{\Delta t}
\arcsin\!\left(\frac{\Delta t\,\Omega_a(k)}2\right)},
\qquad 0<\Delta t\Omega_a(k)<2.
\tag{19}
$$

严格小于二对一般初值的有界振荡是必要的；等号通常出现带线性增长的 Jordan 情况。$\Omega_a=0$ 是自由粒子，非零初速度导致线性位移，不能当作有界振子。本基准 $\omega_0=1$ 排除了这一零模例外。

式 (19) 对线性模式精确；它不直接给出有限振幅非线性解的频率。非线性系统的瞬时切线矩阵包含 $3gq_j^2$，冻结系数步长检查可用于数值诊断，但不能替代完整的非线性稳定性和收敛检查。

展开式 (19) 得

$$
\Omega_{\Delta t}=\Omega_a+\frac{\Delta t^2\Omega_a^3}{24}
+O(\Delta t^4\Omega_a^5).
\tag{20}
$$

而固定低波数的空间频率为

$$
\Omega_a=\Omega_0-\frac{a^2v^2k^4}{24\Omega_0}+O(a^4).
\tag{21}
$$

空间与时间离散的首个相位偏差符号相反，可能偶然抵消。因此必须分别加密网格与时间步长；一次轨迹重合不足以确认两种误差都小。

## 5. 能量恒定与能量数值误差的区别

本轮 $v,\omega_0,g$ 均固定，精确格点流满足

$$
\frac{dH_a}{dt}=\{H_a,H_a\}=0.
\tag{22}
$$

这与之前的非自治驱动场任务不同；本轮没有外部做功项。式 (17) 的辛性不意味着浮点轨迹严格保持原始 $H_a$，应实际记录 $[H_a(t)-H_a(0)]/H_a(0)$。对单个线性模式，它精确保留一个不同的二次型

$$
I_{\Delta t}=\frac{p^2}{2a}
+\frac{a\Omega_a^2}{2}\left(1-\frac{z^2}{4}\right)q^2,
\tag{23}
$$

而非一般情况下的物理模式能量。对非线性系统，修正 Hamilton 量和长期近似守能有相应条件，不能由一个有限时间实验宣称全局保证。

本初值具有一个可直接核对的能量归一化。对 $N>4$，离散初始能量精确为

$$
H_a(0)=\frac{LA^2}{4}
\left[\omega_0^2+\frac{4v^2}{a^2}\sin^2(a/2)\right]
+\frac{3LgA^4}{32}.
\tag{24}
$$

连续初始能量为

$$
H_0(0)=\frac{LA^2}{4}(\omega_0^2+v^2)
+\frac{3LgA^4}{32}=\frac{67\pi}{256}.
\tag{25}
$$

两者的首差为 $H_a(0)-H_0(0)=-Lv^2A^2a^2/48+O(a^4)$，这是不同空间模型的能量差，不是数值积分造成的能量漂移。

## 6. 误差预测与可检验的非线性信号

在本有限时间区间、足够光滑并控制未解析高频的前提下，标准连续模型与数值格点轨迹的预期差异为

$$
E_0=O(a^2)+O(\Delta t^2)+E_{\rm ref}+E_{\rm project},
\tag{26}
$$

其中后两项分别是参考求解误差和投影／截断误差。对共同固定的低频带，加入式 (10) 的首修正后预期为

$$
E_{(2)}=O(a^4)+O(\Delta t^2)+E_{\rm ref}+E_{\rm project}.
\tag{27}
$$

这些阶数不对 $k\sim1/a$ 一致成立，也不保证对任意长时间保持同一误差常数。若取 $\Delta t\propto a^2$，时间误差为 $O(a^4)$：足够分辨标准模型的二阶空间差异，却与首修正的四阶目标同阶。验证首修正的空间机制，还需独立缩小时步、估计时间误差或使用足够精确的格点 ODE 参考。

从初始 $\cos^3x=(3\cos x+\cos3x)/4$ 得到

$$
\phi_{tt}(x,0)=
-\left[A(v^2+\omega_0^2)+\frac{3gA^3}{4}\right]\cos x
-\frac{gA^3}{4}\cos3x.
\tag{28}
$$

若 $C_3(t)$ 定义为实余弦展开中 $\cos3x$ 的系数，则

$$
C_3(t)=-\frac{gA^3}{8}t^2+O(t^4)
=-\frac{t^2}{64}+O(t^4).
\tag{29}
$$

复 Fourier 系数 $\widehat\phi_{\pm3}$ 各为 $C_3/2$；对照时不要混用这两个归一化。初始的三次谐波来自四次势导致的模式耦合，在 $g=0$ 的独立线性对照中不存在。看到该谐波支持实现确实含有相应非线性，不能将其视为新物理或黎曼结构的证据。

## 7. 本基准能支持的结论

辛结构不唯一决定 Hamilton 量；本报告明确选择了式 (3)，再考察它的后果。周期格点的连续极限也不等于统计粗粒化：本轮没有引入热平衡、熵、耗散、随机闭合或量子化，尚未从微观粒子推导热力学或流体理论。

若实验满足误差预测，合适结论是：**对所选平滑初值和有限时间，局部非线性正则格点与连续场在长波段形成了可复现、误差可控制的对应，首个格距修正能进一步解释其差异。** 它为未来引入结构参数、非自治驱动或算术输入提供一个经过检查的基准；不表明真实世界的 Hamilton 量或宏观定律已从黎曼零点获得推导。
