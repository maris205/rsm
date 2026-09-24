# 第4章图注与来源

## 图4.1：同一相互作用：从格点到长波场

文件：`chapter04_micro_macro_bridge_cn.{pdf,svg,png}`。

算术研究启发的逆对数平方响应及可选残差耦合作为候选输入，进入正则格点—时钟系统的共同 Hamilton 量。场端回复项与时钟端反馈均由同一耦合能求导。保持 $L$、$I$ 及模型系数固定，采用 $p_\ell=a\dot q_\ell$ 的正则归一化，在平滑长波条件下获得场方程及相同的积分反馈；加入算术残差时，图中的 $[M^2]'(R)$ 保留完整导数。下排分别列出空间格距修正、时间积分误差和仍需闭合的统计平均。连续描述保留共享均匀 $R$ 的全局性，并未建立局域驱动场，也未从黎曼结构唯一推得真实物理相互作用。

这是新增的方法示意，原第4章没有对应图片。图不展示新实验或拟合；依据 [自主时钟—格点动力学](../../manuscript_v01/chapters/03_autonomous_dynamics.md)、[连续有效场推导](../../manuscript_v01/chapters/04_continuum_limit.md) 与 [研究方向备忘](../../RESEARCH_DIRECTION.md) 绘制。长波方程为主阶描述；有限格距的首修正及其频带限制在下排与正文单列。算子展开本身不等于任意长时间的解收敛证明。

复现命令：`python manuscript_cn/figures/make_chapter04_figures.py`。
提供矢量 PDF、SVG 与 240 dpi PNG；沿用既有中文字体和配色。
SVG 在保存后移除行尾空白，便于版本管理检查。

## 图4.2：微观迭代的可计算宏观读数

文件：`chapter04_wave_spectrum_cn.{pdf,svg,png}`。

左图采用独立实现的45步模态递推，比较逆对数平方系数（p=2，A2=0.45，c0=10）与常系数（p=0，A0=0.040348613564280385）的模态功率传递。两条历史匹配累计长波传播尺度Σsqrt(c_n²)/sqrt(6)=3.690210980139839格点，没有匹配每个模态的精确相位；两者都产生峰谷。横轴λ为归一化三维六邻域算子的非负本征值，范围0—2。

右图展示长波近似λ2=k²/6及四阶修正λ4=k²/6−Σk_i⁴/72相对于精确本征值的误差，分别取轴向和体对角方向；横轴k为波矢模长，采用格点单位。图中二阶、四阶指相对误差的主阶。连续展开只在小k范围受控，大k曲线用于显示其边界，不将展开视为整个频带的物理定律。

来源为[预先固定的波动协议](../../experiments/macro_laws_v01/protocol_wave.md)、[独立代码](../../experiments/macro_laws_v01/code/verify_wave_transfer.py)及[曲线NPZ](../../experiments/macro_laws_v01/results/wave_transfer_curves.npz)。NPZ SHA256为`61c38ac59e6b8af538c5125df51383b24cd4f7c173395c2185f0d17c4ed98ab7`。无天空观测、无参数拟合，不能由峰谷断言p=2独有或优于p=0。

复现命令：`python manuscript_cn/figures/make_chapter04_macro_figure.py`。脚本只读已保存曲线，复用既有中文字体与配色，提供PDF/SVG矢量版及PNG预览；不重跑动力学。
