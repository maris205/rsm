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

## 图4.3：从CMB约束初态到非线性场域

文件：`chapter04_timeline_cn.{pdf,svg,png}`。

上排采用CMB项目固定版本的实际存档状态，显示n=0、45、70、95的96³体积在z=48处的有符号中央切片；共用[-3.2,3.2]色标，不作逐面板归一化。第70步是真实积分快照。下面是完整三维体积在每一步的方差，比较四个共用前45步历史和q44/q45锚点的分支；“仅晚期常系数”只替换46—95步。纵轴为对数刻度，方差不等于熵或空间膨胀率。

初态受同一SMICA反演约束，非线性阶段由给定的双阱势及阻尼驱动。图内没有物理年龄、真实大爆炸奇点、星系或暴涨成功的标签，也没有用这张图替代CMB声学谱检验。此前旧Figure 16的全过程研究构想保留，物理宇宙史的映射作为下一步任务。

[最小来源快照与说明](../../reports/chapter04_timeline_sources/README.md)包含四张真实切片、每步体积统计及原结果/验证JSON副本，约326KB；其六份源输入已与CMB固定提交`f42db51748003fc83be064dade9bb52de5eb8319`逐字节校验。小NPZ SHA256：`ea7635e5a870e1abb0fffccc0614fbc3b50f69fdd39bd0554b79645af6d99351`。这些材料支持图形复现，不替代完整三维状态和原拟合输入。

复现命令：`python manuscript_cn/figures/make_chapter04_timeline_figure.py`。只读取本仓库快照，不需要相邻项目或重新演化。原项目52项验证报告作为历史来源复制，本次未重跑这52项检查。

## 图4.4：具有物理年龄的背景与物质引力演化

文件：`chapter04_cosmic_bridge_cn.{pdf,svg,png}`。来源为[cosmic_bridge_v01](../../experiments/cosmic_bridge_v01/README.md)新计算，使用同一协议下的标准背景与指数势正则时钟背景。上排将物理年龄、尺度因子、χ⁻²响应和共同初动量的线性增长接起来；下排展示主64³/256步时钟背景的z=9、1、0密度投影，共用对数色标。没有将SMICA反演场作为物质，也没有重新缩放初谱制造团块。

完整功率及时间/空间差异见实验目录的`figures/matter_power_diagnostics.pdf`。该图的CAMB z=0曲线是外部线性集合参考，未经过同一有限盒、带限和CIC窗口；差异不等于观测拟合误差或纯非线性修正。

复现：`python manuscript_cn/figures/make_chapter04_cosmic_bridge.py`。脚本只读取本仓库保存结果，输入输出身份与显示读数记录于`experiments/cosmic_bridge_v01/figures/plot_inputs.json`；不重新积分或拟合。

## 图4.5：受控塌缩的数值精度与条件冷却

文件：`chapter04_halo_cooling_cn.{pdf,png}`。来源为[halo_cooling_v01](../../experiments/halo_cooling_v01/README.md)的保存结果，四面板分别展示物理半径—宇宙年龄、六配置的事件及半径误差、质量—红移网格上的瞬时冷却／动力学时间比，以及两个连续球形事件的时间尺度。PM半径使用固定拉格朗日标签的中位缩放，是同调代理；与独立ODE的差异和未通过门限均保留，不以调整初态改善图片。

冷却网格中同一质量、红移的两背景微物理严格相同，因此使用一张共用热图；两条事件读数不同仅因达到阈值的红移不同。12个低于1e4 K、1个高于1e8 K的状态标为范围外，30个低温敏感状态及局部电离松弛警示另作标记。特征温度和均匀气体密度均是条件输入，图不表示已模拟气体热化、恒星或星系。

复现：`python manuscript_cn/figures/make_chapter04_halo_cooling.py`。只读保存数据，不重积分。输入输出哈希、失败标记与显示规则记录于`experiments/halo_cooling_v01/figures/plot_inputs.json`。正式数值及限制以正文和实验完整报告为准。

## 图4.6：匹配粒子与力网格尺度的受控检验

文件：`chapter04_matched_pm_cn.{pdf,png}`。来源为[pm_matched_scale_v01](../../experiments/pm_matched_scale_v01/README.md)冻结协议的两个参考背景配置，以及预先指定的旧64/256半步参照。左上比较物理半径与宇宙年龄，连续球显示至其密度比200事件；右上显示事件尺度因子的绝对相对误差及共同事件前保存采样的最大半径误差，2%／3%为新例门限。旧六配置完整结果仍见图4.5，不以这里的三例替代。

左下保留两个新例的全轨迹同调散布和10%门限。右下按原始拉格朗日半径q分箱，展示各自实际终态的径向缩放y=r/(a q)中位数及16—84百分位带；r为物理半径，a为宇宙尺度因子，各自终态a明确标注，空箱留空，浅蓝带表示固定统计核心。三例采用三项步长上限减半的控制，实际接受步不同。

128/128的最大半径误差3.126624%仍高于3%，两档事件差−1.833988%仍未过1%；所有失败保留。末态横向／径向力比没有随径向精度同步改善，完整角向诊断见实验报告。本图不是一般三维收敛、时钟微小差值或平方指数的确认。

复现：`python manuscript_cn/figures/make_chapter04_matched_pm.py`。只读本仓库保存结果，不重积分或拟合。输入、输出哈希与显示规则保存于`experiments/pm_matched_scale_v01/figures/plot_inputs.json`；图在完整稿中优先使用配对矢量PDF。
