# 自主驱动候选 v0.2：背景延拓与能量核验

**已完成：** 指数势标量在尘埃＋Λ背景上的正向积分，纯物质解析基准、能量账目、数值收敛和另一变量／算法的独立复核。29项主检查及24项独立比较通过。这里没有读取观测响应或执行拟合。

本轮为[v0.2主稿](../../draft_v02/action_framework_v02.md)提供一个可复现的例子：逆对数平方历史可以在标度极限成立，但同一个作用量延拓到晚期后会产生偏离。主实现和独立实现的结果一致。

## 主要结果

- 纯物质标度支恢复$\chi=\ln(t/t_*)$、$\Omega_\chi=3\epsilon/4$与$w_\chi=0$。
- 加入Λ后，在示例$\epsilon=(F_\chi/M_{\rm Pl})^2=10^{-4}$中，今天$\dot\chi_0t_0=0.817013$；固定同一当前漂移，并用同一$t(z),t_0,t_*$比较，在$z=4.2$得到$\mathcal T_\chi/\mathcal T_{\log}=1.167431$。
- 约16.7%是两个模型的漂移归一化传递之差。若示例$D_0=2.5\times10^{-19}$年$^{-1}$，相应α预测只相差约$0.00133$ ppm；既不是α已测得16.7%的变化，也不是拟合效果提高。
- 独立实现更换积分变量、Friedmann计算式、传递公式和积分算法，24项输出一致性比较通过。它读取了主实现的Λ闭合根，未独立重求此根。

![驱动背景与观测传递](figures/background_completion.png)

## 复现与阅读

运行环境：Python 3、NumPy、SciPy、Matplotlib；结果JSON记录实际版本。以下命令在上级工作区 `dsc_world` 执行，顺序先主后独立：

```bash
python riemann_model/experiments/action_completion_v02/code/run_action_background.py
python riemann_model/experiments/action_completion_v02/code/independent_background_check.py
```

命令会覆盖本实验目录中的对应数值结果和图。

| 文件 | 用途 |
| --- | --- |
| [实验设定](protocol.md) | 作用量、变量、初值、比较方法与范围。 |
| [中文结果报告](reports/numerical_results_cn.md) | 方程、数值表、能量收支、初值敏感性和限制。 |
| [独立复核](reports/independent_background_review.md) | 另一变量组与Radau积分；误差与独立性边界。 |
| [主代码](code/run_action_background.py) / [独立代码](code/independent_background_check.py) | 生成及复核结果。 |
| [主JSON](results/background_summary.json) / [独立JSON](results/independent_background_checks.json) | 全部检查值、软件版本及代码／输入哈希。 |
| [摘要CSV](results/background_summary.csv) / [红移网格](results/background_grid.csv) / [初值探测](results/initial_condition_probes.csv) | 机器可读结果。红移网格为3组参数各841点。 |
| [PNG图](figures/background_completion.png) / [PDF图](figures/background_completion.pdf) / [图注](figures/captions.md) | 便于查看、排版和引用。 |

## 模型边界

本轮忽略辐射和平均电磁／物质标量源，尘埃单独守恒，初始$a_i=10^{-3}$只是理想系统的数值起点。所取$H_0,\Omega_{m0},t_*,\epsilon,D_0$是示例输入，不是新测量或已允许参数区间。每例调整Λ以闭合今天的Friedmann方程，没有预测哈勃常数。初始速度扰动也会重新闭合Λ，不能称为固定所有模型参数的吸引子证明。

α对$\chi$的平方倒数响应仍是公设；势函数本身不唯一选择平方指数。物质源、第五力、环境及量子修正须另行计算。本轮也未改变[上一轮模拟恢复](../alpha_recovery_v01/README.md)对现有设计难以辨认平方指数的结论。
