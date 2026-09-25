# 独立数值复核

本检查不导入主求解器；以宇宙时间和正则动量重新求解同一冻结模型，然后反解 N 将两个求解器的读数对齐。

- 状态：8 个案例全部通过。
- 阈值：轨迹与能量误差使用明确非零尺度归一，均要求小于 10⁻⁶。
- 检查：χ 增量、H、年龄、低能 α 阈值、标量能量、时钟速度，以及 ∆ρ + ∫3H Aχ̇²dt=0；CW 斜率和动能修正斜率另作五点差分复核。
- 这些检查验证两种数值实现一致，不验证新带电粒子的现实可行性，也不代表谱线观测支持模型。

| 案例 | 通过 | 最大误差 |
|---|---|---:|
| chi_probe | True | 9.653e-11 |
| chi_loop_1 | True | 9.680e-11 |
| chi_loop_100 | True | 1.062e-10 |
| chi_loop_10000 | True | 5.121e-10 |
| R_probe | True | 1.774e-10 |
| R_loop_1 | True | 1.735e-10 |
| R_loop_100 | True | 1.750e-10 |
| R_loop_10000 | True | 4.682e-10 |

复核方程使用 $\tau=H_{\rm ref}t$、$p=A\,d\chi/d\tau$：

$$
\frac{dN}{d\tau}=E,\qquad
\frac{d\chi}{d\tau}=\frac pA,\qquad
\frac{dp}{d\tau}=-3Ep+\frac{A_{,\chi}}2\left(\frac pA\right)^2-V_{,\chi}.
$$

$E$ 直接由流体密度和 $p^2/(2A)+V$ 的 Friedmann 方程求得；本实现没有复制主程序以 $N$ 为自变量的速度方程，也没有导入主程序。初始 $\chi_i,q_i,\tau_i$、树势幅度与 $\Lambda$ 从冻结输入读取，两求解器面对同一初值问题。CW 与二导数系数的解析依据另见[作用量审计](../../../reports/charged_threshold_action_audit_20260925.md)。

两条 $\ell=10^4$ 分支均独立复现了 $q<0$ 的时钟反转。其初始低能响应分别为 $-0.168392$ ppm（常 $\chi$ 动能）与 $-0.275815$ ppm（常 $R$ 动能）。这些是受限玩具模型的计算输出。复核没有把反转当作错误而排除，也没有为获得单调曲线修改参数。

材料范围：仅使用本项目冻结协议、解析模型与数值输入，无新观测或第三方私有数据，无外部数据传输。执行命令（仓库根目录）：

```bash
timeout 180 python experiments/charged_threshold_v01/code/independent_check.py > experiments/charged_threshold_v01/results/independent_run.log 2>&1
```

输入和实现 SHA256、逐项误差在 [independent_checks.json](../results/independent_checks.json)；独立轨迹在 [independent_trajectories.csv](../results/independent_trajectories.csv)，执行记录在 [independent_run.log](../results/independent_run.log)。运行退出码为 0，未出现求解失败或数值域终止。上述文字说明是运行后补充的审阅记录；表格及 JSON 由复核脚本生成。
