# 三维线性格点的微观—宏观连接：小型确定性验证

运行状态：**PASS，28/28 项通过**。协议先记录后运行；本地检查不构成外部注册或同行评审。

## 模型与检验范围

这是与主 RSM 自主时钟模型分开的三维线性候选。周期六邻居算子除以6，主递推的初态为 q[-1] = q[0]，执行45步；不加载天空数据、不作观测拟合、没有阻尼、非线性或裁剪。单位均为格点/迭代单位。

平方对数分支 A2=0.45、c0=10；常系数分支 A0=0.0403486135643 匹配累计长波距离代理量，双方均为 3.69021098014 个格点。该条件不匹配每个 Fourier 模态的精确相位。

坐标/FFT检查网格为12³；完整协方差检查网格为4³，直接演化确定性协方差平方根的全部列，无大样本ensemble及径向平均。

## 检查结果

| 检查 | 状态 | 数值或说明 |
| --- | --- | --- |
| matched_long_wave_distance | PASS | 8.88178e-16 < 1e-13 |
| coordinate_fft_relative_l2:p=0 | PASS | 9.51218e-15 < 1e-11 |
| linear_mean_conservation:p=0 | PASS | 2.11636e-16 < 1e-13 |
| coordinate_fft_relative_l2:p=2 | PASS | 3.84983e-15 < 1e-11 |
| linear_mean_conservation:p=2 | PASS | 1.59595e-16 < 1e-13 |
| coordinate_fourier_eigenvalue | PASS | 1.20873e-14 < 1e-13 |
| constant_transfer_closed_form | PASS | 2.43972e-14 < 2e-11 |
| single_mode_symplectic | PASS | 0 < 1e-14 |
| 45_step_mode_symplectic | PASS | 2.86438e-14 < 1e-10 |
| 45_step_mode_determinant | PASS | 5.68434e-14 < 1e-10 |
| selected_frozen_oscillatory_upper_bound | PASS | max(c²λ) = 0.156524 < 4；仅冻结系数条件 |
| long_wave_coefficient_one_sixth | PASS | 8.33331e-06 < 1e-05 |
| quartic_direction_coefficient | PASS | 4.62959e-08 < 1e-07 |
| relative_long_wave_error_order_2 | PASS | 2.000541, 2.000271, 2.000180 |
| relative_corrected_error_order_4 | PASS | 4.000709, 4.000354, 4.000236 |
| sine_cosine_eigenvalue_identity | PASS | 4.44089e-16 < 1e-14 |
| dense_laplacian_symmetry | PASS | 0 < 1e-13 |
| dense_eigendecomposition | PASS | 8.32667e-16 < 1e-12 |
| coordinate_modal_operator:p=0 | PASS | 7.22525e-15 < 1e-10 |
| complete_covariance_propagation:p=0 | PASS | 9.26727e-15 < 1e-10 |
| coordinate_modal_operator:p=2 | PASS | 4.59285e-15 < 1e-10 |
| complete_covariance_propagation:p=2 | PASS | 5.27709e-15 < 1e-10 |
| retained_mode_subspace_nonempty | PASS | 保留 64/64 个模式 |
| prior_absorption_target_covariance:p=0 | PASS | 1.79814e-14 < 1e-10 |
| prior_absorption_target_covariance:p=2 | PASS | 1.19569e-14 < 1e-10 |
| prior_absorption_equal_final_covariances | PASS | 2.15744e-14 < 1e-10 |
| finite_propagation_sum_diagnostics | PASS | 有限且为正；只作有限求和诊断 |
| constant_schedule_sum_ratio | PASS | 0 < 2e-14 |

## 初谱退化与有限求和

按事先固定的双方 |T|>0.1 规则保留 64/64 个模式，排除 0 个；没有除以传递零点。各历史输入功率由同一目标功率除以该历史 T² 给出。保留子空间上可得到同一完整终态协方差；这显示初态自由度，不能解释成观测对任一历史的支持。

有限和与渐近主项之比：

| p | N=45 | N=100 | N=1,000 | N=10,000 | N=100,000 | N=1,000,000 |
| --- | --- | --- | --- | --- | --- | --- |
| 0 | 1.00000000 | 1.00000000 | 1.00000000 | 1.00000000 | 1.00000000 | 1.00000000 |
| 1 | 1.06492288 | 1.09019372 | 1.08752707 | 1.06592646 | 1.05082171 | 1.04109973 |
| 2 | 1.13986239 | 1.19735623 | 1.19329651 | 1.14290745 | 1.10804307 | 1.08620208 |
| 3 | 1.22660915 | 1.32563498 | 1.32335429 | 1.23440010 | 1.17325027 | 1.13605123 |

此处所有p单独固定A=0.45，A在比值中抵消；不是45步匹配设计的新增模型比较。有限和不能证明渐近等价或给出有限声学视界。p=2的数学主项为N/ln N，仍发散。

## 解释边界

- 长波展开的1/6与四阶方向项来自三维归一化邻居几何；现有参数尚未定标为物理光速、格距或宇宙年龄。
- 初始相等位置是零后向差分，不是连续零初速度的二阶起步；常系数闭式保留(N+1/2)相位及分母。
- 辛性与固定模态振荡条件没有证明所有时变系统的稳定性或自主能量守恒。
- 平均场保持是本线性周期模型的性质，尚未识别为物质质量守恒。
- 功率是场的二阶统计。不同历史都可产生振荡特征，不能由峰形单独推出平方对数规律、原子能级或CMB光子能谱。
- 协方差吸收允许每个保留模态自由调整初始方差；固定低维初谱、传递零点或额外共同观测可能限制此自由度。

## 失败记录与来源

本次协议内检查无失败。全部数值、阈值与条件保存在JSON中；没有在运行后调整阈值。

来源提交：`f42db51748003fc83be064dade9bb52de5eb8319`。实现前以 `git show` 核对三份来源SHA256；本脚本不导入相邻项目代码。

协议 SHA256：`57108d82f7676be910a8691b80c242e10f4e41862a04f9580994cac36c01fb56`。

脚本 SHA256：`fb220f741b2846c10a44cd2202f00c43833041183fbb0e62816aec4356f7352a`。

曲线 NPZ SHA256：`61c38ac59e6b8af538c5125df51383b24cd4f7c173395c2185f0d17c4ed98ab7`。

曲线文件同时保存轴、方向、单位、参数JSON；读取时无需允许pickle。完整协方差检查是有限矩阵恒等式核对，不是新观测或物理机制发现。
