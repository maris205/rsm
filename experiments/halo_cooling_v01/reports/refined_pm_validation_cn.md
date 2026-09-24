# 细化 PM 网格的独立验证

状态：**FAIL，102/103 项通过**。
第 1 次记录；历次尝试及失败保留。

实现检查：**PASS，95/95 项通过**；物理离散诊断：**FAIL，7/8 项通过**。

冻结协议 SHA256：`f3a0cc5abf21a4374439db756cbf8c4a3103bcabfa2bee6f9bd1921b28e608f2`。

本脚本只检查力算子与初始条件，不执行晕坍缩生产计算。RFFT 与旧 fullFFT、显式 DFT 的相符检查，验证的是声明的离散算子；它们不代表连续物理误差达到机器精度。

质量、自力、净力和规则粒子处力使用预先规定阈值。补偿球的连续力误差、粒子格相位误差单独记录，不据结果重设精度门槛。

每项检查以 scope 区分 implementation_protocol 与 physical_discretization_diagnostic。总状态保留所有失败；生产脚本仅以顶层 implementation_failures 是否为空判断力算法实现检查是否阻断。物理离散诊断的失败仍是失败，必须与随后六个固定配置的收敛结果一并解释。算法实现通过不能视为物理收敛。

| 检查 | 类别 | 状态 | 数值或说明 |
| --- | --- | --- | --- |
| rfft_fullfft_mesh:nf=8 | implementation_protocol | PASS | 9.09686e-17 |
| rfft_fullfft_particles:nf=8 | implementation_protocol | PASS | 9.477795e-17 |
| random_raw_mass:nf=8 | implementation_protocol | PASS | 0 |
| random_net_force:nf=8 | implementation_protocol | PASS | 2.3216964e-17 |
| single_self_force:nf=8:case=0 | implementation_protocol | PASS | 0 |
| single_self_force:nf=8:case=1 | implementation_protocol | PASS | 1.6653345e-16 |
| single_self_force:nf=8:case=2 | implementation_protocol | PASS | 5.5511151e-17 |
| rfft_fullfft_mesh:nf=16 | implementation_protocol | PASS | 1.7285439e-16 |
| rfft_fullfft_particles:nf=16 | implementation_protocol | PASS | 1.7483713e-16 |
| random_raw_mass:nf=16 | implementation_protocol | PASS | 0 |
| random_net_force:nf=16 | implementation_protocol | PASS | 2.7605861e-17 |
| single_self_force:nf=16:case=0 | implementation_protocol | PASS | 0 |
| single_self_force:nf=16:case=1 | implementation_protocol | PASS | 8.8817842e-16 |
| single_self_force:nf=16:case=2 | implementation_protocol | PASS | 2.220446e-16 |
| rfft_fullfft_mesh:nf=32 | implementation_protocol | PASS | 2.0496023e-16 |
| rfft_fullfft_particles:nf=32 | implementation_protocol | PASS | 2.5896436e-16 |
| random_raw_mass:nf=32 | implementation_protocol | PASS | 0 |
| random_net_force:nf=32 | implementation_protocol | PASS | 1.6484689e-17 |
| single_self_force:nf=32:case=0 | implementation_protocol | PASS | 0 |
| single_self_force:nf=32:case=1 | implementation_protocol | PASS | 1.7763568e-15 |
| single_self_force:nf=32:case=2 | implementation_protocol | PASS | 0 |
| explicit_DFT_mesh:nf=8 | implementation_protocol | PASS | 8.8324263e-16 |
| explicit_DFT_and_CIC_particles:nf=8 | implementation_protocol | PASS | 9.2029942e-16 |
| uniform_raw_mass:np=8:nf=8:shift=(0.0, 0.0, 0.0) | implementation_protocol | PASS | 0 |
| uniform_particle_force:np=8:nf=8:shift=(0.0, 0.0, 0.0) | implementation_protocol | PASS | 0 |
| uniform_raw_mass:np=8:nf=8:shift=(0.25, 0.375, 0.5) | implementation_protocol | PASS | 0 |
| uniform_particle_force:np=8:nf=8:shift=(0.25, 0.375, 0.5) | implementation_protocol | PASS | 0 |
| uniform_raw_mass:np=8:nf=8:shift=(0.5, 0.5, 0.5) | implementation_protocol | PASS | 0 |
| uniform_particle_force:np=8:nf=8:shift=(0.5, 0.5, 0.5) | implementation_protocol | PASS | 0 |
| uniform_raw_mass:np=8:nf=16:shift=(0.0, 0.0, 0.0) | implementation_protocol | PASS | 0 |
| uniform_particle_force:np=8:nf=16:shift=(0.0, 0.0, 0.0) | implementation_protocol | PASS | 0 |
| uniform_raw_mass:np=8:nf=16:shift=(0.25, 0.375, 0.5) | implementation_protocol | PASS | 0 |
| uniform_particle_force:np=8:nf=16:shift=(0.25, 0.375, 0.5) | implementation_protocol | PASS | 0 |
| uniform_raw_mass:np=8:nf=16:shift=(0.5, 0.5, 0.5) | implementation_protocol | PASS | 0 |
| uniform_particle_force:np=8:nf=16:shift=(0.5, 0.5, 0.5) | implementation_protocol | PASS | 0 |
| uniform_raw_mass:np=8:nf=32:shift=(0.0, 0.0, 0.0) | implementation_protocol | PASS | 0 |
| uniform_particle_force:np=8:nf=32:shift=(0.0, 0.0, 0.0) | implementation_protocol | PASS | 0 |
| uniform_raw_mass:np=8:nf=32:shift=(0.25, 0.375, 0.5) | implementation_protocol | PASS | 0 |
| uniform_particle_force:np=8:nf=32:shift=(0.25, 0.375, 0.5) | implementation_protocol | PASS | 8.6736174e-19 |
| uniform_raw_mass:np=8:nf=32:shift=(0.5, 0.5, 0.5) | implementation_protocol | PASS | 0 |
| uniform_particle_force:np=8:nf=32:shift=(0.5, 0.5, 0.5) | implementation_protocol | PASS | 1.7347235e-18 |
| phase_response_finite:nf=8:shift=(0.0, 0.0, 0.0):mode=(1, 0, 0):amplitude=2.5e-06 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=8:shift=(0.0, 0.0, 0.0):mode=(1, 0, 0):amplitude=1.25e-06 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=8:shift=(0.0, 0.0, 0.0):mode=(1, 1, 0):amplitude=2.5e-06 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=8:shift=(0.0, 0.0, 0.0):mode=(1, 1, 0):amplitude=1.25e-06 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=8:shift=(0.0, 0.0, 0.0):mode=(1, 1, 1):amplitude=2.5e-06 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=8:shift=(0.0, 0.0, 0.0):mode=(1, 1, 1):amplitude=1.25e-06 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=8:shift=(0.25, 0.375, 0.5):mode=(1, 0, 0):amplitude=2.5e-06 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=8:shift=(0.25, 0.375, 0.5):mode=(1, 0, 0):amplitude=1.25e-06 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=8:shift=(0.25, 0.375, 0.5):mode=(1, 1, 0):amplitude=2.5e-06 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=8:shift=(0.25, 0.375, 0.5):mode=(1, 1, 0):amplitude=1.25e-06 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=8:shift=(0.25, 0.375, 0.5):mode=(1, 1, 1):amplitude=2.5e-06 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=8:shift=(0.25, 0.375, 0.5):mode=(1, 1, 1):amplitude=1.25e-06 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=8:shift=(0.5, 0.5, 0.5):mode=(1, 0, 0):amplitude=2.5e-06 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=8:shift=(0.5, 0.5, 0.5):mode=(1, 0, 0):amplitude=1.25e-06 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=8:shift=(0.5, 0.5, 0.5):mode=(1, 1, 0):amplitude=2.5e-06 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=8:shift=(0.5, 0.5, 0.5):mode=(1, 1, 0):amplitude=1.25e-06 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=8:shift=(0.5, 0.5, 0.5):mode=(1, 1, 1):amplitude=2.5e-06 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=8:shift=(0.5, 0.5, 0.5):mode=(1, 1, 1):amplitude=1.25e-06 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=16:shift=(0.0, 0.0, 0.0):mode=(1, 0, 0):amplitude=1.25e-06 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=16:shift=(0.0, 0.0, 0.0):mode=(1, 0, 0):amplitude=6.25e-07 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=16:shift=(0.0, 0.0, 0.0):mode=(1, 1, 0):amplitude=1.25e-06 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=16:shift=(0.0, 0.0, 0.0):mode=(1, 1, 0):amplitude=6.25e-07 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=16:shift=(0.0, 0.0, 0.0):mode=(1, 1, 1):amplitude=1.25e-06 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=16:shift=(0.0, 0.0, 0.0):mode=(1, 1, 1):amplitude=6.25e-07 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=16:shift=(0.25, 0.375, 0.5):mode=(1, 0, 0):amplitude=1.25e-06 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=16:shift=(0.25, 0.375, 0.5):mode=(1, 0, 0):amplitude=6.25e-07 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=16:shift=(0.25, 0.375, 0.5):mode=(1, 1, 0):amplitude=1.25e-06 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=16:shift=(0.25, 0.375, 0.5):mode=(1, 1, 0):amplitude=6.25e-07 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=16:shift=(0.25, 0.375, 0.5):mode=(1, 1, 1):amplitude=1.25e-06 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=16:shift=(0.25, 0.375, 0.5):mode=(1, 1, 1):amplitude=6.25e-07 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=16:shift=(0.5, 0.5, 0.5):mode=(1, 0, 0):amplitude=1.25e-06 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=16:shift=(0.5, 0.5, 0.5):mode=(1, 0, 0):amplitude=6.25e-07 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=16:shift=(0.5, 0.5, 0.5):mode=(1, 1, 0):amplitude=1.25e-06 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=16:shift=(0.5, 0.5, 0.5):mode=(1, 1, 0):amplitude=6.25e-07 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=16:shift=(0.5, 0.5, 0.5):mode=(1, 1, 1):amplitude=1.25e-06 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=16:shift=(0.5, 0.5, 0.5):mode=(1, 1, 1):amplitude=6.25e-07 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=32:shift=(0.0, 0.0, 0.0):mode=(1, 0, 0):amplitude=6.25e-07 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=32:shift=(0.0, 0.0, 0.0):mode=(1, 0, 0):amplitude=3.125e-07 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=32:shift=(0.0, 0.0, 0.0):mode=(1, 1, 0):amplitude=6.25e-07 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=32:shift=(0.0, 0.0, 0.0):mode=(1, 1, 0):amplitude=3.125e-07 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=32:shift=(0.0, 0.0, 0.0):mode=(1, 1, 1):amplitude=6.25e-07 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=32:shift=(0.0, 0.0, 0.0):mode=(1, 1, 1):amplitude=3.125e-07 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=32:shift=(0.25, 0.375, 0.5):mode=(1, 0, 0):amplitude=6.25e-07 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=32:shift=(0.25, 0.375, 0.5):mode=(1, 0, 0):amplitude=3.125e-07 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=32:shift=(0.25, 0.375, 0.5):mode=(1, 1, 0):amplitude=6.25e-07 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=32:shift=(0.25, 0.375, 0.5):mode=(1, 1, 0):amplitude=3.125e-07 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=32:shift=(0.25, 0.375, 0.5):mode=(1, 1, 1):amplitude=6.25e-07 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=32:shift=(0.25, 0.375, 0.5):mode=(1, 1, 1):amplitude=3.125e-07 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=32:shift=(0.5, 0.5, 0.5):mode=(1, 0, 0):amplitude=6.25e-07 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=32:shift=(0.5, 0.5, 0.5):mode=(1, 0, 0):amplitude=3.125e-07 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=32:shift=(0.5, 0.5, 0.5):mode=(1, 1, 0):amplitude=6.25e-07 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=32:shift=(0.5, 0.5, 0.5):mode=(1, 1, 0):amplitude=3.125e-07 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=32:shift=(0.5, 0.5, 0.5):mode=(1, 1, 1):amplitude=6.25e-07 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| phase_response_finite:nf=32:shift=(0.5, 0.5, 0.5):mode=(1, 1, 1):amplitude=3.125e-07 | implementation_protocol | PASS | 布尔条件，详见 JSON |
| compensated_core_force_finite:np=32:nf=32:shift=(0.0, 0.0, 0.0) | physical_discretization_diagnostic | PASS | 布尔条件，详见 JSON |
| compensated_core_mean_force_inward:np=32:nf=32:shift=(0.0, 0.0, 0.0) | physical_discretization_diagnostic | PASS | -0.0023294162 |
| compensated_core_force_finite:np=32:nf=64:shift=(0.0, 0.0, 0.0) | physical_discretization_diagnostic | PASS | 布尔条件，详见 JSON |
| compensated_core_mean_force_inward:np=32:nf=64:shift=(0.0, 0.0, 0.0) | physical_discretization_diagnostic | PASS | -0.00091158547 |
| compensated_core_force_finite:np=32:nf=128:shift=(0.0, 0.0, 0.0) | physical_discretization_diagnostic | PASS | 布尔条件，详见 JSON |
| compensated_core_mean_force_inward:np=32:nf=128:shift=(0.0, 0.0, 0.0) | physical_discretization_diagnostic | FAIL | 0.00096490816 |
| compensated_core_force_finite:np=32:nf=128:shift=(0.25, 0.375, 0.5) | physical_discretization_diagnostic | PASS | 布尔条件，详见 JSON |
| compensated_core_mean_force_inward:np=32:nf=128:shift=(0.25, 0.375, 0.5) | physical_discretization_diagnostic | PASS | -0.0050146022 |

## 规则粒子格与力网格

规则粒子格在更细网格上通常产生非零沉积密度，并非加入了真实物理扰动。粒子处的力仍应由对称性抵消。二倍网格的 Bragg 模只有零或 Nyquist 分量，本实现的 Nyquist 安全梯度使基态网格力为零；四倍网格则可以在粒子之间存在网格力。

| 每方向粒子数 | 力网格 | 平移（力网格单元） | 网格密度 RMS | 网格力最大值 | 粒子力最大值 |
| --- | --- | --- | --- | --- | --- |
| 8 | 8 | [0.0, 0.0, 0.0] | 0 | 0 | 0 |
| 8 | 8 | [0.25, 0.375, 0.5] | 0 | 0 | 0 |
| 8 | 8 | [0.5, 0.5, 0.5] | 0 | 0 | 0 |
| 8 | 16 | [0.0, 0.0, 0.0] | 2.64575 | 0 | 0 |
| 8 | 16 | [0.25, 0.375, 0.5] | 0.572822 | 0 | 0 |
| 8 | 16 | [0.5, 0.5, 0.5] | 0 | 0 | 0 |
| 8 | 32 | [0.0, 0.0, 0.0] | 7.93725 | 0.21928 | 0 |
| 8 | 32 | [0.25, 0.375, 0.5] | 3.10242 | 0.0785828 | 8.67362e-19 |
| 8 | 32 | [0.5, 0.5, 0.5] | 2.64575 | 0.0464202 | 1.73472e-18 |

## CIC 相位与别名响应

此处未把粒子居中且与网格一一对应时的 sinc 公式用于不等分辨率。JSON 保留正负小扰动、两种幅度及三种整体平移下的力投影、非奇对称部分与横向分量。CIC 节点的分段导数和粒子格别名都可能导致离散响应差异；这些结果不能反拟合成解析增长参照后再用于自证。

## 补偿球的连续物理参照

均匀核心使用实际初始非线性过密度 δ=0.05；其严格映射为 r=q/(1+δ)^(1/3)，解析初始力为 g=−δr/3。固定选择初始拉格朗日半径 q/qL∈[0.25,0.75]，避免球心及核心边缘。有限性及平均径向力向内保留为物理诊断的通过／失败条件，不作为算法实现通过线；下表的离散力误差是诊断值，未预设它应达到高精度或随单独细化力网格单调改善。

| 粒子／力网格 | 平移 | 核心粒子数 | 对连续球力相对 L2 | 横向力相对 L2 | 个体力向内比例 | 映射最大误差 |
| --- | --- | --- | --- | --- | --- | --- |
| 32³ / 32³ | [0.0, 0.0, 0.0] | 880 | 0.0163873 | 5.24627e-05 | 1 | 5.55112e-17 |
| 32³ / 64³ | [0.0, 0.0, 0.0] | 880 | 0.631575 | 0.167789 | 0.854545 | 5.55112e-17 |
| 32³ / 128³ | [0.0, 0.0, 0.0] | 880 | 1.51682 | 0.485472 | 0.327273 | 5.55112e-17 |
| 32³ / 128³ | [0.25, 0.375, 0.5] | 880 | 1.65455 | 1.02524 | 1 | 5.55112e-17 |

## 失败记录

- compensated_core_mean_force_inward:np=32:nf=128:shift=(0.0, 0.0, 0.0)：{"name": "compensated_core_mean_force_inward:np=32:nf=128:shift=(0.0, 0.0, 0.0)", "passed": false, "mean_radial_component": 0.0009649081634239301, "particles_with_inward_force_fraction": 0.32727272727272727, "qualification": "Only the mean sign is a gate; individual-force accuracy is diagnostic", "scope": "physical_discretization_diagnostic"}

细化力网格不会增加粒子采样信息。本验证不等同于晕质量、密度剖面、气体冷却或星系形成收敛；后续仍需独立改变粒子数、力网格、时间步与网格相位。
