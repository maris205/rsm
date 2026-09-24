# 周期宇宙学PM的独立验证

状态：**PASS，43/43项通过**；第2次记录，历次失败不覆盖。

冻结协议：`protocol.md`，哈希`74a2cab351d53f98169b885542d23d18b4d425c6fd6b75b91b49f5c53422e2e1`。阈值未依据计算结果调整。

此实现采用网格节点j/N、粒子初值(j+1/2)/N、CIC沉积/回插、谱Poisson及Nyquist安全梯度。

| 检查 | 状态 | 结果 |
| --- | --- | --- |
| reject_nonfinite_power_input | PASS | 见JSON异常记录 |
| reject_nonfinite_mode_input | PASS | 见JSON异常记录 |
| reject_nonhermitian_mode_input | PASS | 见JSON异常记录 |
| uniform_force:n=16:shift=0.0 | PASS | 0 |
| uniform_raw_mass:n=16:shift=0.0 | PASS | 0 |
| uniform_force:n=16:shift=0.137 | PASS | 0 |
| uniform_raw_mass:n=16:shift=0.137 | PASS | 1.110223e-16 |
| random_raw_mass:n=16 | PASS | 1.110223e-16 |
| random_total_force:n=16 | PASS | 6.6896638e-17 |
| single_particle_self_force:n=16 | PASS | 3.5591467e-17 |
| analytic_mesh_fourier_force:n=16:mode=(1, 0, 0) | PASS | 2.8927011e-16 |
| analytic_mesh_fourier_force:n=16:mode=(1, 2, 0) | PASS | 8.7073147e-16 |
| analytic_mesh_fourier_force:n=16:mode=(2, 1, 3) | PASS | 1.6353985e-15 |
| uniform_force:n=32:shift=0.0 | PASS | 0 |
| uniform_raw_mass:n=32:shift=0.0 | PASS | 0 |
| uniform_force:n=32:shift=0.137 | PASS | 0 |
| uniform_raw_mass:n=32:shift=0.137 | PASS | 0 |
| random_raw_mass:n=32 | PASS | 0 |
| random_total_force:n=32 | PASS | 7.2885199e-17 |
| single_particle_self_force:n=32 | PASS | 3.6188537e-17 |
| analytic_mesh_fourier_force:n=32:mode=(1, 0, 0) | PASS | 2.6662686e-16 |
| analytic_mesh_fourier_force:n=32:mode=(1, 2, 0) | PASS | 7.6968318e-16 |
| analytic_mesh_fourier_force:n=32:mode=(2, 1, 3) | PASS | 1.5700524e-15 |
| eds_kick_factor:0.02:0.04 | PASS | 0 |
| eds_drift_factor:0.02:0.04 | PASS | 0 |
| eds_kick_factor:0.02:1.0 | PASS | 1.110223e-16 |
| eds_drift_factor:0.02:1.0 | PASS | 1.110223e-16 |
| eds_kick_factor:0.73:1.0 | PASS | 2.220446e-16 |
| eds_drift_factor:0.73:1.0 | PASS | 3.3306691e-16 |
| eds_kick_factor:0.1:0.02 | PASS | 2.220446e-16 |
| eds_drift_factor:0.1:0.02 | PASS | 2.220446e-16 |
| short_time_reversal_position_absolute | PASS | 0 |
| short_time_reversal_momentum_absolute | PASS | 6.6174449e-24 |
| plane_mass:n=16:steps=128 | PASS | 2.220446e-16 |
| plane_mass:n=16:steps=256 | PASS | 2.220446e-16 |
| plane_mass:n=16:steps=512 | PASS | 2.220446e-16 |
| fine_step_discrete_growth:n=16 | PASS | 1.8500207e-05 |
| KDK_time_difference_order:n=16 | PASS | 1.9995544 |
| plane_mass:n=32:steps=128 | PASS | 1.110223e-16 |
| plane_mass:n=32:steps=256 | PASS | 1.110223e-16 |
| plane_mass:n=32:steps=512 | PASS | 2.220446e-16 |
| fine_step_discrete_growth:n=32 | PASS | 1.9025994e-05 |
| KDK_time_difference_order:n=32 | PASS | 1.9995512 |

## EdS增长的两种参照

小幅平面模式初始密度幅度为2×10⁻⁶，a=.02→1，初速度采用连续增长支f=1。

有限CIC网格的独立解析力响应为F=sin(kh)/(kh)，r±=(−1±sqrt(1+24F))/4。因为初速度不是网格自身的纯增长支，解析参照保留增长与衰减两支。

时间收敛比较该网格解析增长；同时报告对连续EdS增长a/ai的空间偏差，不能把这两个误差混为一项。

| 网格 | 步数 | 增长 | 对网格解析相对差 | 对连续EdS相对差 |
| --- | --- | --- | --- | --- |
| 16³ | 128 | 47.35665278 | 0.000295911 | 0.0528669 |
| 16³ | 256 | 47.36716505 | 7.39961e-05 | 0.0526567 |
| 16³ | 512 | 47.36979393 | 1.85002e-05 | 0.0526041 |
| 32³ | 128 | 49.3129825 | 0.000304319 | 0.0137403 |
| 32³ | 256 | 49.32424015 | 7.60988e-05 | 0.0135152 |
| 32³ | 512 | 49.32705543 | 1.9026e-05 | 0.0134589 |

## 失败与修正记录

本次无未通过检查；既往尝试及其哈希保存在JSON中。

PM计算验证的是引入标准牛顿引力后的无碰撞物质接口，并未从黎曼结构导出引力。背景膨胀与有限网格平滑都保留；没有要求宇宙学物理能量在膨胀中恒定，也不因KDK名称宣称整个CIC粒子力已严格辛。
