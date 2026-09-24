# 外部CAMB线性参考谱输入

本输入提供标准线性物质传递的起点与参考终点，不是RSM从黎曼结构导出的原初谱、CMB角谱或物质谱，也不复用SMICA反演的相位。后续两种背景共用z=49的输入，实现对后期作用的配对比较；标量对更早传递的影响未计算。

## 固定参数与单位

采用H0=67.4 km/s/Mpc、h=.674、Ωb0=.049、Ωc0=.266、Ωm0=.315、平直空间、w=−1、wa=0。CAMB接收`ombh2=.049*.674²`和`omch2=.266*.674²`。中微子总质量及有质量种数均为0，N_eff=3.046；TCMB=2.7255 K、As=2.1e−9、ns=.965、τreio=.054，标量参考波数0.05 Mpc⁻¹，无running或张量初谱。氦丰度采用该CAMB版本的BBN一致性默认值，实际结果及全部底层参数另外保存。

输出2049个对数等距点，k∈[.001,10] h/Mpc；功率单位为(Mpc/h)³，分别保存z=49和z=0的线性总物质自功率。密度变量指定`delta_tot`，此无质量中微子设置下为CDM与重子按密度加权的同步规范密度对比；不能据此把接近视界的模式直接当成Newton粒子读数。

CAMB的内部`set_matter_power(kmax=...)`使用物理Mpc⁻¹，因此计算上限设为`1.2*h*10=8.088 Mpc⁻¹`。输出的最大k/h则为10 h/Mpc。脚本核对原生支持域覆盖，另比较两套单位的同一原生谱，验证`k_physical=h*k_table`和`P_physical=P_table/h³`。没有高k外推。计算取AccuracyBoost=2、lAccuracyBoost=2、k_per_logint=50、关闭CMB角谱/透镜及所有非线性功率修正；这些是固定数值设置，不是对数据调参。

参数与单位约定依据[CAMB参数文档](https://camb.readthedocs.io/en/latest/model.html)、[功率谱及sigma8接口](https://camb.readthedocs.io/en/latest/results.html)和[密度变量定义](https://camb.readthedocs.io/en/latest/transfer_variables.html)。实际安装版本在执行输出中记录，不以文档版本代替运行版本。

## 输出文件

- `inputs/reference_linear_spectrum.csv`：三列`k_h_mpc,pk_z49_mpc_over_h3,pk_z0_mpc_over_h3`。
- `inputs/reference_linear_spectrum.npz`：同名三数组，另有`redshifts=[49,0]`、对应`sigma8`及来源/单位元数据。
- `inputs/camb_parameters.txt`：实际CAMB参数的完整文本表示，保留版本相关默认值。
- `inputs/reference_spectrum_metadata.json`：协议/脚本/本说明SHA256、实际CAMB与NumPy版本、σ8、CAMB的Ωγ0与无质量Ων0、相对背景Ωr0=9.2e−5的差、单位检查、数组顺序及三个输出文件哈希。

`get_sigma8()`的红移顺序与功率表可能相反，脚本按明确红移寻找对应行，不凭位置猜测。输出中的`fsigma8/sigma8`只作CAMB参考诊断；本轮冻结协议不生成k依赖f_i(k)，并且初始粒子动量统一采用背景模块的平滑辐射正则增长约定。CAMB自身标准传递与后续近似背景之间的差别保留在元数据，不通过重归一化谱隐藏。

## 执行约定

生成前必须显式传入已冻结协议的SHA256；不匹配时在调用CAMB之前退出。当前冻结协议不允许`--growth-diagnostic`。使用隔离安装环境，不自动安装软件：

```bash
PYTHONPATH=/tmp/rsm-cosmic-bridge-deps OMP_NUM_THREADS=1 python experiments/cosmic_bridge_v01/code/make_reference_spectrum.py --protocol-sha256 74a2cab351d53f98169b885542d23d18b4d425c6fd6b75b91b49f5c53422e2e1
```

所有谱、单位、覆盖和背景一致性检查均记录通过或失败；失败输出不用于后续生产计算。该输入不拟合观测，不反推H0或调节σ8以让结构更明显。
