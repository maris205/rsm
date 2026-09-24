# 宇宙年龄与无碰撞物质聚类：最小物理接口

本实验把已有的指数势正则时钟候选接入辐射、无压物质和Λ背景，再用标准Newton引力推进粒子。它给出有物理单位的参考宇宙年龄及后期聚类计算，标准引力、外部初始物质谱和物理单位均明确作为新增接口。

**时钟只通过其能量改变背景H(a)。** Newton常数、粒子质量和Poisson力不随χ变化；χ⁻²属于另行假定的响应函数，没有被写成新的引力律。这里也没有把主稿固定盒的全局R重命名为宇宙尺度因子。

## 当前内容

- [执行前冻结的协议](protocol.md)：输入、运行矩阵、阈值及解释边界。
- [中文结果报告](reports/results_cn.md)：方法、年龄、匹配增长和PM结果。
- [背景JSON](results/background_summary.json)：109/109项检查通过；三组ε均为声明的计算案例。
- [独立背景与输入核查](reports/independent_bridge_checks.json)：另用Radau和不同状态变量，33/33项通过。
- [PM核心验证](reports/pm_validation_cn.md)：最终43/43项通过；首次40/40及输入保护后的尝试分别保存，不重复累计。
- [独立PM小诊断](reports/independent_pm_checks.json)：17/17项通过；检查Fourier/CIC接口，不重复运行生产矩阵。
- [保存结果的独立审计](reports/clustering_independent_audit.json)：347/347项通过；读取8次运行的数组并独立重算FFT/CIC读数，不是另一次完整PM积分。
- [PM生产汇总](results/clustering_summary.json)：ε=0与10⁻⁴的4组配对共8次运行完成，24/24项有限性、质量和总力检查通过。按频箱平均波数≤.3 h/Mpc选出的两个低k箱，主例功率分别降低约0.0265%和0.0287%；第二箱最高贡献模式约.3078 h/Mpc，不能称其中所有模式均≤.3。绝对谱仍有明显空间分辨率差，详见结果报告及[执行说明](reports/execution_notes.md)。

## 复算

以下命令均从**仓库根目录riemann_model**执行。背景、PM和独立核查依赖Python、NumPy与SciPy；保存的背景运行环境为Python 3.12.3、NumPy 2.4.4、SciPy 1.16.1。已有`inputs/`可直接用于复算，常规运行无需安装CAMB。

先生成背景及检查结果：

```bash
python experiments/cosmic_bridge_v01/code/background.py
```

然后执行独立背景核查和PM核心验证：

```bash
python experiments/cosmic_bridge_v01/code/independent_bridge_check.py
python experiments/cosmic_bridge_v01/code/independent_pm_checks.py
python experiments/cosmic_bridge_v01/code/validate_pm.py
```

运行冻结的8组PM矩阵：

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
python experiments/cosmic_bridge_v01/code/run_clustering.py \
  --reference experiments/cosmic_bridge_v01/results/background_eps0.csv \
  --clock experiments/cosmic_bridge_v01/results/background_eps1e-4.csv \
  --spectrum experiments/cosmic_bridge_v01/inputs/reference_linear_spectrum.csv \
  --f-initial 0.9903491493236205 \
  --workers 4
```

`--f-initial`是参考背景在a=.02的正则增长率。两个背景使用共同的实际初动量，所以时钟背景的匹配增长初值是f_ref E_ref/E_clock；不能给两个模型分别设相同的f却称为共同初动量。CAMB的z=49功率直接进入初始位移，不再用今日谱或增长比重缩幅度。

4个工人只并行执行相互独立的运行，共享只读初态，不改变冻结的物理运行矩阵。记录的并行wall时间用于说明本次资源消耗，不作为求解器速度基准。

生产结束后可独立核查保存的数组和汇总：

```bash
python experiments/cosmic_bridge_v01/code/audit_saved_clustering.py
```

### 可选：重新生成外部CAMB输入

输入已保存时应跳过此步。若需要从声明参数重生成，可在独立环境安装固定CAMB版本，避免改动背景/PM运行环境：

```bash
python -m venv /tmp/rsm-camb
/tmp/rsm-camb/bin/python -m pip install camb==2.0.4
OMP_NUM_THREADS=1 /tmp/rsm-camb/bin/python \
  experiments/cosmic_bridge_v01/code/make_reference_spectrum.py \
  --protocol experiments/cosmic_bridge_v01/protocol.md \
  --protocol-sha256 74a2cab351d53f98169b885542d23d18b4d425c6fd6b75b91b49f5c53422e2e1
```

[输入说明](REFERENCE_SPECTRUM_INPUT.md)、[实际参数](inputs/camb_parameters.txt)和[输入元数据](inputs/reference_spectrum_metadata.json)记录外部谱的单位、红移、版本、辐射密度差以及数组哈希。它是标准线性传递输入，不是从黎曼结构导出的原初谱，也不是对SMICA相位的反演。

## 输出和接口

| 文件 | 内容 |
| --- | --- |
| `results/background_eps0.csv` | 无时钟参考背景，a=.02至1，共4097点 |
| `results/background_eps1e-4.csv` | 主时钟例，含相同实际动量的匹配增长 |
| `results/background_eps1e-2.csv` | 较大ε的背景敏感性诊断，不参与主PM配对 |
| `results/background_grid.csv` | 三组背景从a=10⁻⁷开始的完整采样表 |
| `results/background_summary.json` | 初始条件、年龄、闭合Λ、所有检查及失败 |
| `results/pm_validation.json` | 最新PM验证；历次记录保留为`pm_validation_attempt_*.json` |
| `results/initial_modes.npz` | 共同带限初始Fourier系数 |
| `results/reference_n*_s*.npz`、`clock_n*_s*.npz` | 五个尺度因子的密度、投影和功率读数 |
| `results/clustering_summary.json` | 全部生产诊断、配对响应、时间/空间比较及哈希 |

背景JSON以`models["eps0"]`、`models["eps1e-4"]`、`models["eps1e-2"]`索引；每项的`initial_pm`含E、f_regular、f_matched。Python类接口见[background.py](code/background.py)的docstring。

## 解释范围

H0、物质/辐射密度、ε和势的时间锚点都是声明输入；调整Λ闭合今日H0也不是预测H0。共同z=49初谱仅隔离后期背景差异，没有计算标量对早期传递的改变。平滑辐射增长支只规定初速度，不等同完整的辐射期扰动演化。

主PM盒为50 Mpc/h，64³粒子和网格，单元宽0.78125 Mpc/h。没有气体、冷却、恒星形成、反馈或星系识别；输出称为无碰撞物质聚类，不把团块称作已经分辨的暗物质晕或星系。图像差异、同一随机种子的配对差以及数值检查均不构成观测发现。
