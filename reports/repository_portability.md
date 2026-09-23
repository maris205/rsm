# 独立仓库的数据路径修复

2026-09-23。为使 `riemann_model` 独立克隆后可以读取必需输入，将五个小型源文件原样纳入实验目录。仅调整两个脚本的读取路径和随后运行时的来源路径记录；算法、参数、实验方案和既有数值结果未修改。本轮没有重跑模拟或拟合。

## 原样复制的输入

下表的目标路径相对本仓库根目录，原路径相对历史工作区。所有副本均与原文件逐字节相同，总计 26,271 字节。

| 本仓库副本 | 历史来源 | 字节数 | SHA-256 |
|---|---|---:|---|
| `experiments/arithmetic_residual_v01/data/raw/mathematical_reference_zeros.csv` | `riemann_clock/data/processed/mathematical_reference_zeros.csv` | 5460 | `24a20835b4c4028c6235d646d31279ad9d0bb49d98780062f1497fb240f5186c` |
| `experiments/alpha_recovery_v01/data/source/raw/King2012_tablea1.dat` | `riemann_constant/data/raw/King2012_tablea1.dat` | 16815 | `8bdec23fffb60800bba6afdcce4e077b639f4f743dc623207e80e38000d95964` |
| `experiments/alpha_recovery_v01/data/source/raw/King2012_provenance.json` | `riemann_constant/data/raw/King2012_provenance.json` | 1706 | `033faeb62c2984543bfda5d99d716dd1bf029de1c6367504a7e36461a512e8ab` |
| `experiments/alpha_recovery_v01/data/source/modern_measurements.csv` | `riemann_constant/data/modern_measurements.csv` | 427 | `43ff770b0f87232b38facfe51350df5cf8d1451480bff76999294f541bde9c96` |
| `experiments/alpha_recovery_v01/data/source/clock_constraints.csv` | `riemann_constant/data/clock_constraints.csv` | 1863 | `0f90032944f01da5138f5e4320a64ee05fbea3193f7357587a0b0d1c4c229929` |

零点参考表记录 mpmath 1.3.0、35 位十进制，含索引 1–81、4200、4201，不是区间认证。King 原表及其来源 JSON、现代测量和钟约束沿用既有来源与单位；它们的纳入不增加观测样本。前作 10,000 点缓存与 Odlyzko 表在修复前已经位于算术实验 `data/raw/`，无需依赖旁项目。

## 路径及来源记录

[prepare_inputs.py](../experiments/arithmetic_residual_v01/code/prepare_inputs.py) 从本实验 `data/raw/mathematical_reference_zeros.csv` 读取 83 点参考。后续重新生成的来源清单使用实验内相对路径，并另存 `original_workspace_path` 表明副本来源。

[run_recovery.py](../experiments/alpha_recovery_v01/code/run_recovery.py) 从本实验 `data/source/` 读取四个输入；后续生成的 `source_hashes` 键相对本实验根目录，方案哈希对应 `protocol.md`。不再要求克隆目录外有 `riemann_constant` 或 `riemann_clock`。

既有结果中的旧源路径和代码哈希描述原始运行，保留不变。此次路径修复的脚本哈希变化如下；不能把原运行记录误读为修改后脚本已完整重跑。

| 脚本 | 修复前 SHA-256（与既有结果一致） | 修复后 SHA-256 |
|---|---|---|
| `experiments/arithmetic_residual_v01/code/prepare_inputs.py` | `32c932cff382383b7c270c96334714fd53352b3db3b77d611ce445405ac481e5` | `1be1a2a367058a7cd478cab0479d7011e0f4cf15e789db57a706d0e0c6dc6ce8` |
| `experiments/alpha_recovery_v01/code/run_recovery.py` | `7f8f56274a8ab42f5dd8fd9b057998cf1cf6639e83f158394f2d253cadd7bfe4` | `2ef734151f5e3e6369bcd66b1498cdd1ddd169c4cefefb9b9a96657da73b91dd` |

## 本轮验证

- 两个脚本均通过 Python 语法解析，并以非主入口方式加载，未执行数值主程序。
- Alpha 输入函数从新副本生成的 293 条 King 设计与 1 条 ESPRESSO 设计，在红移、误差、目标、编号和仪器字段上逐项等于既有 `observation_design.csv`。
- 钟标准差与既有摘要完全相同，浮点值为 `2.5000000000000004`。初次检查误把浮点商要求为字面量 `2.5` 精确相等；随后改为与历史结果精确比较，并按原脚本容差比较 `2.5`，均通过，输入及算法未因此改动。
- 83 个零点参考与已携带 10,000 点缓存的最大 float64 绝对差为 0；五个副本的逐字节与哈希检查均通过。

以上核验只支持数据路径迁移的一致性。历史实验保留的匹配失败、精度诊断和物理解释边界不因仓库发布而改变。
