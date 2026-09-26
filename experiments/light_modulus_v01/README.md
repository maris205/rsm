# 保留Higgs轻模：非线性动能、有限辅助场势与驻点检验

**原链的小端点交叉接触可以保留，但右端源对剩余轻模的作用没有受到链长抑制。在声明的刚性分支中，有限非零源使势沿轻模实方向严格单调，没有有限内部驻点。** 正的局部曲率因此不能作为完整真空稳定的证据。相位仍然无势，下一步优先研究实方向的明确稳定来源。

- [中文结果报告](reports/results_cn.md)：作用量、准确适用范围、推导、基准和后续任务。
- [四联图PNG](figures/light_modulus_diagnostics.png)／[PDF](figures/light_modulus_diagnostics.pdf)：模式分布、松弛重场后的势、非线性动能和有界局部瞬态。
- [计算协议](protocol.md)、[主结果](results/summary.json)、[独立逐行比较](../../reports/light_modulus_output_comparison_20260926.md)。
- 独立专题：[几何](../../reports/light_modulus_geometry_20260926.md)、[有限辅助场驻点](../../reports/light_modulus_stationarity_20260926.md)、[局部有效作用](../../reports/light_modulus_response_20260926.md)。
- [材料和哈希清单](materials_manifest.json)、[归档审计](results/archive_audit.json)。

从仓库根目录复现，依赖Python、NumPy、SciPy、SymPy、mpmath和Matplotlib：

```bash
OPENBLAS_NUM_THREADS=1 python experiments/light_modulus_v01/code/run_light_modulus.py
OPENBLAS_NUM_THREADS=1 python experiments/light_modulus_v01/code/plot_results.py
OPENBLAS_NUM_THREADS=1 python reports/light_modulus_geometry_20260926_checks.py
OPENBLAS_NUM_THREADS=1 python reports/light_modulus_stationarity_20260926_checks.py
OPENBLAS_NUM_THREADS=1 python reports/light_modulus_response_20260926_checks.py
OPENBLAS_NUM_THREADS=1 python reports/light_modulus_output_comparison_20260926.py
python experiments/light_modulus_v01/code/audit_archive.py
```

执行数值脚本会更新相应输出；归档审计只验证当前文件与清单一致。复现环境、时间戳或输出格式造成的哈希差异应先检查，不应直接重写清单来消除诊断。首次归档使用 `audit_archive.py --write-manifest`。

本轮只研究指定作用量的条件结论，没有新增观测拟合或真实宇宙轨迹。图中的约 $2.507\times10^{-17}$ 秒是固定源、平直时空、初始静止的10%局部位移诊断，不是宇宙寿命。86页正式稿与历史实验保留；$1/\ln^2t$ 及真实相互作用的黎曼来源仍未导出。
