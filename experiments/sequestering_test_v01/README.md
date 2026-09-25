# 空间隔离候选与时钟曲率检验

承接 [EFT匹配](../eft_matching_v01/README.md)。本轮导入规范边界五维引力模型的已知有限系数，并有条件地接入本项目四维移位时钟与稳定几何。**小的带电反馈不保证时钟横向稳定。** 同号时钟度规扩展κ=1使原1 eV示例产生负曲率，而较小mG下仍有通过所列必要局部检查的分支。

- [中文结果报告](reports/results_cn.md)
- [参数、来源与附加假设](protocol.md)
- [独立高精度复核](reports/independent_results_cn.md)
- [对称性审计](../../reports/protection_symmetry_20260925.md)
- [五维有限系数来源](../../reports/protection_5d_20260925.md)
- [几何重调与时钟质量推导](../../reports/protection_power_correction_20260925.md)
- [主摘要](results/summary.json)、[条件边界](results/conditional_limits.csv)、[材料清单](materials_manifest.json)

![条件反馈与横向曲率](figures/transfer_and_clock_curvature.png)

代表例：mKK=1TeV、带电质量100GeV、κ=1、ν0=0时，mG=1eV的选定带电力预算约0.422%，但局部横向负曲率对应约20.4秒的快速线性增长尺度。mG=10⁻¹⁸eV时上述局部曲率为正，选定反馈约4.22×10⁻³⁹。后一例通过的是必要条件交集，不是完整UV模型、滚动背景稳定或观测确认。

时钟κ扩展、v⁴稳定来源、同边界接触和完整两环有限匹配仍未确定。原86页正文未改，逆对数平方与真实相互作用的黎曼来源仍待推导。

仓库根目录复算（Python 3、NumPy、SciPy、Matplotlib、SymPy、mpmath）：

```bash
python experiments/sequestering_test_v01/code/run_sequestering.py
python experiments/sequestering_test_v01/code/plot_results.py
python experiments/sequestering_test_v01/code/independent_check.py --derive --compare
python reports/protection_symmetry_20260925_checks.py
python reports/protection_5d_20260925_checks.py
python reports/protection_power_correction_20260925_checks.py
python reports/protection_power_correction_20260925_review.py
```

使用仓库内旧参数和χ坐标，不依赖新观测，不运行旧宇宙ODE。重跑会更新结果时间戳；失败诊断和解释保留，哈希对应归档版本。

不重新计算时，运行 `python experiments/sequestering_test_v01/code/audit_archive.py` 可核对归档哈希、检查记录、文件链接与历史材料未改动。确需归档新一轮输出时，先完成全部科学检查，再用同一脚本的 `--write-manifest` 明确更新材料清单。
