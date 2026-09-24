"""Summarize saved results without rerunning or selecting physical cases."""
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
CASE_LABELS = {
    'ref_force128': '参考，64³粒子／128³力网格',
    'ref_main': '参考主例，64³／256³',
    'ref_halfstep': '参考，64³／256³，三项步长界限减半',
    'ref_particles128': '参考，128³／256³',
    'ref_shift': '参考，64³／256³，共同平移',
    'clock_main': 'clock主例，64³／256³',
}


TEMP_LABELS = {
    'below_atomic_diagnostic_range': '低于原子诊断范围，T < 10⁴ K',
    'low_temperature_sensitive': '低温敏感区，10⁴ ≤ T < 10^4.5 K',
    'declared_main_temperature_range': '主解释温区，10^4.5 ≤ T ≤ 10⁸ K',
    'above_declared_main_temperature_range': '高于主解释温区，T > 10⁸ K',
}


def read(name):
    return json.loads((ROOT / name).read_text())


def main():
    sphere = read('results/sphere_summary.json')
    cooling = read('results/cooling_summary.json')
    pm = read('results/spherical_pm_summary.json')
    force = read('results/refined_pm_validation.json')
    independent = read('results/physics_independent_validation.json')
    pm_audit = read('results/pm_artifact_validation.json')
    if not pm_audit['complete'] or not pm_audit['passed'] or len(pm['runs']) != 6:
        raise RuntimeError('Full six-case saved-artifact audit required before final report')
    s = sphere['paired_event_shifts']['200']
    lines = [
        '# 受控球形塌缩与条件气体冷却：完整结果', '',
        '本轮完成了有物理单位的球形塌缩参照、六组粒子—网格演化以及原初H/He气体的瞬时冷却接口。连续球形与冷却实现得到独立核验，但三维主配置没有达到预设连续精度。三者的验证范围分别报告；数值实现准确不等于离散粒子已接近连续球，也不等于观测支持。协议在生产前冻结，未重新选择初幅或丢弃失败配置。', '',
        '## 球形方程与共同初态', '',
        f"参考方程一次校准初幅 δ_i={sphere['delta_i']:.12f}，目标是首次密度比200事件在a=0.5。这个预设基准不含观测拟合。独立物理半径解、精度收紧及EdS摆线对照共{sphere['summary']['checks_passed']}/{sphere['summary']['checks_total']}项通过。", '',
        '| 量 | 参考背景 | clock背景 |', '| --- | ---: | ---: |',
        f"| Δ200的尺度因子 | {s['reference_a']:.10f} | {s['clock_a']:.10f} |",
        f"| Δ200的年龄/Gyr | {s['reference_t_Gyr']:.10f} | {s['clock_t_Gyr']:.10f} |", '',
        f"共同初态与实际初动量下，示例clock的事件尺度因子增加{100*s['fractional_delta_a']:.7f}%，年龄增加{s['delta_t_Myr']:.6f} Myr。该高精度ODE差不是已经由PM分辨的物理信号。Δ200仍是内落阶段，不称维里平衡或星系。", '',
        '## 六组PM：所有结果与失败保留', '',
        '下面误差以各背景自己的连续球形ODE为参照。半径采用固定拉格朗日标签median(r/q)，仅为同调代理；密度、形状及包围计数另存在逐例JSON。共同事件前的最大半径误差与事件误差是不同指标。', '',
        '| 配置 | 代理达到200的a | 事件相对误差/% | 事件前最大半径误差/% | 同调散布最大值/% | 事件实际包围密度比 | 失败项数 |',
        '| --- | ---: | ---: | ---: | ---: | ---: | ---: |',
    ]
    for name in CASE_LABELS:
        run = pm['runs'][name]
        label = CASE_LABELS[name] + f' (`{name}`)'
        if run.get('status') != 'completed':
            lines.append(f'| {label} | 执行失败 | — | — | — | — | 保留日志 |')
            continue
        d = run['diagnostics']
        a = '未达到' if run['event200'] is None else f"{run['event200']['a']:.8f}"
        err = '—' if d['event_relative_error_vs_ode'] is None else f"{100*d['event_relative_error_vs_ode']:.5f}"
        actual = "—" if run["event200"] is None else f"{run['event200']['actual_enclosed_mean_density_ratio']:.3f}"
        lines.append(f"| {label} | {a} | {err} | {100*d['max_radius_relative_error_vs_ode_common_pre_event']:.4f} | {100*d['max_nonhomology_fraction']:.4f} | {actual} | {len(run['failures'])} |")
    lines += ['', f"全矩阵协调器墙钟为{pm['wall_seconds']/60:.2f}分钟，最多{pm['concurrent_production_jobs']}个并行作业、每作业4个FFT线程。状态依赖步数及粒子数使各例成本不同；这些是执行成本记录，不是受控速度基准，也不是与其他引力代码的同物理比较。", '', '预先固定的主要门限为：主例事件误差≤2%、共同事件前半径误差≤3%、步长减半的事件差≤0.5%、增加粒子／平移的事件差≤1%、同调与轴形偏离≤10%。未通过的项目不能通过放宽门限改写为成功。只有一个clock PM配对，没有对clock重复完整收敛矩阵。', '',
              f"逐例检查失败共{pm['case_failure_count']}项，跨例汇总失败{len(pm['failures'])}项。该计数用于定位协议门限，不是统计显著性。", '',
              '配对／分辨率比较原值：', '', '```json', json.dumps(pm['comparisons'], indent=2, ensure_ascii=False), '```', '',
              '每例失败详情：', '']
    for name in CASE_LABELS:
        run = pm['runs'][name]
        for failure in run.get('failures', []):
            lines.append(f"- {name}: `{failure['name']}`；值={failure.get('value')}，门限={failure.get('threshold')}。")
    for failure in pm['failures']:
        lines.append(f"- 汇总：`{failure['name']}`；值={failure.get('value')}，门限={failure.get('threshold')}。")
    lines += ['', '初力、独立离散DFT、RFFT与旧FFT对照及净力检查见[完整力验证](../results/refined_pm_validation.json)及[分类说明](refined_pm_validation_cn.md)。细力网格不能增加粒子采样信息；其连续球初力偏差另列，不能由离散算子的一致性推出连续物理精度。首次验证记录里的额外方向诊断及后续分类说明全部保留。规则立方格的二阶矩轴比接近1并不排除高阶角向误差，因此同时报告标签径向散布。', '',
              '## 条件冷却，不是气体演化', '',
              f"历史CIE实现检查{cooling['validation_n_checks']}项通过；生产与官方解析参照分别实现。固定质量和红移的微观读数最大绝对差为{cooling['fixed_Mz_microphysics_max_absolute_difference']}。两背景共有175个唯一条件状态，350行不是350个独立对象。", '',
              '| 温区 | 状态数 | 化学弛豫慢于冷却 | 化学弛豫慢于动力学 |', '| --- | ---: | ---: | ---: |']
    for name, label in TEMP_LABELS.items():
        row = cooling['unique_grid_warning_counts']['by_temperature_regime'][name]
        lines.append(f"| {label} | {row['n_unique_Mz_states']} | {row['n_CIE_slower_than_cooling']} | {row['n_CIE_slower_than_dynamics']} |")
    lines += ['', '参考Δ200事件的条件特征温度约4.10×10^6 K，瞬时冷却时间27.02095 Gyr，动力学时间0.913874 Gyr，前者约为后者29.57倍。clock事件对应27.03027 Gyr和0.914081 Gyr，差异来自事件红移不同。气体密度取fb×200ρm，不是PM直接测得的气体密度；温度也不是已经模拟出的热化结果。', '',
              '局部CIE弛豫时间只是给定温度与密度附近的时间尺度诊断。即使通过这一检查，仍未检验光学深度、外部紫外场、金属、分子、真实热历史或反馈。低温及平衡建立较慢的点保持可见。tcool<tdyn不独自推出恒星形成，tcool>tdyn也不排除后来更稠密／富金属环境中的冷却。', '',
              '## 独立复算、来源与复现', '',
              f"另一套宇宙时间物理半径积分及保存冷却数组的独立复算通过{independent['n_passed']}/{independent['n_checks']}项。它不重跑大PM、不验证真实宇宙模型成立；详见[独立报告](physics_independent_review.md)。", '',
              '来源：Pace等球塌缩方法、Katz等原初气体速率和固定版本Grackle解析工具。完整引文与许可见[论文参考58—62](../../../manuscript_cn/references.md#r58)和[来源清单](../sources/grackle_source_metadata.json)。', '',
              '[固定协议](../protocol.md)、[运行命令](../README.md)、[球形摘要](../results/sphere_summary.json)、[PM摘要](../results/spherical_pm_summary.json)、[冷却摘要](../results/cooling_summary.json)均随仓库保存。新计算没有使用观测数据选择初幅、时钟参数或速率。 本轮没有比较响应幂；背景时钟的影响不能直接当作χ⁻²响应形式优于其他形式的证据。', '',
              '后续顺序是先使三维引力误差得到控制，再用守恒气体方程和非平衡化学演化温度／密度；恒星形成、反馈及观测比较需要在此基础上继续。黎曼结构到真实作用的推导与逆对数响应的物理选择仍属独立未完成的问题。', '']
    lines += [f"力算子实现检查{force['implementation_checks_passed']}/{force['implementation_checks_total']}项通过；额外物理离散诊断{force['physical_diagnostic_checks_passed']}/{force['physical_diagnostic_checks_total']}项通过，总状态仍保留FAIL。连续初态和补偿质量另有59项通过检查，详见[sphere_results_cn.md](sphere_results_cn.md)。", '']
    lines += [f"六例保存产物及两个完整粒子末态的独立审计通过{pm_audit['n_passed']}/{pm_audit['n_checks']}项，详见[PM存档核验](pm_artifact_review.md)。它验证保存结果与报告一致，科学精度失败保持原样。", '']
    path = ROOT/'reports/results_cn.md'
    path.write_text('\n'.join(lines))
    print(path)


if __name__ == '__main__':
    main()
