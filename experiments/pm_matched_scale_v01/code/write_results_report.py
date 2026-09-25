#!/usr/bin/env python3
"""Describe the complete frozen matched-scale matrix without rerunning physics."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
OLD = ROOT.parent / "halo_cooling_v01"
PROTOCOL = "4409835c2fd93047d5df7bcb9675b0771e72fc0bfdfb0dd69e8a89c54cae9bf0"
CASES = ("ref_matched64", "ref_matched128")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def number(value, digits=5, percent=False):
    return "缺失" if value is None else f"{value * (100 if percent else 1):.{digits}f}" + ("%" if percent else "")


def main():
    if sha(ROOT / "protocol.md") != PROTOCOL:
        raise ValueError("Frozen protocol changed")
    summary_path = OUT / "matched_pm_summary.json"
    summary = read(summary_path)
    if summary["protocol_sha256"] != PROTOCOL or set(summary["runs"]) != set(CASES):
        raise ValueError("Complete two-case frozen matrix required")
    sources = [ROOT / "protocol.md", summary_path, OLD / "results/spherical_pm_summary.json",
               OLD / "results/sphere_summary.json"]
    runs = {}
    for name in CASES:
        saved = summary["runs"][name]
        suffix = "_execution_failure" if saved["status"] == "execution_failed" else ""
        path = OUT / f"pm_{name}{suffix}.json"
        run = read(path)
        if suffix:
            if run != saved:
                raise ValueError(f"Failure record mismatch: {name}")
        elif sha(path) != saved["json_sha256"] or run["metadata"]["protocol_sha256"] != PROTOCOL:
            raise ValueError(f"Result identity mismatch: {name}")
        sources.append(path)
        runs[name] = run
    old = read(sources[2])
    ode_event = read(sources[3])["models"]["reference"]["events"]["200"]
    endpoint_log_slope = ode_event["y_N"] / ode_event["y"]
    lines = [
        "# 匹配粒子与力网格尺度：完整结果",
        "", "## Material Passport", "",
        "- 类型：本项目计算实验的完整矩阵描述；非观测数据、无显著性检验。",
        "- 方法：固定初态、参考背景与力核的两档匹配尺度；由已存产物生成，不重积分。",
        "- 完整性：所有预设案例、缺失事件及未通过门限均纳入。独立审计另列，不将统计自洽当作连续精度。",
        f"- 协议SHA256：`{PROTOCOL}`。",
        "", "## 两个预定配置", "",
        "下表半径误差为两方法共同事件前实际保存采样的最大值，同调散布为全保存轨迹最大值；非连续时间上的严格误差上界。事件为固定标签的中位缩放密度代理200，不是真实维里化事件。",
        "", "|案例 Np/Nf|状态|事件 a|相对ODE事件误差|最大半径误差|最大同调散布|末态实际包围密度比|接受步|墙钟/min|",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    values = {}
    for name, run in runs.items():
        if run["status"] == "execution_failed":
            lines.append(f"|{name}|执行失败|缺失|缺失|缺失|缺失|缺失|—|—|")
            values[name] = {"status": run["status"], "error": run["error"]}
            continue
        d, e = run["diagnostics"], run["event200"]
        values[name] = {
            "status": run["status"], "a200": None if e is None else e["a"],
            "event_relative_error": d["event_relative_error_vs_ode"],
            "max_radius_relative_error": d["max_radius_relative_error_vs_ode_common_pre_event"],
            "max_nonhomology": d["max_nonhomology_fraction"],
            "final_actual_enclosed_density": run["history"][-1]["actual_enclosed_mean_density_ratio"],
        }
        crossing = next((row for row in run["history"]
                         if row["relative_y_error_vs_ode"] is not None
                         and row["a"] <= run["diagnostics"]["common_pre_event_end_a"]*(1+1e-13)
                         and abs(row["relative_y_error_vs_ode"]) > .03), None)
        values[name]["first_sample_exceeding_radius_gate"] = None if crossing is None else {
            "a": crossing["a"], "t_Gyr": crossing["t_Gyr"],
            "Delta_proxy": crossing["Delta_proxy"], "ODE_Delta": crossing["ode_y"]**-3,
            "relative_y_error": crossing["relative_y_error_vs_ode"],
        }
        v = values[name]
        label = f"{run['metadata']['nparticle']}/{run['metadata']['nmesh']}"
        lines.append(f"|{label}|{run['status']}|{number(v['a200'],8)}|{number(v['event_relative_error'],4,True)}|{number(v['max_radius_relative_error'],4,True)}|{number(v['max_nonhomology'],4,True)}|{number(v['final_actual_enclosed_density'],4)}|{d['steps']}|{run['wall_seconds']/60:.2f}|")
    difference = summary["comparisons"]["matched128_event_relative_to_matched64"]
    lines += ["", f"两档事件相对差（128/128除以64/64再减1）为 **{number(difference,5,True)}**，事前门限1%。完整矩阵墙钟{summary['wall_seconds']/60:.2f}分钟，两个并发案例，每例最多4个FFT线程。逐例门限失败数{summary['case_failure_count']}；矩阵门限失败数{len(summary['failures'])}。这些是固定阈值判断，不是统计置信度。", "", "## 初态与末态形状诊断", "", "Q4及横向力仅作描述，不设置事后通过线。连续均匀方向Q4期望为零，有限粒子初态不必为零。力比的分母为核心径向力整体L2；它与初始相对连续力L2误差的分母不同。", "", "|案例|初始核心力相对L2误差|Q4初态|Q4末态|Q4变化|末态横向/径向力L2|最大轴偏离|", "|---|---:|---:|---:|---:|---:|---:|"]
    for name, run in runs.items():
        if run["status"] == "execution_failed":
            continue
        initial, final = (run["force_shape_diagnostics"][key] for key in ("initial", "final"))
        lines.append(f"|{name}|{number(run['initial_force_diagnostic']['relative_L2_error_vs_uniform_core'],4,True)}|{initial['angular_cubic_Q4']:.7g}|{final['angular_cubic_Q4']:.7g}|{final['angular_cubic_Q4_difference_from_initial']:.7g}|{number(final['core_transverse_to_radial_force_L2_ratio'],5)}|{number(run['diagnostics']['max_axis_deviation'],4,True)}|")
    lines += ["", "本次两例的径向精度与同调散布改善，但角向读数并非全部单调改善。末态|Q4|由约0.0861降至0.0624，横向／径向力L2比却由约6.49%增至19.86%；不能以轴比接近1或单一径向指标宣布三维收敛。", "", "## 所有未通过项目", ""]
    failures = []
    for name, run in runs.items():
        if run["status"] == "execution_failed":
            failures.append(f"- {name}：执行错误 `{run['error']}`。")
        else:
            for item in run["failures"]:
                failures.append(f"- {name} / `{item['name']}`：值 `{item['value']}`，门限 `{item['threshold']}`。")
    for item in summary["failures"]:
        failures.append(f"- 矩阵 / `{item['name']}`：值 `{item['value']}`，门限 `{item['threshold']}`。")
    lines += failures or ["预先声明的逐例和矩阵门限全部通过；该结论仅限本协议，不自动扩展到未执行的对照。"]
    lines += ["", "## 半径误差何时超过既定界限", "",
              "以下是从完整轨迹追加提取的描述性范围读数：首次保存采样超过原有3%半径门限。它不是新设验收标准，没有提前终止或截短生产结果，也不证明采样之间的严格误差上界。", "",
              "|案例|首次超线的采样a|年龄/Gyr|PM代理密度比|同一时刻连续球密度比|半径误差|",
              "|---|---:|---:|---:|---:|---:|"]
    for name, value in values.items():
        item = value.get("first_sample_exceeding_radius_gate")
        if item is None:
            detail = "无完整轨迹" if value["status"] == "execution_failed" else "共同事件前采样未超线"
            lines.append(f"|{name}|{detail}|—|—|—|—|")
        else:
            lines.append(f"|{name}|{item['a']:.7f}|{item['t_Gyr']:.6f}|{item['Delta_proxy']:.5f}|{item['ODE_Delta']:.5f}|{number(item['relative_y_error'],4,True)}|")
    lines += ["", "## 与上一轮的完整背景对照", "", "旧六例原始材料与结论保持不变。只有ref_halfstep与本轮同属三项限制减半；其他旧例采用基准限制。因此下表不用于拟合统一收敛阶，也不将最接近ODE的例子选为唯一结果。clock_main是旧单一时钟配对，不属于本轮reference矩阵。", "", "|旧案例|事件a|事件误差|最大共同事件前半径误差|最大同调散布|", "|---|---:|---:|---:|---:|"]
    for name in ("ref_force128", "ref_main", "ref_halfstep", "ref_particles128", "ref_shift", "clock_main"):
        run = old["runs"][name]
        d, e = run["diagnostics"], run["event200"]
        lines.append(f"|{name}|{number(None if e is None else e['a'],8)}|{number(d['event_relative_error_vs_ode'],4,True)}|{number(d['max_radius_relative_error_vs_ode_common_pre_event'],4,True)}|{number(d['max_nonhomology_fraction'],4,True)}|")
    lines += ["", "## 可支持的判断与下一步", "", "必须同时使用事件、全轨迹半径及同调散布评价改善程度。粒子/力网格同步加密改变了粒子质量与数值赋值平滑两方面；两档测试只能说明这一固定配比序列的表现，不能唯一确定上一轮误差的全部来源。达到代理200时，理想核心半径仅占约2.74/5.47个力格胞，不能因此证明内部密度结构充分分辨。", "", f"事件和半径误差并非同一个指标。保存的连续解在阈值200处给出d ln y/d ln a={endpoint_log_slope:.5f}，塌缩末期半径随时间变化很快；局部小时间偏移可以产生更大的同一时刻半径偏差。这是连续参考的局部灵敏度说明，并未平移任何曲线来改善误差，也没有证明PM失配只由单一时间延迟造成。", "", "本轮仅epsilon=0，无clock配对，无新独立平移与时间步矩阵。未分辨的微小时钟响应、逆对数平方指数的物理选择及黎曼结构如何导出引力，仍是独立问题。后续应依据全部误差读数选择数值方法；气体守恒、非平衡化学、恒星形成和反馈尚未加入。未执行的[气体接口设计](gas_interface_next_design.md)明确了质量分账、膨胀项、热化及冷却检验顺序；不作为当前已验证成果。", "", "连续密度积分与33壳参考见[shell报告](shell_reference_cn.md)；独立存档审计见[审阅及审计记录](independent_review.md)。各报告的通过只对各自列明的范围有效。", "", "## 复现与产物身份", "", "本报告由 `python experiments/pm_matched_scale_v01/code/write_results_report.py` 只读生成。机器可读数值和逐文件身份位于 `results/reported_values.json`。", ""]
    report = ROOT / "reports/results_cn.md"
    report.parent.mkdir(exist_ok=True)
    report.write_text("\n".join(lines))
    payload = {"protocol_sha256": PROTOCOL, "reported_values": values,
               "pair_fractional_difference": difference, "ode_endpoint_dln_y_dln_a": endpoint_log_slope,
               "inputs": {str(p.relative_to(ROOT.parent.parent)): sha(p) for p in sources},
               "generator_sha256": sha(Path(__file__)), "report_sha256": sha(report)}
    (OUT / "reported_values.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"report": str(report), "case_failures": summary["case_failure_count"], "matrix_failures": len(summary["failures"])}))


if __name__ == "__main__":
    main()
