# 匹配粒子与力网格运行器：小规模一致性验证

状态：**PASS，95/95 项通过**。

本验证只使用 8³ 粒子的短步测试及临时文件，不运行 64³／128³ 生产轨迹，也不重复巨大 DFT 或占用生产初力预检。实现一致性通过不代表连续物理或坍缩事件已经收敛。

冻结协议 SHA256：`4409835c2fd93047d5df7bcb9675b0771e72fc0bfdfb0dd69e8a89c54cae9bf0`。

| 检查 | 状态 | 数值或说明 |
| --- | --- | --- |
| archived_source_hash:halo_cooling_v01/code/pm_refined.py | PASS | 见 JSON 详情 |
| archived_source_hash:halo_cooling_v01/code/spherical_collapse.py | PASS | 见 JSON 详情 |
| archived_source_hash:halo_cooling_v01/code/run_spherical_pm.py | PASS | 见 JSON 详情 |
| archived_source_hash:cosmic_bridge_v01/code/pm.py | PASS | 见 JSON 详情 |
| frozen_protocol_matches_runner | PASS | 见 JSON 详情 |
| exact_two_case_matched_matrix | PASS | 见 JSON 详情 |
| archived_PM_import_identity | PASS | 见 JSON 详情 |
| archived_sphere_import_identity | PASS | 见 JSON 详情 |
| old_canonical_KDK_identity | PASS | 见 JSON 详情 |
| runner_dependency:halo_cooling_v01/code/pm_refined.py | PASS | 见 JSON 详情 |
| runner_dependency:halo_cooling_v01/code/spherical_collapse.py | PASS | 见 JSON 详情 |
| runner_dependency:halo_cooling_v01/code/run_spherical_pm.py | PASS | 见 JSON 详情 |
| runner_dependency:halo_cooling_v01/results/sphere_summary.json | PASS | 见 JSON 详情 |
| runner_dependency:halo_cooling_v01/results/refined_pm_validation.json | PASS | 见 JSON 详情 |
| runner_dependency:cosmic_bridge_v01/code/pm.py | PASS | 见 JSON 详情 |
| runner_dependency:cosmic_bridge_v01/results/background_summary.json | PASS | 见 JSON 详情 |
| runner_dependency:cosmic_bridge_v01/results/background_eps0.csv | PASS | 见 JSON 详情 |
| runner_dependency:cosmic_bridge_v01/results/background_eps1e-4.csv | PASS | 见 JSON 详情 |
| frozen_constant:EVENT_THRESHOLD | PASS | 见 JSON 详情 |
| frozen_constant:RADIUS_THRESHOLD | PASS | 见 JSON 详情 |
| frozen_constant:PAIR_EVENT_THRESHOLD | PASS | 见 JSON 详情 |
| frozen_constant:SHAPE_THRESHOLD | PASS | 见 JSON 详情 |
| frozen_constant:MAX_TOTAL_WALL_SECONDS | PASS | 见 JSON 详情 |
| frozen_constant:MAX_STEPS | PASS | 见 JSON 详情 |
| frozen_constant:MAX_ARCHIVE_BYTES | PASS | 见 JSON 详情 |
| both_cases_apply_frozen_case_thresholds | PASS | 见 JSON 详情 |
| no_name_condition_around_ODE_thresholds | PASS | 见 JSON 详情 |
| actual_tdyn_uses_minimum_of_both_endpoints | PASS | 见 JSON 详情 |
| unchanged_calibrated_initial_amplitude | PASS | 0.07295182597218097 |
| independent_initial_mapping:positions | PASS | 5.551115123125783e-17 |
| independent_initial_mapping:momenta | PASS | 7.792703114739563e-19 |
| independent_initial_mapping:q | PASS | 0.0 |
| independent_initial_mapping:q_radius | PASS | 0.0 |
| independent_fixed_core_labels | PASS | 见 JSON 详情 |
| initial_physics_independent_of_force_mesh:positions | PASS | 见 JSON 详情 |
| initial_physics_independent_of_force_mesh:momenta | PASS | 见 JSON 详情 |
| initial_physics_independent_of_force_mesh:q | PASS | 见 JSON 详情 |
| initial_physics_independent_of_force_mesh:q_radius | PASS | 见 JSON 详情 |
| initial_physics_independent_of_force_mesh:core_mask | PASS | 见 JSON 详情 |
| independent_EdS_KDK_positions:a=0.0201 | PASS | 5.551115123125783e-17 |
| independent_EdS_KDK_momenta:a=0.0201 | PASS | 2.1560468178160837e-16 |
| independent_EdS_KDK_positions:a=0.0202 | PASS | 5.551115123125783e-17 |
| independent_EdS_KDK_momenta:a=0.0202 | PASS | 1.338262048315174e-16 |
| independent_EdS_KDK_positions:a=0.0203 | PASS | 5.551115123125783e-17 |
| independent_EdS_KDK_momenta:a=0.0203 | PASS | 6.387725063135929e-16 |
| halfstep_independent_endpoint:ordinary | PASS | 0.0 |
| halfstep_independent_halvings:ordinary | PASS | 见 JSON 详情 |
| halfstep_maximum_delta_ln_a:ordinary | PASS | 0.00647301954037597 |
| halfstep_drift_bound:ordinary | PASS | 0.00017037792694998084 |
| halfstep_dynamical_bound:ordinary | PASS | 0.004764223590933273 |
| halfstep_independent_endpoint:large_momentum | PASS | 0.0 |
| halfstep_independent_halvings:large_momentum | PASS | 见 JSON 详情 |
| halfstep_maximum_delta_ln_a:large_momentum | PASS | 0.00647301954037597 |
| halfstep_drift_bound:large_momentum | PASS | 0.05540487481830872 |
| halfstep_dynamical_bound:large_momentum | PASS | 0.004764223590933273 |
| halfstep_independent_endpoint:short_dynamical_time | PASS | 0.0 |
| halfstep_independent_halvings:short_dynamical_time | PASS | 见 JSON 详情 |
| halfstep_maximum_delta_ln_a:short_dynamical_time | PASS | 0.0004045637212733886 |
| halfstep_drift_bound:short_dynamical_time | PASS | 0.0 |
| halfstep_dynamical_bound:short_dynamical_time | PASS | 0.012797310782003015 |
| halfstep_test_exercises_reduction | PASS | 见 JSON 详情 |
| partial_event_closed_form_scale_factor | PASS | 1.509903313490213e-14 |
| partial_event_prediction_matches_complete_KDK_position | PASS | 0.0 |
| partial_event_reaches_density_proxy_200 | PASS | 1.3322676295501878e-15 |
| partial_event_reports_root_calls | PASS | 见 JSON 详情 |
| unbracketed_event_is_not_fabricated | PASS | 见 JSON 详情 |
| projection_nodal_coordinate:projection_axis_x_box | PASS | 见 JSON 详情 |
| projection_nodal_coordinate:projection_axis_y_box | PASS | 见 JSON 详情 |
| coordinate_correction_preserves_all_density_and_profile_values | PASS | 见 JSON 详情 |
| projected_mass_normalization | PASS | 0.0 |
| missing_event_is_null_and_failed:ref_matched64 | PASS | 见 JSON 详情 |
| missing_event_is_null_and_failed:ref_matched128 | PASS | 见 JSON 详情 |
| missing_cases_explicitly_fail | PASS | 见 JSON 详情 |
| pair_event_threshold:a=0.503 | PASS | 见 JSON 详情 |
| pair_event_threshold:a=0.506 | PASS | 见 JSON 详情 |
| missing_scalar_check_fails | PASS | 见 JSON 详情 |
| separate_position_momentum_archives | PASS | 见 JSON 详情 |
| state_archive_hash:positions | PASS | 见 JSON 详情 |
| state_archive_bytes:positions | PASS | 见 JSON 详情 |
| state_archive_shape_dtype_manifest:positions | PASS | 见 JSON 详情 |
| state_archive_exact_roundtrip:positions | PASS | 见 JSON 详情 |
| state_archive_event_metadata:positions | PASS | 见 JSON 详情 |
| state_archive_hash:momenta | PASS | 见 JSON 详情 |
| state_archive_bytes:momenta | PASS | 见 JSON 详情 |
| state_archive_shape_dtype_manifest:momenta | PASS | 见 JSON 详情 |
| state_archive_exact_roundtrip:momenta | PASS | 见 JSON 详情 |
| state_archive_event_metadata:momenta | PASS | 见 JSON 详情 |
| state_archive_overwrite_rejected_without_changes | PASS | 见 JSON 详情 |
| non_float64_final_state_rejected | PASS | 见 JSON 详情 |
| nonfinite_final_state_rejected | PASS | 见 JSON 详情 |
| missing_events_remain_JSON_null | PASS | 见 JSON 详情 |
| angular_Q4_independent_definition | PASS | 0.0 |
| radial_force_has_zero_transverse_component | PASS | 1.1365894944938756e-16 |
| zero_radial_force_ratio_is_missing_not_NaN_or_zero | PASS | 见 JSON 详情 |
| archived_directories_unchanged | PASS | 见 JSON 详情 |

旧实验目录在导入运行器之前及测试结束后逐文件校验；禁止隐式生成旧目录字节码。分片测试只在临时目录进行，检查顺序、形状、float64、范围、逐片哈希及完整往返。旧实验的物理离散失败不因新运行器一致性通过而改变。
