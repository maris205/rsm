#!/usr/bin/env python3
"""Independent saved-artifact audit for the two matched-scale reference cases.

Does not integrate particle trajectories. Initial forces reuse the saved,
hashed continuum-diagnostic force vectors; per completed case, one final force
evaluation uses the archived complex-FFT kernel. Audit evaluations are counted
separately from the continuum diagnostic and production.
Previously audited pure geometry/CSV-audit helpers are reused read-only.
"""
from __future__ import annotations

import argparse
import ast
import importlib.util
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np
import scipy
from scipy.interpolate import CubicSpline, PchipInterpolator

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results"
OLD=ROOT.parent/"halo_cooling_v01"
BRIDGE=ROOT.parent/"cosmic_bridge_v01"
PROTOCOL="4409835c2fd93047d5df7bcb9675b0771e72fc0bfdfb0dd69e8a89c54cae9bf0"
HELPER_SHA="3734d34eb02f7c2798dbbbdcd1f62f8d911120d7861c92fc1a8fe96b40717cdb"
LEGACY_SHA="1b2607ac0b7ba35c5df9ceb7375910ba75b7bed65a1edec79c0583f067a5d982"
CASES={"ref_matched64":64,"ref_matched128":128}


def module_at(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


common=module_at("archived_independent_artifact_helpers",OLD/"code/validate_pm_artifacts.py")
if common.sha(OLD/"code/validate_pm_artifacts.py")!=HELPER_SHA:
    raise RuntimeError("Independent audit helper identity changed")
legacy=module_at("matched_audit_legacy_complexFFT",BRIDGE/"code/pm.py")
if common.sha(BRIDGE/"code/pm.py")!=LEGACY_SHA:
    raise RuntimeError("Legacy fullFFT source identity changed")


def read_json(audit,path):
    return json.loads(audit.track(path).read_text())


def runner_constants(audit):
    tree=ast.parse(audit.track(ROOT/"code/run_matched_pm.py").read_text())
    result={}
    for node in tree.body:
        if isinstance(node,ast.Assign):
            for target in node.targets:
                if isinstance(target,ast.Name) and target.id in ("SOURCE_SHA256","MATRIX","FROZEN_PROTOCOL"):
                    result[target.id]=ast.literal_eval(node.value)
    audit.check("runner_literal_protocol",result["FROZEN_PROTOCOL"]==PROTOCOL)
    expected=(("ref_matched64","reference",64,64,True,(0,0,0)),
              ("ref_matched128","reference",128,128,True,(0,0,0)))
    audit.check("runner_exact_two_reference_halfstep_configurations",result["MATRIX"]==expected)
    for relative,expected_sha in result["SOURCE_SHA256"].items():
        audit.check("archived_input_hash_"+relative,common.sha(audit.track(ROOT.parent/relative))==expected_sha)
    return result


def shape_force(x,force,qr,mask,center,q4_initial=None):
    vectors=common.minimum(x[mask],center)
    radius=np.linalg.norm(vectors,axis=1)
    unit=vectors/radius[:,None]
    g=force[mask]
    radial=np.einsum("ij,ij->i",g,unit)
    transverse=g-radial[:,None]*unit
    radial_norm=np.linalg.norm(radial)
    tangent_norm=np.linalg.norm(transverse)
    q4=np.mean(np.sum(unit**4,axis=1))-.6
    return {"angular_cubic_Q4":float(q4),
            "angular_cubic_Q4_difference_from_initial":0. if q4_initial is None else float(q4-q4_initial),
            "core_radial_force_L2":float(radial_norm),"core_transverse_force_L2":float(tangent_norm),
            "core_transverse_to_radial_force_L2_ratio":None if radial_norm==0 else float(tangent_norm/radial_norm)}


def check_state_archive(audit,name,result,quantity,count,center,final_a):
    meta=result["state_archives"][quantity]
    path=audit.track(OUT/meta["file"])
    audit.check(name+"_"+quantity+"_file_identity",path.name==f"pm_{name}_final_{quantity}.npz"
                and common.sha(path)==meta["sha256"] and path.stat().st_size==meta["bytes"])
    audit.check(name+"_"+quantity+"_size_ceiling",path.stat().st_size<100_000_000)
    with np.load(path,allow_pickle=False) as archive:
        array=archive[quantity]
        audit.close(name+"_"+quantity+"_archive_a",archive["a"],final_a,rtol=0,atol=0)
        audit.close(name+"_"+quantity+"_archive_center",archive["center_box"],center,rtol=0,atol=0)
    audit.check(name+"_"+quantity+"_float64_shape",array.dtype==np.float64 and array.shape==(count,3)
                and meta["dtype"]=="float64" and meta["shape"]==[count,3] and np.all(np.isfinite(array)))
    if quantity=="positions":
        audit.check(name+"_positions_periodic_range",np.all((array>=0)&(array<1)))
    return array


def archived_initial_force(audit,name,n,x0,p0):
    manifest=read_json(audit,OUT/"shell_reference_manifest.json")
    audit.check(name+"_initial_force_manifest_protocol",manifest["protocol_sha256"]==PROTOCOL)
    audit.check(name+"_initial_force_manifest_code",manifest["code_sha256"]==common.sha(audit.track(ROOT/"code/shell_reference.py")))
    for relative,digest in manifest["input_sha256"].items():
        audit.check(name+"_initial_force_source_"+relative,common.sha(audit.track(ROOT.parent.parent/relative))==digest)
    metadata_path=OUT/f"shell_reference_initial_force{n}.json"
    force_path=OUT/f"shell_reference_initial_force{n}.npy"
    for path in (metadata_path,force_path):
        relative=str(path.relative_to(ROOT.parent.parent))
        audit.check(name+"_initial_force_manifest_artifact_"+path.name,
                    manifest["artifact_sha256"][relative]==common.sha(audit.track(path)))
    metadata=json.loads(metadata_path.read_text())
    force=np.load(force_path,allow_pickle=False)
    audit.check(name+"_initial_force_source_configuration",metadata["nparticle_per_side"]==metadata["nmesh_per_side"]==n
                and metadata["force_calls"]==1 and metadata["label"]==name)
    audit.check(name+"_initial_force_matches_particle_order",common.array_sha(x0)==metadata["initial_positions_array_sha256"]
                and common.array_sha(p0)==metadata["initial_momenta_array_sha256"])
    audit.check(name+"_initial_force_bytes_shape",force.dtype==np.float64 and force.shape==(n**3,3)
                and np.all(np.isfinite(force)) and common.array_sha(force)==metadata["force_array_sha256"]
                and common.sha(force_path)==metadata["force_file_sha256"])
    # Force moments below are reconstructed from the archived vector. Density
    # diagnostics remain a cross-file record comparison, not a new deposition.
    diag=dict(metadata["operator_diagnostics"])
    norm=np.linalg.norm(force,axis=1)
    diag.update(mean_force_magnitude=float(norm.mean()),max_force_magnitude=float(norm.max()),
                total_force_normalized=float(np.linalg.norm(force.mean(axis=0))/norm.mean()))
    return force,diag


def audit_case(audit,name,n,source_constants,calibration,force_counts):
    path=OUT/f"pm_{name}.json"
    if not path.exists():
        failure_path=OUT/f"pm_{name}_execution_failure.json"
        if failure_path.exists():
            failure=read_json(audit,failure_path)
            audit.check(name+"_execution_failure_record",failure["status"]=="execution_failed"
                        and failure["protocol_sha256"]==PROTOCOL)
            audit.cases[name]={"status":"execution_failed","error":failure["error"]}
        else:
            audit.unavailable.append(name)
        return
    result=read_json(audit,path)
    meta=result["metadata"]
    audit.check(name+"_case_status",result["status"] in ("completed","stopped_resource_limit"))
    audit.check(name+"_reference_matched_halfstep",meta["background"]=="reference"
                and meta["nparticle"]==meta["nmesh"]==n and meta["halfstep"] is True)
    audit.check(name+"_frozen_protocol",meta["protocol_sha256"]==PROTOCOL)
    audit.check(name+"_source_manifest",meta["source_sha256"]==source_constants["SOURCE_SHA256"])
    for field,path in (("runner_sha256",ROOT/"code/run_matched_pm.py"),
                       ("legacy_pm_sha256",BRIDGE/"code/pm.py"),
                       ("pm_refined_sha256",OLD/"code/pm_refined.py"),
                       ("spherical_calibration_sha256",OLD/"results/sphere_summary.json"),
                       ("background_source_sha256",BRIDGE/"results/background_eps0.csv")):
        audit.check(name+"_"+field,meta[field]==common.sha(audit.track(path)))
    bg=np.genfromtxt(BRIDGE/"results/background_eps0.csv",delimiter=",",names=True)
    loge=PchipInterpolator(np.log(bg["a"]),np.log(bg["E"]))
    age=PchipInterpolator(np.log(bg["a"]),bg["t_Gyr"])
    init=calibration["initial"]["reference"]
    x0,p0,q,qr,mask,center=common.initial_state(n,n,(0,0,0),calibration["delta_i"],init["f_reference_i"],init["E_i"])
    audit.check(name+"_initial_position_bytes",common.array_sha(x0)==meta["initial"]["positions_sha256"])
    audit.check(name+"_initial_momentum_bytes",common.array_sha(p0)==meta["initial"]["momenta_sha256"])
    audit.close(name+"_unretuned_delta",meta["initial"]["delta_nl_core"],.07295182597218097,rtol=0,atol=0)
    audit.check(name+"_fixed_core_label_count",int(mask.sum())==meta["initial"]["core_label_count"])
    history=result["history"]
    aa=np.array([r["a"] for r in history]); yy=np.array([r["y_median"] for r in history])
    audit.check(name+"_ordered_time",np.all(np.diff(aa)>0) and np.all(np.diff([r["t_Gyr"] for r in history])>0))
    audit.close(name+"_a_initial",aa[0],.02,rtol=0,atol=1e-15)
    audit.close(name+"_a_final",meta["a_final"],aa[-1],rtol=0,atol=0)
    audit.close(name+"_background_age",[r["t_Gyr"] for r in history],age(np.log(aa)))
    audit.close(name+"_proxy_density_definition",[r["Delta_proxy"] for r in history],yy**-3)
    audit.close(name+"_scatter_definition",[r["nonhomology_fraction"] for r in history],
                (np.array([r["y_p84"] for r in history])-np.array([r["y_p16"] for r in history]))/(2*yy))
    audit.close(name+"_axis_definition",[r["axis_deviation"] for r in history],1-np.array([r["axis_ratio_min_max"] for r in history]))
    dln=np.diff(np.log(aa))
    if dln.size:
        audit.close(name+"_saved_dln_a",[r["delta_ln_a"] for r in history[1:]],dln)
        audit.check(name+"_halfstep_ln_bound",np.max(dln)<=np.log(.55/.02)/512+2e-14)
        audit.check(name+"_halfstep_predicted_drift_bound",all(r["predicted_max_drift_cells"]<=.075*(1+2e-12) for r in history[1:]))
        audit.check(name+"_halfstep_predicted_time_bound",all(r["predicted_dt_over_tdyn"]<=.015*(1+2e-12) for r in history[1:]))
        nodes,weights=np.polynomial.legendre.leggauss(32)
        mid=(np.log(aa[1:])+np.log(aa[:-1]))/2
        tau=dln/2*np.sum(weights/np.exp(loge(mid[:,None]+dln[:,None]/2*nodes)),axis=1)
        tdyn=np.sqrt(2/.315)*aa**1.5*yy**1.5
        audit.close(name+"_independent_cosmic_interval",[r["dt_H0"] for r in history[1:]],tau,rtol=2e-10)
        audit.close(name+"_actual_interval_over_dynamical_time",[r["actual_dt_over_tdyn"] for r in history[1:]],
                    tau/np.minimum(tdyn[:-1],tdyn[1:]),rtol=2e-10)
    else:
        audit.check(name+"_zero_interval_resource_stop",result["status"]=="stopped_resource_limit"
                    and result["diagnostics"]["steps"]==0,
                    scope="No accepted intervals exist; no historical interval calculation claimed.")
    expansion=np.array([r["physical_expansion_ratio_median"] for r in history])
    changes=np.flatnonzero((expansion[:-1]>0)&(expansion[1:]<=0))
    turn=result["turnaround"]
    if len(changes):
        j=int(changes[0]); fraction=expansion[j]/(expansion[j]-expansion[j+1])
        ta=np.exp(np.log(aa[j])+fraction*np.log(aa[j+1]/aa[j]))
        audit.close(name+"_turnaround_first_bracket",turn["bracket_a"],aa[j:j+2],rtol=0,atol=0)
        audit.close(name+"_turnaround_interpolation",turn["a"],ta)
        audit.close(name+"_turnaround_age",turn["t_Gyr"],age(np.log(ta)))
    else:
        audit.check(name+"_no_false_turnaround",turn is None)
    event=result["event200"]; diag=result["diagnostics"]
    reason=diag["resource_stop_reason"]
    audit.check(name+"_status_resource_reason_consistent",(result["status"]=="stopped_resource_limit")==(reason is not None))
    audit.check(name+"_accepted_step_ceiling",0<=diag["steps"]<=10000)
    if result["status"]=="stopped_resource_limit":
        audit.check(name+"_resource_stop_has_no_event",event is None)
        audit.check(name+"_resource_stop_reason_known",reason in ("total_wall_limit","per_case_step_limit"))
        if reason=="per_case_step_limit": audit.check(name+"_step_limit_exact_count",diag["steps"]==10000)
    if event is not None:
        for field,value in event.items():
            if field!="estimator": audit.close(name+"_event_final_"+field,value,history[-1][field])
        audit.close(name+"_event_target_radius",yy[-1],200**(-1/3),rtol=2e-8)
        audit.check(name+"_first_saved_event",np.all(yy[:-1]**-3<200) and history[-1].get("terminal_partial_KDK") is True)
    else:
        audit.check(name+"_missing_event_not_crossed",np.max(yy**-3)<200)
        if result["status"]=="completed":
            audit.close(name+"_missing_event_stops_at_a_max",aa[-1],.55,rtol=0,atol=2e-14)
        else:
            audit.check(name+"_resource_stop_explicit",diag["resource_stop_reason"] in ("total_wall_limit","per_case_step_limit"))
    for saved,field in {"max_mass_relative_error":"mass_relative_error","max_total_force_normalized":"total_force_normalized",
        "max_total_force_absolute":"total_force_absolute","max_momentum_change_absolute":"momentum_change_absolute",
        "max_drift_cells":"max_drift_cells","max_actual_dt_over_tdyn":"actual_dt_over_tdyn",
        "max_nonhomology_fraction":"nonhomology_fraction","max_axis_deviation":"axis_deviation",
        "maximum_step_halvings":"step_halvings"}.items():
        audit.close(name+"_aggregate_"+saved,diag[saved],max(row[field] for row in history))
    audit.check(name+"_step_count",diag["steps"]==len(history)-1)
    end=min(aa[-1],result["ode_event200"]["a"])
    common_rows=[row for row in history if row["a"]<=end*(1+1e-13) and row["relative_y_error_vs_ode"] is not None]
    audit.close(name+"_radius_error_aggregate",diag["max_radius_relative_error_vs_ode_common_pre_event"],
                max(abs(row["y_median"]/row["ode_y"]-1) for row in common_rows))
    audit.check(name+"_common_pre_event_samples",diag["common_pre_event_samples"]==len(common_rows))
    audit.close(name+"_common_pre_event_end",diag["common_pre_event_end_a"],end)
    reference=np.genfromtxt(audit.track(OLD/"results/sphere_reference.csv"),delimiter=",",names=True)
    interp=CubicSpline(reference["N"],reference["y"])
    subset=[row for row in history if row["a"]<=reference["a"][-1]]
    audit.close(name+"_archived_ODE_radius",[row["ode_y"] for row in subset],interp(np.log([row["a"] for row in subset])),rtol=1e-6)
    for field in ("a","y","Delta","t_Gyr"):
        audit.close(name+"_unretuned_ODE_event_"+field,result["ode_event200"][field],calibration["models"]["reference"]["events"]["200"][field])
    event_error=None if event is None else event["a"]/result["ode_event200"]["a"]-1
    audit.check(name+"_event_error_aggregate",diag["event_relative_error_vs_ode"]==event_error)
    criteria={"finite":(diag["finite"],None),"event200_present":(event is not None,None),
              "mass":(diag["max_mass_relative_error"],1e-12),"total_force":(diag["max_total_force_normalized"],1e-11),
              "actual_drift_cells":(diag["max_drift_cells"],.2),"actual_dt_over_tdyn":(diag["max_actual_dt_over_tdyn"],.05),
              "nonhomology":(diag["max_nonhomology_fraction"],.1),"axis_deviation":(diag["max_axis_deviation"],.1),
              "resource_limits_not_reached":(diag["resource_stop_reason"] is None,None),
              "PM_vs_ODE_event":(event_error,.02),"PM_vs_ODE_radius":(diag["max_radius_relative_error_vs_ode_common_pre_event"],.03)}
    audit.check(name+"_frozen_criteria_complete",{row["name"] for row in result["checks"]}==set(criteria))
    for row in result["checks"]:
        value,threshold=criteria[row["name"]]
        audit.check(name+"_frozen_criterion_"+row["name"],row["value"]==value and row["threshold"]==threshold)
    common.verify_reported_checks(audit,name,result["checks"],result["failures"])
    x=check_state_archive(audit,name,result,"positions",n**3,center,aa[-1])
    p=check_state_archive(audit,name,result,"momenta",n**3,center,aa[-1])
    for stage,xx,pp,scale,row in (("initial",x0,p0,.02,history[0]),("final",x,p,aa[-1],history[-1])):
        geometric,_=common.geometry(xx,pp,scale,loge,qr,mask,center)
        for field,value in geometric.items():
            audit.close(name+"_independent_"+stage+"_"+field,row[field],value,atol=3e-13)
    audit.close(name+"_final_momentum_mean",history[-1]["mean_momentum"],np.mean(p,axis=0),atol=1e-15)
    audit.close(name+"_final_momentum_change",history[-1]["momentum_change_absolute"],np.linalg.norm(p.mean(axis=0)-p0.mean(axis=0)),atol=1e-15)
    arrays_path=audit.track(OUT/f"pm_{name}.npz")
    audit.check(name+"_profile_archive_hash",common.sha(arrays_path)==result["array_sha256"])
    with np.load(arrays_path,allow_pickle=False) as arrays:
        for field in ("a","t_Gyr","y_median","Delta_proxy","nonhomology_fraction","axis_deviation","ode_y","relative_y_error_vs_ode"):
            audit.close(name+"_npz_trajectory_"+field,arrays["trajectory_"+field],[np.nan if row[field] is None else row[field] for row in history])
        audit.close(name+"_npz_final_a",arrays["final_a"],aa[-1],rtol=0,atol=0)
        audit.close(name+"_npz_center",arrays["center_box"],center,rtol=0,atol=0)
        for axis in ("x","y"):
            audit.close(name+"_physical_CIC_node_axis_"+axis,arrays[f"projection_axis_{axis}_box"],np.arange(128)/128,rtol=0,atol=0)
        common.audit_profiles(audit,name,arrays,x,qr,center)
        audit.close(name+"_projection_total_mass",np.mean(arrays["projection_density_over_mean"]),1,rtol=1e-12)
    # Existing independent-diagnostic initial vector; one new final fullFFT.
    mesh=legacy.PMGrid(n)
    q4_initial=None
    force_detail={}
    for stage,xx in (("initial",x0),("final",x)):
        if stage=="initial":
            force,force_diag=archived_initial_force(audit,name,n,x0,p0)
            force_counts["reused_initial_arrays"]+=1
            method="archived_initial_force"
        else:
            force,force_diag=mesh.force(xx)
            force_counts["final"]+=1
            method="fullFFT_final"
        shape=shape_force(xx,force,qr,mask,center,q4_initial)
        if stage=="initial": q4_initial=shape["angular_cubic_Q4"]
        saved=result["force_shape_diagnostics"][stage]
        for field,value in shape.items():
            if value is None: audit.check(name+"_force_shape_"+stage+"_"+field,saved[field] is None)
            else: audit.close(name+"_force_shape_"+stage+"_"+field,saved[field],value,rtol=2e-10,atol=3e-13)
        row=history[0] if stage=="initial" else history[-1]
        for field in ("raw_mass","density_contrast_mean","mean_force_magnitude","max_force_magnitude"):
            audit.close(name+"_"+method+"_"+field,row[field],force_diag[field],rtol=2e-10,atol=2e-13)
        audit.check(name+"_"+method+"_net_force",force_diag["total_force_normalized"]<=1e-11,
                    value=force_diag["total_force_normalized"])
        force_detail[stage]=shape
        if stage=="initial":
            vec=common.minimum(x0[mask],center); radius=np.linalg.norm(vec,axis=1)
            gr=np.sum(force[mask]*vec,axis=1)/radius
            exact=-calibration["delta_i"]*vec/3
            expected={"mean_radial_g":np.mean(gr),"median_radial_g":np.median(gr),
                      "fraction_particles_with_inward_g":np.mean(gr<0),
                      "relative_L2_error_vs_uniform_core":np.linalg.norm(force[mask]-exact)/np.linalg.norm(exact)}
            for field,value in expected.items():
                audit.close(name+"_initial_continuous_force_"+field,result["initial_force_diagnostic"][field],value,rtol=2e-10,atol=2e-13)
    audit.cases[name]={"status":result["status"],"event_a":None if event is None else event["a"],
                       "event_error_vs_ODE":event_error,"scientific_failures":[r["name"] for r in result["failures"]],
                       "final_actual_enclosed_density":history[-1]["actual_enclosed_mean_density_ratio"],
                       "independent_force_shape":force_detail,"steps":diag["steps"],
                       "descriptive_actual_drift_le_point1":diag["max_drift_cells"]<=.1,
                       "descriptive_actual_dt_le_point025":diag["max_actual_dt_over_tdyn"]<=.025}


def audit_summary(audit,constants):
    path=OUT/"matched_pm_summary.json"
    if not path.exists():
        audit.unavailable.append(path.name)
        return None
    summary=read_json(audit,path)
    audit.check("summary_frozen_protocol",summary["protocol_sha256"]==PROTOCOL)
    audit.check("summary_source_manifest",summary["archived_source_sha256"]==constants["SOURCE_SHA256"])
    for name,expected in summary["code_sha256"].items():
        audit.check("summary_code_hash_"+name,common.sha(audit.track(ROOT/"code"/name))==expected)
    for name,expected in summary["input_sha256"].items():
        audit.check("summary_input_hash_"+name,common.sha(audit.track(OUT/name))==expected)
    audit.check("summary_two_case_coverage",set(summary["runs"])==set(CASES)==set(audit.cases))
    failure_count=0
    required={"status","event200","turnaround","ode_event200","initial_force_diagnostic","initial_force_seconds",
              "force_shape_diagnostics","diagnostics","checks","failures","wall_seconds","array_sha256","state_archives","json_sha256"}
    for name,compact in summary["runs"].items():
        if compact["status"]=="execution_failed":
            full=read_json(audit,OUT/f"pm_{name}_execution_failure.json")
            audit.check(name+"_summary_failed_execution_preserved",compact==full)
            failure_count+=1
        else:
            path=OUT/f"pm_{name}.json"; full=read_json(audit,path)
            audit.check(name+"_summary_required_fields",set(compact)==required)
            audit.check(name+"_summary_full_json_hash",compact["json_sha256"]==common.sha(path))
            for field,value in compact.items():
                if field!="json_sha256": audit.check(name+"_summary_"+field,value==full[field])
            failure_count+=len(full["failures"])
    audit.check("summary_independent_scientific_failure_count",failure_count==summary["case_failure_count"])
    low,high=(summary["runs"][name] for name in CASES)
    available=all(row["status"]=="completed" and row["event200"] is not None for row in (low,high))
    difference=high["event200"]["a"]/low["event200"]["a"]-1 if available else None
    audit.check("summary_pair_comparison",summary["comparisons"]=={"matched128_event_relative_to_matched64":difference})
    expected={"matched128_event_vs_matched64":(difference,.01),
              "both_cases_completed":(all(row["status"]=="completed" for row in (low,high)),None)}
    audit.check("summary_all_global_checks_present",{row["name"] for row in summary["checks"]}==set(expected))
    for row in summary["checks"]:
        value,limit=expected[row["name"]]
        audit.check("summary_frozen_global_"+row["name"],row["value"]==value and row["threshold"]==limit)
    common.verify_reported_checks(audit,"summary",summary["checks"],summary["failures"])
    return {"scientific_case_failures":failure_count,"scientific_global_failures":[r["name"] for r in summary["failures"]],
            "high_to_low_event_fractional_difference":difference}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--available",action="store_true",help="Audit available completed cases and label pending artifacts.")
    args=parser.parse_args()
    audit=common.Audit(); counts={"initial":0,"final":0,"reused_initial_arrays":0}
    audit.track(Path(__file__)); audit.track(OLD/"code/validate_pm_artifacts.py"); audit.track(BRIDGE/"code/pm.py")
    audit.check("protocol_bytes",common.sha(audit.track(ROOT/"protocol.md"))==PROTOCOL)
    constants=runner_constants(audit)
    validation_path=OUT/"matched_runner_validation.json"
    if validation_path.exists():
        validation=read_json(audit,validation_path)
        audit.check("matching_runner_implementation_validated",validation["passed"] is True
                    and validation["protocol_sha256"]==PROTOCOL
                    and validation["runner_sha256"]==common.sha(ROOT/"code/run_matched_pm.py"))
    else:
        audit.unavailable.append(validation_path.name)
    calibration=read_json(audit,OLD/"results/sphere_summary.json")
    for name,n in CASES.items():
        try: audit_case(audit,name,n,constants,calibration,counts)
        except Exception as exc: audit.check(name+"_audit_completed",False,exception=repr(exc))
    summary=None
    try: summary=audit_summary(audit,constants)
    except Exception as exc: audit.check("summary_audit_completed",False,exception=repr(exc))
    if not args.available: audit.check("all_artifacts_available",not audit.unavailable,pending=audit.unavailable)
    changed=[relative for relative,digest in audit.hashes.items() if common.sha(ROOT.parent.parent/relative)!=digest]
    audit.check("inputs_unchanged_during_audit",not changed,changed=changed)
    result={"created_utc":datetime.now(timezone.utc).isoformat(),"protocol_sha256":PROTOCOL,
            "environment":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__},
            "passed":all(row["passed"] for row in audit.checks),"complete":not audit.unavailable and len(audit.cases)==2,
            "n_checks":len(audit.checks),"n_passed":sum(row["passed"] for row in audit.checks),
            "checks":audit.checks,"cases":audit.cases,"summary_audit":summary,"pending":audit.unavailable,
            "input_sha256":audit.hashes,"independent_fullFFT_calls":counts,
            "scope":"Artifact consistency and two full-state reconstructions; no PM trajectory integration or clock calculation.",
            "limitations":["Scientific gate failures remain scientific results; artifact audit PASS does not override them.",
                           "Historical drift/force diagnostics are checked as saved records; intermediate particle states are not archived.",
                           "Initial force vectors are reused from hashed continuum-diagnostic artifacts; no new initial force call. Their force moments/Q4 are independently reconstructed, while raw deposition diagnostics are cross-file comparisons.",
                           "One final fullFFT per case is an audit evaluation separate from the initial diagnostic and production trajectory.",
                           "Legacy fullFFT shares the fixed CIC operator and checks discrete execution, not continuum gravitational accuracy.",
                           "Q4 supplements second-moment axes but does not exhaust angular anisotropy or establish three-dimensional convergence.",
                           "Two matched scales and epsilon=0 do not establish phase/time convergence, resolved halo profiles, galaxies, or clock response."]}
    output=ROOT/"reports/independent_artifact_audit.json"
    if output.exists():
        old=json.loads(output.read_text()); failed=old.get("prior_failed_attempts",[])
        if not old.get("passed",False): failed.append({k:v for k,v in old.items() if k!="prior_failed_attempts"})
        if failed: result["prior_failed_attempts"]=failed
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    lines=["# 匹配尺度 PM 产物独立审计", "",f"审计检查 {result['n_passed']}/{result['n_checks']} 通过；完整：{result['complete']}。",
           "", "审计通过表示保存结果及失败报告一致，不表示科学门限通过。", "",
           "| 案例 | 事件 a | 科学门限失败 |", "| --- | ---: | --- |"]
    for name,row in audit.cases.items():
        if row["status"]=="execution_failed": lines.append(f"| {name} | 缺失 | execution_failed |")
        else: lines.append(f"| {name} | {row['event_a']} | {', '.join(row['scientific_failures']) or '无'} |")
    lines += ["",f"独立旧 fullFFT 求值次数：初态 {counts['initial']}，终态 {counts['final']}；复用已存初态力向量 {counts['reused_initial_arrays']} 份，不重跑粒子轨迹。", ""]
    if summary: lines += ["汇总核验：`"+json.dumps(summary,ensure_ascii=False)+"`", ""]
    lines += ["- "+item for item in result["limitations"]]
    if audit.unavailable: lines += ["", "待完成："+", ".join(audit.unavailable)]
    failures=[row for row in audit.checks if not row["passed"]]
    if failures: lines += ["", "审计失败："]+["- "+repr(row) for row in failures]
    lines += ["", "完整数值和哈希见 `independent_artifact_audit.json`。复现：`python code/audit_matched_artifacts.py`。", ""]
    (ROOT/"reports/independent_artifact_audit.md").write_text("\n".join(lines))
    print(json.dumps({"passed":result["passed"],"complete":result["complete"],"checks":result["n_checks"],
                      "force_calls":counts,"pending":audit.unavailable,"failed":[row["name"] for row in failures]}))
    if not result["passed"]: raise SystemExit(1)


if __name__=="__main__": main()
