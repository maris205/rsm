#!/usr/bin/env python3
"""Audit saved PM artifacts independently of the production runner.

Audit failure means inconsistent/corrupt reporting, not a failed scientific
resolution criterion. Does not integrate a PM trajectory. Main final states
permit independent geometry, counts, 2-D CIC projection, and legacy full-FFT
force-statistic checks; the other four cases expose only saved aggregates.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy
from scipy.interpolate import CubicSpline, PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
BRIDGE = ROOT.parent / "cosmic_bridge_v01"
PROTOCOL = "f3a0cc5abf21a4374439db756cbf8c4a3103bcabfa2bee6f9bd1921b28e608f2"
LEGACY_SHA = "1b2607ac0b7ba35c5df9ceb7375910ba75b7bed65a1edec79c0583f067a5d982"
MATRIX = {
    "ref_force128": ("reference",64,128,False,(0,0,0)),
    "ref_main": ("reference",64,256,False,(0,0,0)),
    "ref_halfstep": ("reference",64,256,True,(0,0,0)),
    "ref_particles128": ("reference",128,256,False,(0,0,0)),
    "ref_shift": ("reference",64,256,False,(.5,.25,.375)),
    "clock_main": ("clock",64,256,False,(0,0,0)),
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def array_sha(array):
    return hashlib.sha256(np.ascontiguousarray(array).view(np.uint8)).hexdigest()


def minimum(x, center):
    return np.remainder(x-center+.5,1)-.5


class Audit:
    def __init__(self):
        self.checks, self.hashes, self.cases, self.unavailable = [], {}, {}, []

    def track(self, path):
        key=str(path.relative_to(ROOT.parent.parent))
        current=sha(path)
        if key in self.hashes and self.hashes[key]!=current:
            self.check("tracked_input_changed_on_reread",False,path=key)
        else:
            self.hashes[key]=current
        return path

    def check(self, name, condition, **details):
        self.checks.append({"name":name,"passed":bool(condition),**details})

    def close(self, name, value, expected, rtol=2e-11, atol=2e-14):
        x,y=np.broadcast_arrays(np.asarray(value,dtype=float),np.asarray(expected,dtype=float))
        both_nan=np.isnan(x)&np.isnan(y)
        finite=np.isfinite(x)&np.isfinite(y)
        error=np.where(finite,np.abs(x-y),0)
        relative=np.divide(error,np.abs(y),out=np.zeros_like(error),where=finite&(y!=0))
        self.check(name,np.all(finite|both_nan) and np.all(both_nan|(error<=atol+rtol*np.abs(y))),
                   n_values=int(x.size),n_expected_nan=int(both_nan.sum()),
                   max_absolute_error=float(error.max()),max_relative_error_nonzero_expected=float(relative.max()),
                   relative_tolerance=rtol,absolute_tolerance=atol)


def initial_state(nparticle, nmesh, shift_cells, delta, fref, eref):
    axis=(np.arange(nparticle)+.5)/nparticle
    lattice=np.stack(np.meshgrid(axis,axis,axis,indexing="ij"),axis=-1).reshape(-1,3)
    q=minimum(lattice,np.full(3,.5))
    qr=np.linalg.norm(q,axis=1)
    ss=np.clip((qr-.25)/(.4-.25),0,1)
    window=1-10*ss**3+15*ss**4-6*ss**5
    window[qr>=.4]=0
    contrast=delta*window
    yy=(1+contrast)**(-1/3)
    yp=-fref*contrast*yy/(3*(1+contrast))
    shift=np.asarray(shift_cells)/nmesh
    center=np.remainder(np.full(3,.5)+shift,1)
    x=np.remainder(center+q*yy[:,None],1)
    p=.02**2*eref*q*yp[:,None]
    mask=(qr>=.25*.25)&(qr<=.75*.25)
    return x,p,q,qr,mask,center


def geometry(x,p,a,loge,qr,mask,center):
    vectors=minimum(x,center)
    radius=np.linalg.norm(vectors,axis=1)
    corevec,corer=vectors[mask],radius[mask]
    ratio=corer/qr[mask]
    p16,y,p84=np.quantile(ratio,[.16,.5,.84])
    radial=np.einsum("ij,ij->i",p[mask],corevec)/corer
    yp=radial/(a*a*np.exp(loge(np.log(a)))*qr[mask])
    perpendicular=p[mask]-radial[:,None]*corevec/corer[:,None]
    singular=np.linalg.svd(corevec,compute_uv=False)
    axes=singular[-1]/singular[0]
    enclosed=int(np.count_nonzero(radius<=.25*y))
    values={"a":a,"z":1/a-1,"y_median":y,"y_p16":p16,"y_p84":p84,
            "Delta_proxy":y**-3,"nonhomology_fraction":(p84-p16)/(2*y),
            "axis_ratio_min_max":axes,"axis_deviation":1-axes,
            "y_N_median":np.median(yp),"physical_expansion_ratio_median":np.median(ratio+yp),
            "radial_momentum_rms":np.sqrt(np.mean(radial*radial)),
            "transverse_momentum_rms":np.sqrt(np.mean(np.sum(perpendicular**2,axis=1))),
            "radius_proxy_comoving_box":.25*y,
            "actual_enclosed_particle_count_at_proxy_radius":enclosed,
            "actual_enclosed_mean_density_ratio":enclosed/(len(x)*4*np.pi/3*(.25*y)**3)}
    return values,radius


def projected_cic(x,nmesh=128):
    """Independent 2-D assignment, equal to a complete z-sum of 3-D CIC."""
    coordinate=np.remainder(x[:,:2],1)*nmesh
    low=np.floor(coordinate).astype(np.int64)
    fraction=coordinate-low
    flat=np.zeros(nmesh*nmesh)
    for dx in (0,1):
        for dy in (0,1):
            weight=(fraction[:,0] if dx else 1-fraction[:,0])*(fraction[:,1] if dy else 1-fraction[:,1])
            index=((low[:,0]+dx)%nmesh)*nmesh+(low[:,1]+dy)%nmesh
            flat+=np.bincount(index,weights=weight,minlength=nmesh*nmesh)
    return flat.reshape(nmesh,nmesh)*nmesh**2/len(x)


def audit_profiles(audit,name,arrays,x,qr,center):
    radius=np.linalg.norm(minimum(x,center),axis=1)
    edges=np.linspace(0,.45,91)
    count=np.histogram(radius,bins=edges)[0]
    volume=4*np.pi/3*np.diff(edges**3)
    enclosed=np.cumsum(count)
    for field,expected in (("radius_edges_box",edges),("radial_particle_count",count),
                           ("radial_density_over_mean",count/(len(x)*volume)),
                           ("enclosed_particle_count",enclosed),
                           ("enclosed_density_over_mean",enclosed/(len(x)*4*np.pi/3*edges[1:]**3))):
        audit.close(name+"_independent_"+field,arrays[field],expected)
    qedges=np.linspace(0,.4,33)
    qcounts,percentiles=[],[]
    for lo,hi in zip(qedges[:-1],qedges[1:]):
        selected=(qr>=lo)&(qr<hi)
        qcounts.append(int(selected.sum()))
        percentiles.append(np.quantile(radius[selected]/qr[selected],[.16,.5,.84]) if selected.any() else np.zeros(3))
    audit.close(name+"_independent_q_edges",arrays["q_edges_box"],qedges)
    audit.close(name+"_independent_q_shell_counts",arrays["q_shell_count"],qcounts,rtol=0,atol=0)
    audit.close(name+"_independent_q_shell_percentiles",arrays["q_shell_y_percentiles"],percentiles)
    audit.close(name+"_independent_full_depth_2D_CIC_projection",arrays["projection_density_over_mean"],
                projected_cic(x),rtol=2e-11,atol=2e-12)


def verify_reported_checks(audit,name,checks,failures):
    for row in checks:
        value,threshold=row["value"],row["threshold"]
        expected=(False if value is None else bool(value) if threshold is None
                  else bool(np.isfinite(value) and abs(value)<=threshold))
        audit.check(name+"_criterion_verdict_"+row["name"],row["passed"] == expected)
    audit.check(name+"_all_scientific_failures_retained",failures==[r for r in checks if not r["passed"]])


def audit_case(audit,name,args,sphere_summary):
    label,nparticle,nmesh,halfstep,shift=MATRIX[name]
    path=OUT/f"pm_{name}.json"
    if not path.exists():
        failure_path=OUT/f"pm_{name}_execution_failure.json"
        if failure_path.exists():
            failure=json.loads(audit.track(failure_path).read_text())
            audit.check(name+"_execution_failure_retained",failure["status"]=="execution_failed"
                        and failure["protocol_sha256"]==PROTOCOL)
            audit.cases[name]={"status":"execution_failed","error":failure["error"]}
        else:
            audit.unavailable.append(name)
        return
    result=json.loads(audit.track(path).read_text())
    meta=result["metadata"]
    audit.check(name+"_completed_status",result["status"]=="completed")
    audit.check(name+"_frozen_matrix_metadata",(meta["background"],meta["nparticle"],meta["nmesh"],meta["halfstep"])
                ==(label,nparticle,nmesh,halfstep))
    audit.check(name+"_protocol",meta["protocol_sha256"]==PROTOCOL)
    for field,filepath in (("legacy_pm_sha256",BRIDGE/"code/pm.py"),
                           ("pm_refined_sha256",ROOT/"code/pm_refined.py"),
                           ("runner_sha256",ROOT/"code/run_spherical_pm.py"),
                           ("spherical_calibration_sha256",OUT/"sphere_summary.json")):
        audit.check(name+"_hash_"+field,meta[field]==sha(audit.track(filepath)))
    bgpath=BRIDGE/f"results/background_{'eps0' if label=='reference' else 'eps1e-4'}.csv"
    audit.check(name+"_hash_background",meta["background_source_sha256"]==sha(audit.track(bgpath)))
    bg=np.genfromtxt(bgpath,delimiter=",",names=True)
    loge=PchipInterpolator(np.log(bg["a"]),np.log(bg["E"]))
    age=PchipInterpolator(np.log(bg["a"]),bg["t_Gyr"])
    init=sphere_summary["initial"]["reference"]
    x0,p0,q,qr,mask,center=initial_state(nparticle,nmesh,shift,sphere_summary["delta_i"],init["f_reference_i"],init["E_i"])
    # Raw hashes demand bitwise reproduction of the prescribed IC, while all
    # physical readouts below use independent formulas and tolerances.
    audit.check(name+"_initial_positions_bytes",array_sha(x0)==meta["initial"]["positions_sha256"])
    audit.check(name+"_initial_momenta_bytes",array_sha(p0)==meta["initial"]["momenta_sha256"])
    audit.check(name+"_fixed_labels_count",int(mask.sum())==meta["initial"]["core_label_count"])
    audit.close(name+"_center_shift",meta["initial"]["center_box"],center,rtol=0,atol=0)
    audit.close(name+"_particle_mass_Msun_h",meta["initial"]["particle_mass_msun_h"],
                2.77536627e11*.315*10**3/nparticle**3)
    audit.check(name+"_discrete_core_initial_count",meta["initial"]["discrete_core_particle_count"]==int((qr<=.25).sum()))
    history=result["history"]
    aa=np.array([r["a"] for r in history])
    yy=np.array([r["y_median"] for r in history])
    tt=np.array([r["t_Gyr"] for r in history])
    audit.check(name+"_strict_time_order",np.all(np.diff(aa)>0) and np.all(np.diff(tt)>0))
    audit.close(name+"_initial_a",aa[0],.02,rtol=0,atol=1e-15)
    audit.close(name+"_final_a_metadata",aa[-1],meta["a_final"],rtol=0,atol=0)
    audit.close(name+"_saved_background_age",tt,age(np.log(aa)),rtol=2e-12)
    audit.close(name+"_saved_Delta_proxy",[r["Delta_proxy"] for r in history],yy**-3)
    audit.close(name+"_saved_nonhomology_definition",[r["nonhomology_fraction"] for r in history],
                (np.array([r["y_p84"] for r in history])-np.array([r["y_p16"] for r in history]))/(2*yy))
    audit.close(name+"_saved_axis_deviation",[r["axis_deviation"] for r in history],
                1-np.array([r["axis_ratio_min_max"] for r in history]))
    audit.close(name+"_saved_delta_ln_a",[r["delta_ln_a"] for r in history[1:]],np.diff(np.log(aa)))
    nodes,weights=np.polynomial.legendre.leggauss(32)
    mids=(np.log(aa[:-1])+np.log(aa[1:]))/2
    halfwidth=np.diff(np.log(aa))/2
    samples=mids[:,None]+halfwidth[:,None]*nodes
    dtime=halfwidth*np.sum(weights/np.exp(loge(samples)),axis=1)
    tdyn=np.sqrt(2/.315)*aa**1.5*yy**1.5
    audit.close(name+"_independent_interval_cosmic_time",[r["dt_H0"] for r in history[1:]],dtime,rtol=2e-10)
    audit.close(name+"_actual_dt_over_tdyn_saved_readout",[r["actual_dt_over_tdyn"] for r in history[1:]],
                dtime/np.minimum(tdyn[:-1],tdyn[1:]),rtol=2e-10)
    expansion=np.array([r["physical_expansion_ratio_median"] for r in history])
    crossing=np.flatnonzero((expansion[:-1]>0)&(expansion[1:]<=0))
    turnaround=result["turnaround"]
    if len(crossing):
        j=int(crossing[0])
        fraction=expansion[j]/(expansion[j]-expansion[j+1])
        turn_a=np.exp(np.log(aa[j])+fraction*np.log(aa[j+1]/aa[j]))
        audit.check(name+"_physical_turnaround_record_exists",turnaround is not None)
        if turnaround is not None:
            audit.close(name+"_physical_turnaround_first_bracket",turnaround["bracket_a"],aa[j:j+2],rtol=0,atol=0)
            audit.close(name+"_physical_turnaround_interpolation",turnaround["a"],turn_a)
            audit.close(name+"_physical_turnaround_age",turnaround["t_Gyr"],age(np.log(turn_a)))
    else:
        audit.check(name+"_no_false_turnaround",turnaround is None)
    first_geometry,_=geometry(x0,p0,aa[0],loge,qr,mask,center)
    for field,value in first_geometry.items():
        audit.close(name+"_independent_initial_"+field,history[0][field],value,atol=3e-13)
    event=result["event200"]
    if event is None:
        audit.close(name+"_missing_event_stops_at_a_max",aa[-1],.55,atol=2e-14)
        audit.check(name+"_missing_event_no_saved_crossing",np.max(yy**-3)<200)
    else:
        for key,value in event.items():
            if key!="estimator":
                audit.close(name+"_event_matches_final_"+key,value,history[-1][key])
        audit.close(name+"_event_radius_threshold",yy[-1],200**(-1/3),rtol=2e-8)
        audit.check(name+"_event_first_saved_crossing",np.all(yy[:-1]**-3<200)
                    and history[-1].get("terminal_partial_KDK") is True)
    diagnostic=result["diagnostics"]
    mappings={"max_mass_relative_error":"mass_relative_error","max_total_force_normalized":"total_force_normalized",
              "max_total_force_absolute":"total_force_absolute","max_momentum_change_absolute":"momentum_change_absolute",
              "max_drift_cells":"max_drift_cells","max_actual_dt_over_tdyn":"actual_dt_over_tdyn",
              "max_nonhomology_fraction":"nonhomology_fraction","max_axis_deviation":"axis_deviation",
              "maximum_step_halvings":"step_halvings"}
    for saved,field in mappings.items():
        audit.close(name+"_aggregate_"+saved,diagnostic[saved],max(r[field] for r in history))
    audit.check(name+"_step_count",diagnostic["steps"]==len(history)-1)
    common_end=min(aa[-1],result["ode_event200"]["a"])
    common=[r for r in history if r["a"]<=common_end*(1+1e-13) and r["relative_y_error_vs_ode"] is not None]
    audit.close(name+"_aggregate_ODE_radius_error",diagnostic["max_radius_relative_error_vs_ode_common_pre_event"],
                max(abs(r["y_median"]/r["ode_y"]-1) for r in common))
    audit.check(name+"_aggregate_common_sample_count",diagnostic["common_pre_event_samples"]==len(common))
    if event is not None:
        audit.close(name+"_aggregate_ODE_event_error",diagnostic["event_relative_error_vs_ode"],
                    event["a"]/result["ode_event200"]["a"]-1)
    else:
        audit.check(name+"_missing_ODE_event_error_is_null",diagnostic["event_relative_error_vs_ode"] is None)
    ode_csv=np.genfromtxt(audit.track(OUT/f"sphere_{label}.csv"),delimiter=",",names=True)
    ode_interp=CubicSpline(ode_csv["N"],ode_csv["y"])
    before=[r for r in history if r["a"]<=ode_csv["a"][-1]]
    audit.close(name+"_saved_ODE_y_against_archived_sphere",[r["ode_y"] for r in before],
                ode_interp(np.log([r["a"] for r in before])),rtol=1e-6)
    for field in ("a","t_Gyr","y","Delta"):
        audit.close(name+"_ODE_event_against_calibration_"+field,result["ode_event200"][field],
                    sphere_summary["models"][label]["events"]["200"][field])
    criterion_map={"finite":("finite",None),"event200_present":(None,None),
                   "mass":("max_mass_relative_error",1e-12),
                   "total_force":("max_total_force_normalized",1e-11),
                   "actual_drift_cells":("max_drift_cells",.2),
                   "actual_dt_over_tdyn":("max_actual_dt_over_tdyn",.05),
                   "nonhomology":("max_nonhomology_fraction",.1),
                   "axis_deviation":("max_axis_deviation",.1)}
    if name=="ref_main":
        criterion_map.update(PM_vs_ODE_event=("event_relative_error_vs_ode",.02),
                             PM_vs_ODE_radius=("max_radius_relative_error_vs_ode_common_pre_event",.03))
    audit.check(name+"_all_frozen_criteria_present",{r["name"] for r in result["checks"]}==set(criterion_map))
    for row in result["checks"]:
        field,limit=criterion_map[row["name"]]
        expected_value=event is not None if field is None else diagnostic[field]
        audit.check(name+"_frozen_criterion_"+row["name"],row["threshold"]==limit and row["value"]==expected_value)
    verify_reported_checks(audit,name,result["checks"],result["failures"])
    arraypath=audit.track(OUT/f"pm_{name}.npz")
    audit.check(name+"_array_hash",sha(arraypath)==result["array_sha256"])
    arrays=np.load(arraypath,allow_pickle=False)
    for field in ("a","t_Gyr","y_median","Delta_proxy","nonhomology_fraction","axis_deviation","ode_y","relative_y_error_vs_ode"):
        audit.close(name+"_npz_matches_history_"+field,arrays["trajectory_"+field],
                    [np.nan if r.get(field) is None else r[field] for r in history])
    audit.close(name+"_npz_final_a",arrays["final_a"],aa[-1],rtol=0,atol=0)
    audit.close(name+"_npz_center",arrays["center_box"],center,rtol=0,atol=0)
    audit.close(name+"_projection_mean_mass",np.mean(arrays["projection_density_over_mean"]),1,rtol=1e-12)
    audit.check(name+"_projection_nonnegative",np.all(arrays["projection_density_over_mean"]>=0))
    audit.close(name+"_projection_display_axis",arrays["projection_axis_x_box"],(np.arange(128)+.5)/128,rtol=0,atol=0)
    audit.close(name+"_projection_display_axis_y",arrays["projection_axis_y_box"],(np.arange(128)+.5)/128,rtol=0,atol=0)
    audit.close(name+"_radial_count_cumulative",arrays["enclosed_particle_count"],np.cumsum(arrays["radial_particle_count"]),rtol=0,atol=0)
    audit.check(name+"_radial_count_window",arrays["enclosed_particle_count"][-1]<=nparticle**3)
    audit.close(name+"_declared_radial_edges",arrays["radius_edges_box"],np.linspace(0,.45,91),rtol=0,atol=0)
    counts=arrays["radial_particle_count"]
    audit.check(name+"_radial_counts_integer_nonnegative",np.issubdtype(counts.dtype,np.integer) and np.all(counts>=0))
    audit.close(name+"_radial_density_count_consistency",arrays["radial_density_over_mean"],
                counts/(nparticle**3*4*np.pi/3*np.diff(arrays["radius_edges_box"]**3)))
    audit.close(name+"_enclosed_density_count_consistency",arrays["enclosed_density_over_mean"],
                np.cumsum(counts)/(nparticle**3*4*np.pi/3*arrays["radius_edges_box"][1:]**3))
    audit.close(name+"_declared_q_edges",arrays["q_edges_box"],np.linspace(0,.4,33),rtol=0,atol=0)
    expected_qcounts=[np.count_nonzero((qr>=lo)&(qr<hi)) for lo,hi in
                      zip(arrays["q_edges_box"][:-1],arrays["q_edges_box"][1:])]
    audit.close(name+"_fixed_q_shell_label_counts",arrays["q_shell_count"],expected_qcounts,rtol=0,atol=0)
    qpercent=arrays["q_shell_y_percentiles"]
    audit.check(name+"_q_shell_percentiles_finite_ordered",qpercent.shape==(32,3)
                and np.all(np.isfinite(qpercent)) and np.all(qpercent>=0) and np.all(np.diff(qpercent,axis=1)>=0))
    full_state=name in ("ref_main","clock_main")
    force_status="not_available_no_full_state"
    if full_state:
        statepath=audit.track(OUT/f"pm_{name}_final_state.npz")
        audit.check(name+"_state_hash",sha(statepath)==result["state_sha256"])
        state=np.load(statepath,allow_pickle=False)
        x,p=state["positions"],state["momenta"]
        audit.check(name+"_state_float64_shape",x.dtype==np.float64 and p.dtype==np.float64
                    and x.shape==p.shape==(nparticle**3,3))
        audit.check(name+"_state_finite_periodic_positions",np.all(np.isfinite(x)) and np.all(np.isfinite(p))
                    and np.all((x>=0)&(x<1)))
        audit.close(name+"_state_a",state["a"],aa[-1],rtol=0,atol=0)
        audit.close(name+"_state_center",state["center_box"],center,rtol=0,atol=0)
        geo,_=geometry(x,p,aa[-1],loge,qr,mask,center)
        for field,value in geo.items():
            audit.close(name+"_independent_final_"+field,history[-1][field],value,atol=3e-13)
        audit.close(name+"_independent_final_mean_momentum",history[-1]["mean_momentum"],p.mean(axis=0),atol=1e-15)
        audit.close(name+"_independent_final_momentum_change",history[-1]["momentum_change_absolute"],
                    np.linalg.norm(p.mean(axis=0)-p0.mean(axis=0)),atol=1e-15)
        audit_profiles(audit,name,arrays,x,qr,center)
        if not args.skip_fullfft:
            legacy_path=BRIDGE/"code/pm.py"
            audit.check(name+"_legacy_oracle_pinned",sha(legacy_path)==LEGACY_SHA)
            spec=importlib.util.spec_from_file_location("legacy_pm_artifact_oracle",legacy_path)
            legacy=importlib.util.module_from_spec(spec)
            spec.loader.exec_module(legacy)
            force,force_diag=legacy.PMGrid(nmesh).force(x)
            for field in ("raw_mass","density_contrast_mean","mean_force_magnitude","max_force_magnitude"):
                audit.close(name+"_legacy_fullFFT_final_"+field,history[-1][field],force_diag[field],rtol=2e-10,atol=2e-13)
            audit.check(name+"_legacy_fullFFT_finite",np.all(np.isfinite(force)))
            # Nearly zero net-force ratios are cancellation diagnostics: their
            # raw bit patterns need not match between real and complex FFTs.
            audit.check(name+"_legacy_fullFFT_net_force",force_diag["total_force_normalized"]<=1e-11,
                        value=force_diag["total_force_normalized"],threshold=1e-11)
            force_status="one_final_fullFFT_evaluation_no_trajectory_integration"
        else:
            force_status="explicitly_skipped"
    audit.cases[name]={"status":"completed","event200_a":None if event is None else event["a"],
                       "final_a":aa[-1],"final_Delta_proxy":float(yy[-1]**-3),
                       "scientific_failures":[r["name"] for r in result["failures"]],
                       "full_final_state":full_state,"force_audit":force_status,
                       "history_samples":len(history)}


def audit_summary(audit):
    path=OUT/"spherical_pm_summary.json"
    if not path.exists():
        audit.unavailable.append("spherical_pm_summary.json")
        return
    summary=json.loads(audit.track(path).read_text())
    audit.check("summary_protocol",summary["protocol_sha256"]==PROTOCOL)
    for filename,expected in summary["code_sha256"].items():
        audit.check("summary_code_hash_"+filename,sha(audit.track(ROOT/"code"/filename))==expected)
    for filename,expected in summary["input_sha256"].items():
        audit.check("summary_input_hash_"+filename,sha(audit.track(OUT/filename))==expected)
    audit.check("summary_run_coverage_matches_saved_cases",set(summary["runs"])==set(audit.cases))
    actual_failure_count=0
    for name,compact in summary["runs"].items():
        if compact["status"]=="completed":
            path=OUT/f"pm_{name}.json"
            full=json.loads(path.read_text())
            required={"status","event200","turnaround","ode_event200","initial_force_diagnostic",
                      "diagnostics","checks","failures","wall_seconds","array_sha256","state_sha256","json_sha256"}
            audit.check(name+"_summary_required_fields",set(compact)==required)
            audit.check(name+"_summary_json_hash",sha(path)==compact["json_sha256"])
            actual_failure_count+=len(full["failures"])
            for field,value in compact.items():
                if field!="json_sha256":
                    audit.check(name+"_summary_matches_"+field,value==full[field])
        else:
            failure_path=OUT/f"pm_{name}_execution_failure.json"
            failure=json.loads(audit.track(failure_path).read_text())
            audit.check(name+"_summary_execution_failure_matches_archive",compact==failure)
            actual_failure_count+=1
    verify_reported_checks(audit,"summary",summary["checks"],summary["failures"])
    failure_count=actual_failure_count
    audit.check("summary_scientific_case_failure_count",failure_count==summary["case_failure_count"])
    reference=summary["runs"].get("ref_main",{})
    expected_comparisons={}
    if reference.get("status")=="completed" and reference["event200"] is not None:
        base=reference["event200"]
        for other in ("ref_halfstep","ref_particles128","ref_shift"):
            row=summary["runs"].get(other,{})
            if row.get("status")=="completed" and row["event200"] is not None:
                expected_comparisons[other+"_event_relative_to_main"]=row["event200"]["a"]/base["a"]-1
        clock=summary["runs"].get("clock_main",{})
        if clock.get("status")=="completed" and clock["event200"] is not None:
            expected_comparisons["single_clock_pair_fractional_a200_shift"]=clock["event200"]["a"]/base["a"]-1
            expected_comparisons["single_clock_pair_delta_t_Myr"]=1000*(clock["event200"]["t_Gyr"]-base["t_Gyr"])
    audit.check("summary_comparison_availability",set(expected_comparisons)==set(summary["comparisons"]))
    for name,value in expected_comparisons.items():
        audit.close("summary_independent_"+name,summary["comparisons"][name],value)
    expected_global={"all_six_cases_completed":(all(summary["runs"].get(name,{}).get("status")=="completed"
                                                   for name in MATRIX),None)}
    for name,limit in (("ref_halfstep",.005),("ref_particles128",.01),("ref_shift",.01)):
        key=name+"_event_relative_to_main"
        if key in expected_comparisons:
            expected_global[name+"_event_vs_main"]=(expected_comparisons[key],limit)
    audit.check("summary_frozen_global_criteria_present",{r["name"] for r in summary["checks"]}==set(expected_global))
    for row in summary["checks"]:
        expected_value,threshold=expected_global[row["name"]]
        audit.check("summary_frozen_global_criterion_"+row["name"],row["threshold"]==threshold and row["value"]==expected_value)
    audit.summary={"scientific_case_failure_count":failure_count,
                   "scientific_global_failures":[r["name"] for r in summary["failures"]],
                   "available_comparisons":expected_comparisons,
                   "unavailable_comparison_scope":"Missing threshold events yield unavailable convergence/clock comparisons, not successful zero differences."}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--available",action="store_true",help="Audit only currently completed artifacts; label result partial.")
    parser.add_argument("--skip-fullfft",action="store_true",help="Explicitly omit final-state force re-evaluation.")
    args=parser.parse_args()
    audit=Audit()
    audit.track(Path(__file__))
    audit.check("frozen_protocol_bytes",sha(audit.track(ROOT/"protocol.md"))==PROTOCOL)
    sphere_summary=json.loads(audit.track(OUT/"sphere_summary.json").read_text())
    for name in MATRIX:
        try:
            audit_case(audit,name,args,sphere_summary)
        except Exception as exc:
            audit.check(name+"_artifact_audit_completed",False,exception=repr(exc))
    try:
        audit_summary(audit)
    except Exception as exc:
        audit.check("summary_audit_completed",False,exception=repr(exc))
    if not args.available:
        audit.check("all_requested_artifacts_available",not audit.unavailable,missing=audit.unavailable)
    changed=[]
    for relative,expected in audit.hashes.items():
        if sha(ROOT.parent.parent/relative)!=expected:
            changed.append(relative)
    audit.check("audited_inputs_unchanged",not changed,changed=changed)
    result={"created_utc":datetime.now(timezone.utc).isoformat(),"protocol_sha256":PROTOCOL,
            "environment":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__},
            "execution":{"available_subset_mode":args.available,"skip_fullfft":args.skip_fullfft,
                         "fullfft_workers":1,"particle_trajectory_reintegrated":False},
            "complete":not audit.unavailable and len(audit.cases)==6,
            "passed":all(r["passed"] for r in audit.checks),"n_checks":len(audit.checks),
            "n_passed":sum(r["passed"] for r in audit.checks),"checks":audit.checks,
            "cases":audit.cases,"pending_artifacts":audit.unavailable,"input_sha256":audit.hashes,
            "summary_audit":getattr(audit,"summary",None),
            "scope":"Saved-artifact consistency and two full final-state reconstructions; scientific criterion failures are retained, not audit failures.",
            "limitations":["No particle trajectory was reintegrated; historical drifts/forces use saved diagnostics only.",
                           "Only ref_main and clock_main archive full final states; other cases permit aggregate checks only.",
                           "Final fullFFT uses the pinned legacy CIC operator and an independent complex FFT path; it is not an independent continuum gravitational solver.",
                           "Radial profile window is 0<=r<=.45 box, initial-label shells 0<=q<.4; projection sums the entire periodic z depth.",
                           "Saved projection_axis coordinates (j+.5)/128 are display pixel centres; physical CIC nodes are j/128. No physical phase inference uses display axes.",
                           "Equal second-moment axes do not exclude higher-order cubic anisotropy; regular cubic symmetry can preserve axis ratio one despite nonhomologous radial collapse.",
                           "Implementation/statistic audit success does not imply resolved halos, converged collapse, virialization, stars, or a resolved clock response."]}
    target=OUT/"pm_artifact_validation.json"
    if target.exists():
        old=json.loads(target.read_text())
        history=old.get("prior_failed_attempts",[])
        if not old.get("passed",False):
            history.append({k:v for k,v in old.items() if k!="prior_failed_attempts"})
        if history:
            result["prior_failed_attempts"]=history
    target.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    lines=["# PM 保存产物独立审计", "", f"审计检查：{result['n_passed']}/{result['n_checks']} 通过；完整性：{'全六例' if result['complete'] else '当前可用子集'}。", "",
           "这里的通过表示保存产物、独立重建和失败报告一致。科学精度门限未通过不被改写为通过，也不因此自动视为审计错误。", "",
           "| 案例 | 终点 a | Δ200 事件 a | 科学门限失败 | 完整终态 |", "| --- | ---: | ---: | --- | --- |"]
    for name,row in audit.cases.items():
        if row["status"]=="completed":
            event="缺失" if row["event200_a"] is None else f"{row['event200_a']:.10g}"
            lines.append(f"| {name} | {row['final_a']:.10g} | {event} | {', '.join(row['scientific_failures']) or '无'} | {'有' if row['full_final_state'] else '无'} |")
        else:
            lines.append(f"| {name} | — | — | execution_failed：{row['error']} | 无 |")
    lines += ["", "独立重建包括初始粒子标签、共同动量、核心半径分位数、奇异值轴比、实际包围计数、固定窗口径向剖面及完整深度二维 CIC 投影。时间顺序、保存步长的宇宙时间积分、聚合极值、阈值事件与所有失败列表逐项核对。", ""]
    lines += ["- "+item for item in result["limitations"]]
    if audit.unavailable:
        lines += ["", "尚未可用："+", ".join(audit.unavailable)]
    failed=[r for r in audit.checks if not r["passed"]]
    if failed:
        lines += ["", "审计失败："]+["- "+str(row) for row in failed]
    lines += ["", "全部数值、源哈希和历史审计失败保存在 `results/pm_artifact_validation.json`。", "",
              "复现：`python code/validate_pm_artifacts.py`。运行中的可用子集可用 `--available`；`--skip-fullfft` 会明确记录未做终态力重算。", ""]
    (ROOT/"reports/pm_artifact_review.md").write_text("\n".join(lines))
    print(json.dumps({"passed":result["passed"],"complete":result["complete"],"checks":result["n_checks"],
                      "cases":list(audit.cases),"pending":audit.unavailable,
                      "failed":[r["name"] for r in failed]}))
    if not result["passed"]:
        raise SystemExit(1)


if __name__=="__main__":
    main()
