#!/usr/bin/env python3
"""One-loop matching of a specified minimal U(1); archived backgrounds only."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import platform
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
PARENT = REPO / 'experiments/sugra_shift_v01'
OUT = ROOT / 'results'
PI = np.pi


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def write_json(path, data):
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def read_axes():
    with (PARENT/'results/trajectories.csv').open() as f:
        rows = list(csv.DictReader(f))
    result = {}
    for model in ('global', 'canonical', 'shift'):
        selected = [r for r in rows if r['case'] == model+'_axis']
        result[model] = {k: np.array([float(r[k]) for r in selected])
                         for k in selected[0] if k not in ('case', 'model')}
    return result


def main():
    start = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    inputs = json.loads((PARENT/'results/inputs.json').read_text())
    alpha = inputs['alpha_reference']
    eps, xi = inputs['epsilon'], inputs['xi']
    all_rows, summaries, checks = [], [], []

    def check(name, value, limit=1e-11):
        checks.append(dict(name=name, value=float(value), limit=limit,
                           passed=bool(np.isfinite(value) and value < limit)))

    for model, r in read_axes().items():
        chi = r['chi']; chi0 = chi[-1]
        k = eps*chi**2/2 if model == 'canonical' else np.zeros_like(chi)
        dk = k-k[-1]
        # Relative holomorphic mass logarithm; no absolute large-log subtraction.
        X = xi*(chi0-chi)*(chi0+chi)/(chi**2*(chi0**2+xi))
        lm = np.log1p(X)
        ls = lm+dk
        dr, hr2 = r['d']/r['s'], r['h2']/r['s']**2
        Rminus1 = 2*dr+dr**2-hr2
        lr = np.log1p(Rminus1); dlr = lr-lr[-1]
        qh0 = alpha/(2*PI)*lm
        qr = alpha/(12*PI)*dlr
        qh = qh0+qr
        qp = alpha/(4*PI)*(2*ls+dlr/3)
        hol = qh/(1-qh)
        physical = qp/(1-qp)
        split_delta = qr/((1-qh0)*(1-qh0-qr))
        mass_Y = -(2*ls+dlr/3)/(16*PI**2)
        anomaly_Y = dk/(8*PI**2)
        holo_Y = mass_Y+anomaly_Y
        simplified_Y = -(2*lm+dlr/3)/(16*PI**2)
        physical_boundary_f = -anomaly_Y
        scale = max(float(np.max(abs(simplified_Y))), 1e-30)
        check(model+'_high_to_low_vs_holomorphic', np.max(abs(holo_Y-simplified_Y))/scale)
        check(model+'_constant_physical_uv_recovery',
              np.max(abs(physical-r['delta_alpha']))/max(float(np.max(abs(physical))), 1e-30), 1e-10)
        check(model+'_inverse_coupling_reconstruction',
              np.max(abs(hol-(-4*PI*alpha*holo_Y)/(1+4*PI*alpha*holo_Y)))/max(float(np.max(abs(hol))),1e-30))
        check(model+'_positive_R', 0. if np.all(Rminus1 > -1) else 1., .5)
        check(model+'_positive_scalar_spectrum', 0. if np.all(1+dr > np.sqrt(hr2)) else 1., .5)
        check(model+'_UV_above_pair', float(np.max(np.sqrt(r['s']))/1e13), 1.)
        check(model+'_IR_coupling_positive', 0. if np.all(1-qh > 0) else 1., .5)

        # Same physical theory in a true Kahler-transformed frame, not the other model.
        j = -eps*chi**2/4 if model != 'global' else np.zeros_like(chi)
        dj = j-j[-1]
        lm_K = lm-2*dj
        dk_K = dk+2*dj
        df_K = -dj/(4*PI**2)
        Y_K = df_K+dk_K/(8*PI**2)-(2*(lm_K+dk_K)+dlr/3)/(16*PI**2)
        check(model+'_kahler_physical_mass_invariant', np.max(abs(lm_K+dk_K-ls)))
        check(model+'_kahler_coupling_invariant', np.max(abs(Y_K-simplified_Y))/scale)
        kahler_wrong = -2*lm_K/(16*PI**2)-dlr/(48*PI**2)

        # Holomorphic Q' = exp(u) Q; anomalies must accompany the field redefinition.
        uplus = .07*(chi-inputs['chi_i'])
        uminus = -.03*(chi-inputs['chi_i'])
        du = uplus+uminus-(uplus[-1]+uminus[-1])
        lnZ_new = -2*du
        lm_Q = lm-2*du
        df_Q = -du/(4*PI**2)
        ls_Q = dk+lm_Q-lnZ_new
        anomaly_Q = dk/(8*PI**2)-lnZ_new/(8*PI**2)
        Y_Q = df_Q+anomaly_Q-(2*ls_Q+dlr/3)/(16*PI**2)
        check(model+'_konishi_physical_mass_invariant', np.max(abs(ls_Q-ls)))
        check(model+'_konishi_coupling_invariant', np.max(abs(Y_Q-simplified_Y))/scale)
        Y_Q_wrong = anomaly_Q-(2*ls_Q+dlr/3)/(16*PI**2)
        # Negative controls: omitting the transformed UV term must matter when nonzero.
        if model != 'global':
            check(model+'_kahler_omission_detected',
                  0. if np.max(abs(kahler_wrong-simplified_Y))/scale > .01 else 1., .5)
        check(model+'_konishi_omission_detected',
              0. if np.max(abs(Y_Q_wrong-simplified_Y))/scale > .01 else 1., .5)

        # Shape diagnostics have no fitted amplitude: compare normalized endpoints.
        chi_response = (chi0-chi)*(chi0+chi)/(chi**2*chi0**2)
        t_sec = r['t_Gyr']*(365.25*86400*1e9)
        logt = np.log(t_sec/inputs['tstar_s'])
        t_response = (logt[-1]-logt)*(logt[-1]+logt)/(logt**2*logt[-1]**2)
        shape_chi = chi_response/chi_response[0]
        shape_time = t_response/t_response[0]
        shape_hol = hol/hol[0]
        q_taylor = alpha/(2*PI)*X
        approx = q_taylor/(1-q_taylor)
        shape_gap = np.max(abs(shape_chi-shape_time))
        taylor_error = np.max(abs(approx-hol))/np.max(abs(hol))
        lrprime = (2*(1+dr)*(r['dp']/r['s']-dr*r['sp']/r['s'])
                   -r['h2p']/r['s']**2+2*hr2*r['sp']/r['s'])/(1+Rminus1)
        lmprime = -2*xi/(chi*(chi**2+xi))
        yr = 365.25*86400
        speed_today = r['wx'][-1]*inputs['H_ref_s_inv']*yr
        drift_hol = alpha/(4*PI)*(2*lmprime[-1]+lrprime[-1]/3)*speed_today
        drift_phys = alpha/(4*PI)*(2*r['sp'][-1]/r['s'][-1]+lrprime[-1]/3)*speed_today
        summary = dict(model=model,points=len(chi),holomorphic_uv_zi_ppm=float(hol[0]*1e6),
            physical_uv_zi_ppm=float(physical[0]*1e6),
            mass_piece_delta_inverse_alpha_zi=float(4*PI*mass_Y[0]),
            anomaly_piece_delta_inverse_alpha_zi=float(4*PI*anomaly_Y[0]),
            total_holo_delta_inverse_alpha_zi=float(4*PI*simplified_Y[0]),
            split_delta_alpha_zi=float(split_delta[0]),max_abs_split_delta_alpha=float(np.max(abs(split_delta))),
            holomorphic_uv_today_drift_per_year=float(drift_hol),
            physical_uv_today_drift_per_year=float(drift_phys),
            max_mass_log_taylor_relative_peak_error=float(taylor_error),
            max_normalized_chi_vs_time_shape_gap=float(shape_gap),
            max_normalized_matched_vs_chi_shape_gap=float(np.max(abs(shape_hol-shape_chi))),
            max_abs_X=float(np.max(abs(X))),max_abs_log_R=float(np.max(abs(lr))),
            min_one_minus_q=float(np.min(1-qh)),
            max_kahler_wrong_boundary_relative_Y=float(np.max(abs(kahler_wrong-simplified_Y))/scale),
            max_konishi_wrong_boundary_relative_Y=float(np.max(abs(Y_Q_wrong-simplified_Y))/scale))
        summaries.append(summary)
        columns = dict(N=r['N'],z=r['z'],chi=chi,t_Gyr=r['t_Gyr'],k=k,delta_k=dk,
                       delta_log_Mhol2=lm,delta_log_s=ls,log_R=lr,delta_log_R=dlr,
                       mass_delta_g_inverse2=mass_Y,anomaly_delta_g_inverse2=anomaly_Y,
                       physical_uv_delta_f=physical_boundary_f,matched_delta_g_inverse2=simplified_Y,
                       holomorphic_uv_delta_alpha=hol,physical_uv_delta_alpha=physical,
                       split_delta_alpha=split_delta,taylor_delta_alpha=approx,
                       shape_chi=shape_chi,shape_time=shape_time,shape_matched=shape_hol,
                       kahler_delta_g_inverse2=Y_K,konishi_delta_g_inverse2=Y_Q)
        for n in range(len(chi)):
            all_rows.append(dict(model=model,**{k:float(v[n]) for k,v in columns.items()}))
        print(json.dumps(summary),flush=True)
    with (OUT/'matched_curves.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(all_rows[0]),lineterminator='\n')
        w.writeheader();w.writerows(all_rows)
    result = dict(created_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=time.time()-start,
                  scope='Specified local one-loop minimal U(1), not a complete visible-sector or quantum gravity prediction',
                  python=platform.python_version(),numpy=np.__version__,
                  source_hashes={str(p.relative_to(REPO)):sha(p) for p in
                      [PARENT/'results/inputs.json',PARENT/'results/trajectories.csv',ROOT/'protocol.md',Path(__file__)]},
                  summaries=summaries,checks=checks,checks_passed=sum(c['passed'] for c in checks),check_count=len(checks))
    write_json(OUT/'summary.json',result)
    print('Checks',result['checks_passed'],'/',result['check_count'],flush=True)
    if not all(c['passed'] for c in checks):raise SystemExit(1)


if __name__ == '__main__':
    main()
