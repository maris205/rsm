#!/usr/bin/env python3
"""Independent KL audit without importing or reading the primary implementation.

The scalar spectrum is reconstructed from the archived model definitions. All
3,075 curve points use 180-digit arithmetic; 18 points also use direct individual
mass logarithms. The output is a deterministic consistency audit, not a fit.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import platform
import time

import mpmath as mp
import numpy as np

EXPERIMENT = Path(__file__).resolve().parents[1]
REPO = EXPERIMENT.parents[1]
ARCHIVE = REPO / 'experiments/sugra_shift_v01/results'
OUT = EXPERIMENT / 'results'
REPORTS = EXPERIMENT / 'reports'
INPUTS = ARCHIVE / 'inputs.json'
OLD = ARCHIVE / 'trajectories.csv'
MAIN = OUT / 'matched_curves.csv'
MAIN_SUMMARY = OUT / 'summary.json'
PROTOCOL = EXPERIMENT / 'protocol.md'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    started = time.time()
    mp.mp.dps = 180
    raw = json.loads(INPUTS.read_text())
    p = {key: mp.mpf(str(value)) for key, value in raw.items()
         if isinstance(value, (int, float))}
    eps, ci, alpha = p['epsilon'], p['chi_i'], p['alpha_reference']
    mi, fs, mp2 = p['holomorphic_mass_i_eV'], p['F_squared_eV2'], p['Mpl_eV']**2
    pi = mp.pi
    old = list(csv.DictReader(OLD.open()))
    primary = list(csv.DictReader(MAIN.open()))
    summary = json.loads(MAIN_SUMMARY.read_text())
    checks, independent, endpoints, precision_samples = [], [], [], []

    def check(name, error, tolerance=mp.mpf('1e-10')):
        checks.append({'name':name,'error':mp.nstr(error,20),
                       'tolerance':mp.nstr(tolerance,10),
                       'passed':bool(mp.isfinite(error) and abs(error)<tolerance)})

    def boolean(name, condition):
        checks.append({'name':name,'passed':bool(condition)})

    def spectrum(model, chi):
        g = (1+(ci/chi)**2)/2
        mh2 = mi**2*g
        rr = -ci**2/(2*chi**3*g)
        U = p['potential_unit_eV4']*p['W_i']*mp.exp(-2*(chi-ci))
        if model == 'canonical':
            k = eps*chi**2/2
            a = 1-eps*chi/2
            dj = a*a-eps
            B = a*(rr+eps*chi/2)+eps/2
        else:
            k = mp.mpf(0)
            dj = 1-eps
            B = rr+(eps/2 if model == 'shift' else 0)
        s = mp.exp(k)*mh2
        d = mp.exp(k)*U*dj/mp2 if model != 'global' else mp.mpf(0)
        h = mp.exp(k)*mp.sqrt(2*U*mh2/fs)*abs(B)
        # Independent stable form: sum the two scalar logarithms, rather
        # than the primary's single polynomial log1p argument.
        logR = mp.log1p((d+h)/s)+mp.log1p((d-h)/s)
        return {'chi':chi,'k':k,'mh2':mh2,'s':s,'d':d,'h':h,'logR':logR}

    def mass_piece(v):
        return -(mp.log(v['s']+v['d']+v['h'])/3
                 +mp.log(v['s']+v['d']-v['h'])/3
                 +mp.mpf(4)/3*mp.log(v['s']))/(16*pi*pi)

    def transformed(v, kind):
        if kind == 'kahler':
            j = -eps*v['chi']**2/4 if model != 'global' else mp.mpf(0)
            kt = v['k']+2*j
            lm2 = mp.log(v['mh2'])-2*j
            zlog = mp.mpf(0)
            ft = -2*j/(8*pi*pi)
            wrong_shift = 2*j/(8*pi*pi)
        else:
            up = mp.mpf('.07')*(v['chi']-ci)
            um = -mp.mpf('.03')*(v['chi']-ci)
            u = up+um
            kt = v['k']
            lm2 = mp.log(v['mh2'])-2*u
            zlog = -2*u
            ft = -u/(4*pi*pi)
            wrong_shift = u/(4*pi*pi)
        st = mp.exp(kt+lm2-zlog)
        vt = dict(v, s=st)
        low = ft+kt/(8*pi*pi)-zlog/(8*pi*pi)+mass_piece(vt)
        wrong = kt/(8*pi*pi)-zlog/(8*pi*pi)+mass_piece(vt)
        return low, wrong, st, wrong_shift

    def stringify(row):
        return {key:mp.nstr(value,30) if isinstance(value,mp.mpf) else value
                for key,value in row.items()}

    for model in ('global','canonical','shift'):
        axis = [r for r in old if r['case']==model+'_axis']
        main_axis = [r for r in primary if r['model']==model]
        boolean(model+'_complete_1025_rows',len(axis)==len(main_axis)==1025)
        end = spectrum(model,mp.mpf(axis[-1]['chi']))
        beginning = spectrum(model,mp.mpf(axis[0]['chi']))
        final_mass_piece = mass_piece(end)
        final_transforms = {kind:transformed(end,kind) for kind in ('kahler','konishi')}
        chi_shape_norm = beginning['chi']**-2-end['chi']**-2
        time_end = mp.mpf(axis[-1]['t_Gyr'])*mp.mpf('1e9')*mp.mpf('365.25')*86400
        time_begin = mp.mpf(axis[0]['t_Gyr'])*mp.mpf('1e9')*mp.mpf('365.25')*86400
        time_end_log = mp.log(time_end/p['tstar_s'])
        time_norm = mp.log(time_begin/p['tstar_s'])**-2-time_end_log**-2
        rows, worst_transforms = [], {}
        qbegin = alpha/(2*pi)*mp.log(beginning['mh2']/end['mh2'])+alpha/(12*pi)*(beginning['logR']-end['logR'])
        abegin = qbegin/(1-qbegin)
        min_scalar_ratio = mp.mpf(1)
        for index,row in enumerate(axis):
            v = spectrum(model,mp.mpf(row['chi']))
            dk = v['k']-end['k']
            lm = mp.log(v['mh2']/end['mh2'])
            ls = mp.log(v['s']/end['s'])
            lr = v['logR']-end['logR']
            mass = mass_piece(v)-final_mass_piece
            anomaly = dk/(8*pi*pi)
            matched = -lm/(8*pi*pi)-lr/(48*pi*pi)
            q0 = alpha*lm/(2*pi)
            qr = alpha*lr/(12*pi)
            qphysical = alpha*ls/(2*pi)+qr
            aholo = (q0+qr)/(1-q0-qr)
            aphysical = qphysical/(1-qphysical)
            # Exact algebraic difference of two inversions of the same
            # one-loop inverse coupling, not an all-orders loop result.
            split = qr/((1-q0)*(1-q0-qr))
            X = (ci**2/v['chi']**2-ci**2/end['chi']**2)/(1+ci**2/end['chi']**2)
            qtaylor = alpha*X/(2*pi)+qr
            tt = mp.mpf(row['t_Gyr'])*mp.mpf('1e9')*mp.mpf('365.25')*86400
            result = {'model':model,'index':index,'N':row['N'],'z':row['z'],
                      'k':v['k'],'delta_k':dk,'delta_log_Mhol2':lm,'delta_log_s':ls,
                      'log_R':v['logR'],'delta_log_R':lr,
                      'mass_delta_g_inverse2':mass,'anomaly_delta_g_inverse2':anomaly,
                      'physical_uv_delta_f':-anomaly,
                      'matched_delta_g_inverse2':matched,
                      'holomorphic_uv_delta_alpha':aholo,
                      'physical_uv_delta_alpha':aphysical,'split_delta_alpha':split,
                      'taylor_delta_alpha':qtaylor/(1-qtaylor),
                      'shape_chi':(v['chi']**-2-end['chi']**-2)/chi_shape_norm,
                      'shape_time':(mp.log(tt/p['tstar_s'])**-2-time_end_log**-2)/time_norm,
                      'shape_matched':aholo/abegin}
            for kind in ('kahler','konishi'):
                trans, wrong, st, wrong_shift = transformed(v,kind)
                trans0, wrong0, st0, wrong_shift0 = final_transforms[kind]
                result[kind+'_delta_g_inverse2'] = trans-trans0
                result[kind+'_wrong_boundary_delta_g_inverse2'] = wrong-wrong0
                result[kind+'_negative_control_shift'] = wrong_shift-wrong_shift0
                result[kind+'_mass_invariance_error'] = (st-v['s'])/v['s']
                result[kind+'_negative_control_residual'] = (wrong-wrong0)-matched-(wrong_shift-wrong_shift0)
            result['high_low_identity_residual'] = mass+anomaly-matched
            result['scalar_min_over_s'] = (v['s']+v['d']-v['h'])/v['s']
            min_scalar_ratio = min(min_scalar_ratio,result['scalar_min_over_s'])
            rows.append(result)
            if index in (0,1,256,512,768,1024):
                directR = mp.log((v['s']+v['d']+v['h'])/v['s'])+mp.log((v['s']+v['d']-v['h'])/v['s'])
                relR = abs(directR-v['logR'])/max(abs(v['logR']),mp.mpf('1e-170'))
                check(f'{model}_{index}_direct_R_vs_stable',relR,mp.mpf('1e-80'))
                main_r = mp.mpf(main_axis[index]['log_R'])
                check(f'{model}_{index}_primary_R',abs(main_r-v['logR'])/max(abs(v['logR']),mp.mpf('1e-170')))
                main_split = mp.mpf(main_axis[index]['split_delta_alpha'])
                split_error = abs(main_split-split)/max(abs(split),mp.mpf('1e-170'))
                check(f'{model}_{index}_primary_exact_inversion_split',split_error)
                precision_samples.append(stringify({'model':model,'index':index,
                    'direct_log_R':directR,'stable_log_R':v['logR'],
                    'split_exact_inversion':split,'split_first_order':qr,
                    'exact_vs_first_relative':(split/qr-1) if qr else mp.mpf(0)}))
        keys = ['k','delta_k','delta_log_Mhol2','delta_log_s','log_R','delta_log_R',
                'mass_delta_g_inverse2','anomaly_delta_g_inverse2','physical_uv_delta_f',
                'matched_delta_g_inverse2','holomorphic_uv_delta_alpha',
                'physical_uv_delta_alpha','split_delta_alpha','taylor_delta_alpha',
                'shape_chi','shape_time','shape_matched',
                'kahler_delta_g_inverse2','konishi_delta_g_inverse2']
        for key in keys:
            expected = [r[key] for r in rows]
            peak = max(max(abs(x) for x in expected),mp.mpf('1e-170'))
            error = max(abs(mp.mpf(r[key])-v) for r,v in zip(main_axis,expected))/peak
            check(model+'_all_points_'+key,error)
        boolean(model+'_positive_scalar_spectrum',min_scalar_ratio>0)
        peakY = max(abs(r['matched_delta_g_inverse2']) for r in rows)
        check(model+'_full_high_low_identity',max(abs(r['high_low_identity_residual']) for r in rows)/peakY,mp.mpf('1e-170'))
        for kind in ('kahler','konishi'):
            check(model+'_'+kind+'_physical_mass',max(abs(r[kind+'_mass_invariance_error']) for r in rows),mp.mpf('1e-170'))
            check(model+'_'+kind+'_coupling',max(abs(r[kind+'_delta_g_inverse2']-r['matched_delta_g_inverse2']) for r in rows)/peakY,mp.mpf('1e-170'))
            check(model+'_'+kind+'_negative_control_formula',max(abs(r[kind+'_negative_control_residual']) for r in rows)/peakY,mp.mpf('1e-170'))
            ratio = max(abs(r[kind+'_wrong_boundary_delta_g_inverse2']-r['matched_delta_g_inverse2']) for r in rows)/peakY
            if kind!='kahler' or model!='global':
                boolean(model+'_'+kind+'_wrong_boundary_detected',ratio>mp.mpf('.1'))
            worst_transforms[kind+'_wrong_boundary_relative_peak'] = ratio
        # Direct differentiation of the analytic spectrum and logarithms.
        def qchi(chi, physical=False):
            v = spectrum(model,chi)
            lm = mp.log((v['s'] if physical else v['mh2'])/
                        (end['s'] if physical else end['mh2']))
            return alpha*lm/(2*pi)+alpha*(v['logR']-end['logR'])/(12*pi)
        chi_dot_year = mp.mpf(axis[-1]['wx'])*p['H_ref_s_inv']*mp.mpf('365.25')*86400
        drift_holo = mp.diff(lambda x:qchi(x),end['chi'])*chi_dot_year
        drift_physical = mp.diff(lambda x:qchi(x,True),end['chi'])*chi_dot_year
        primary_summary = next(s for s in summary['summaries'] if s['model']==model)
        for label,value in [('holomorphic_uv_today_drift_per_year',drift_holo),
                            ('physical_uv_today_drift_per_year',drift_physical)]:
            check(model+'_'+label,abs(mp.mpf(str(primary_summary[label]))-value)/abs(value))
        endpoints.append(stringify({'model':model,'points':len(rows),
            'holomorphic_uv_zi_ppm':rows[0]['holomorphic_uv_delta_alpha']*10**6,
            'physical_uv_zi_ppm':rows[0]['physical_uv_delta_alpha']*10**6,
            'split_exact_inversion_zi':rows[0]['split_delta_alpha'],
            'split_first_order_zi':alpha*(beginning['logR']-end['logR'])/(12*pi),
            'holomorphic_uv_today_drift_per_year':drift_holo,
            'physical_uv_today_drift_per_year':drift_physical,**worst_transforms}))
        independent.extend(stringify(r) for r in rows)
        print(model, 'points',len(rows),'completed',flush=True)
    OUT.mkdir(exist_ok=True,parents=True)
    REPORTS.mkdir(exist_ok=True,parents=True)
    with (OUT/'independent_curves.csv').open('w',newline='') as handle:
        writer = csv.DictWriter(handle,fieldnames=list(independent[0]),lineterminator='\n')
        writer.writeheader()
        writer.writerows(independent)
    result = {'created_utc':datetime.now(timezone.utc).isoformat(),
              'elapsed_seconds':time.time()-started,'precision_digits':mp.mp.dps,
              'method':'Independent analytic spectrum, direct physical mass logs and KL anomalies; no primary module imported or read',
              'scope':'Frozen old tree trajectories; 3075 complete points, 18 high-precision soft samples; no data fit or background rerun',
              'python':platform.python_version(),'mpmath':mp.__version__,'numpy':np.__version__,
              'source_hashes':{str(path.relative_to(REPO)):sha(path) for path in
                              (INPUTS,OLD,MAIN,MAIN_SUMMARY,PROTOCOL,Path(__file__).resolve())},
              'curve_points':len(independent),'precision_samples':precision_samples,
              'endpoints':endpoints,'checks':checks,
              'checks_passed':sum(c['passed'] for c in checks),'check_count':len(checks)}
    (OUT/'independent_checks.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    lines = ['# 规范匹配全曲线独立复核','',
        '本复核直接读取旧树级轨迹及冻结输入，重新构造解析物理谱；未调用或读取主匹配程序。'
        '三个模型各 1025 点全部保留，共 3075 点。全曲线运算使用 180 位精度，另在预先声明的 18 点直接计算三个质量对数和微小分裂。',
        '',
        f"结果：**{result['checks_passed']}/{result['check_count']} 项通过**。"
        '两种 UV 边界、组成项、总量、Kähler 与 Konishi 的真正变换及漏掉规范函数变换的负对照、今天漂移、微小分裂和形状曲线均已核查。',
        '',
        '| 模型 | 固定全纯 UV 端点（ppm） | 固定物理 UV 端点（ppm） | 今天全纯边界漂移（每年） |',
        '|---|---:|---:|---:|']
    for rec in endpoints:
        lines.append(f"| {rec['model']} | {float(rec['holomorphic_uv_zi_ppm']):+.9f} | {float(rec['physical_uv_zi_ppm']):+.9f} | {float(rec['holomorphic_uv_today_drift_per_year']):+.8e} |")
    lines.extend(['',
        'Kähler 非平凡变换只用于局域超对称的 canonical 与 shift；global 为刚性参考，取 j=0，其该项没有负对照识别力。Konishi 重定义在三者中均作非平凡测试。'
        '对所有适用的非平凡变换，正确变换 f 保持物理量不变，故意不变换 f 则产生可识别偏差。',
        '',
        '微小分裂另存两种数：一阶读数 δq，以及对同一个一环逆耦合表达式作代数求逆后得到的精确差 '
        'δq/[(1−q₀)(1−q₀−δq)]。后者包含求逆的代数因子，不称为全阶量子计算。'
        '它与之前独立报告的一阶读数有数个 10⁻⁵ 的相对差，这与 ppm 级基准 q₀ 一致，不能误报为公式或程序不一致。',
        '',
        '今天漂移用高精度直接对谱函数求 χ 导数，再乘旧轨迹今天的物理 χ 速度；包含局部微小分裂项。'
        '模型间的共同今天 α₀ 归一化，不等于同一数值 UV 耦合。',
        '',
        '输入、协议及两条计算代码输出的 SHA-256 见 [独立检查 JSON](../results/independent_checks.json)；'
        '完整独立曲线见 [CSV](../results/independent_curves.csv)。'
        '复现命令：timeout 180 python experiments/gauge_matching_v01/code/independent_check.py。',
        '',
        '这是内部确定性一致性检验，不是观测显著性、外部同行评审或真实标准模型常数的完整计算。'
        '复用的背景与参数未修改，也未扩大 KL 与弱破缺局部近似的适用范围。',''])
    (REPORTS/'independent_results_cn.md').write_text('\n'.join(lines))
    print('Checks',result['checks_passed'],'/',result['check_count'],flush=True)
    if result['checks_passed'] != result['check_count']:
        for check in checks:
            if not check['passed']:
                print('FAILED',check,flush=True)
        raise SystemExit(1)


if __name__=='__main__':
    main()
