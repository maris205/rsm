#!/usr/bin/env python3
"""Independent frozen-background reconstruction; never imports run_transfer.py."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import time
import platform
import numpy as np
import mpmath as mp
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
OUT=ROOT/'results'
REPORT=ROOT/'reports'
PHASES={'zero':(1.,0.),'quadrature':(0.,1.),'pi':(-1.,0.)}
MASSES=(0.,1e-33,1e-17,1e-6,1.,1e3)
INDICES=(0,512,1024)


def read_csv(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def dump_csv(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n')
        w.writeheader();w.writerows(rows)


def main():
    started=time.time(); REPORT.mkdir(parents=True,exist_ok=True)
    inputs_path=REPO/'experiments/sugra_shift_v01/results/inputs.json'
    trajectory_path=REPO/'experiments/sugra_shift_v01/results/trajectories.csv'
    inputs=json.loads(inputs_path.read_text())
    parent=np.genfromtxt(trajectory_path,names=True,delimiter=',',dtype=None,encoding='utf8')
    parent=parent[parent['case']=='shift_axis']
    chi=parent['chi']; ci=inputs['chi_i']; eps=inputs['epsilon']
    fs=inputs['F_squared_eV2']; P=inputs['Mpl_eV']; H=inputs['H_ref_eV']
    us=inputs['potential_unit_eV4']; mi=1e11; b=1-1.5*eps
    rho=3*inputs['omega_lambda']*H*H*P*P
    U=us*inputs['W_i']*np.exp(-2*(chi-ci)); oldforce=2*b*U
    gp=-ci**2/chi**3; g=.5*(1+(ci/chi)**2); gpp=3*ci**2/chi**4
    s=mi*mi*g; sp=mi*mi*gp; L=np.log(g)
    r=gp/(2*g); rp=(gpp/g-(gp/g)**2)/2
    b0=np.sqrt(2*U/fs)*(r+eps/2)
    b0p=np.sqrt(2*U/fs)*(rp-r-eps/2)
    checks=[]

    def check(name,value,tol=1e-8):
        checks.append(dict(name=name,value=float(value),tolerance=tol,
                           passed=bool(np.isfinite(value) and value<tol)))

    # Taylor expansion through total splitting order four, derived independently
    # around s. The main implementation's expression is not loaded or inspected.
    def stable(m,cp,st):
        C=np.sqrt(2*fs*U)*m
        d=m*m+rho/P**2+(1-eps)*U/P**2+2*C*cp/P**2
        dp=-2*(1-eps)*U/P**2-2*C*cp/P**2
        h2=s*((b0-m*cp)**2+(m*st)**2)
        h2p=sp*((b0-m*cp)**2+(m*st)**2)+2*s*(b0-m*cp)*b0p
        t1=4*s*d*(L-1)
        t1p=4*((sp*d+s*dp)*(L-1)+d*sp)
        t2=2*L*(d*d+h2)
        t2p=2*(sp/s*(d*d+h2)+L*(2*d*dp+h2p))
        pol3=d**3+3*d*h2
        pol3p=3*d*d*dp+3*dp*h2+3*d*h2p
        t3=2*pol3/(3*s)
        t3p=2/3*(pol3p/s-pol3*sp/s**2)
        pol4=d**4+6*d*d*h2+h2*h2
        pol4p=4*d**3*dp+12*d*dp*h2+6*d*d*h2p+2*h2*h2p
        t4=-pol4/(6*s*s)
        t4p=-pol4p/(6*s*s)+pol4*sp/(3*s**3)
        cw=(t1+t2+t3+t4)/(32*np.pi**2)
        cwp=(t1p+t2p+t3p+t4p)/(32*np.pi**2)
        tree=C*np.sqrt(9*cp*cp+st*st)/oldforce
        return dict(d=d,h2=h2,cw=cw,cwp=cwp,tree=tree,loop=abs(cwp)/oldforce)

    # Reconstruct every scan case over all 1025 old-background points.
    scan=read_csv(OUT/'scan.csv'); scan_errors={}; rebuilt=[]
    oldbudgets=json.loads((REPO/'experiments/protected_threshold_v01/results/budgets.json').read_text())
    oldsoft=next(v['msoft_for_force_point1_eV'] for v in oldbudgets['records'] if v['mass_eV']==mi and v['mu_over_mi']==1.)
    g2=4*np.pi*inputs['alpha_reference']; kappa=g2/(2*np.pi**2)*np.log(1e13/mi)
    for row in scan:
        phase=row['phase'];m=float(row['mG_eV']); cp,st=PHASES[phase]; a=stable(m,cp,st)
        va=(3*m*m*P*P+rho)**.25
        gl=g2/2*np.sqrt(3*m*m*P*P+rho)/P
        soft=np.sqrt(kappa)*gl
        rec=dict(phase=phase,mG_eV=m,
                 max_initial_axis_tree_force_ratio=float(max(a['tree'])),
                 max_tree_index=int(np.argmax(a['tree'])),
                 max_charged_CW_force_ratio=float(max(a['loop'])),
                 max_CW_index=int(np.argmax(a['loop'])),
                 max_d_over_s=float(max(a['d']/s)),max_h_over_s=float(max(np.sqrt(a['h2'])/s)),
                 min_scalar_ratio=float(min(1+a['d']/s-np.sqrt(a['h2'])/s)),
                 VA_marker_eV=va,VA_over_MI=va/mi,
                 direct_gaugino_cX1_eV=gl,gauge_common_msoft_cX1_eV=soft,
                 gauge_common_old_force_ratio_cX1=.1*(soft/oldsoft)**2)
        if m>0:
            rec['phase_offset_for_axis_longitudinal_budget']=float(min(.1*oldforce/(3*np.sqrt(2*fs*U)*m)))
        for key,value in rec.items():
            if key in ('phase','mG_eV'):continue
            expected=float(row[key])
            err=abs(value-expected)/max(abs(value),abs(expected),1e-300)
            scan_errors[key]=max(scan_errors.get(key,0.),err)
        flags={'heavy_threshold_below_marker':va>=mi,'heavy_threshold_below_marker_by_ten':va>=10*mi,
               'tree_axis_budget_pass':max(a['tree'])<=.1,'CW_budget_pass':max(a['loop'])<=.1}
        for key,value in flags.items():
            check('scan_flag_'+phase+'_'+str(m)+'_'+key,int(value!=(row[key]=='True')),.5)
        rebuilt.append(rec)
    for key,error in scan_errors.items():check('full_scan_'+key,error)
    dump_csv(OUT/'independent_scan.csv',[{k:r.get(k,'') for k in rebuilt[-1]} for r in rebuilt])

    # Fixed representative comparisons, with direct 180/240-digit supertraces.
    reps=read_csv(OUT/'representative_curves.csv')
    selected={(r['phase'],float(r['mG_eV']),int(r['index'])):r for r in reps if float(r['mG_eV']) in MASSES and int(r['index']) in INDICES}
    def constants():
        return {k:mp.mpf(str(v)) for k,v in inputs.items() if isinstance(v,(int,float))}

    def exact_spectrum(x,m,cp,st):
        ii=constants(); ee=ii['epsilon']; pp=ii['Mpl_eV']; ff=mp.sqrt(ii['F_squared_eV2'])
        uu=ii['potential_unit_eV4']*ii['W_i']*mp.exp(-2*(x-ii['chi_i']))
        gg=(1+(ii['chi_i']/x)**2)/2; ss=mp.mpf('1e22')*gg
        rr=-ii['chi_i']**2/(2*x**3*gg)
        rrho=3*ii['omega_lambda']*ii['H_ref_eV']**2*pp**2
        WW=-ff*mp.sqrt(uu)/mp.sqrt(2)+m*pp**2*(cp+1j*st)
        ffhidden2=3*m*m*pp*pp+rrho
        # d reconstructed from the unexpanded SUGRA diagonal expression.
        dd=(uu+ffhidden2)/pp**2-2*abs(WW)**2/pp**4
        MM=mp.sqrt(ss); MZ=mp.sqrt(2)/ff*MM*rr
        BB=mp.sqrt(uu)*MZ-MM*mp.conj(WW)/pp**2
        return ss,dd,abs(BB)

    def exact_cw(x,m,cp,st):
        ss,dd,hh=exact_spectrum(x,m,cp,st)
        def f(v):return v*v*(mp.log(v/mp.mpf('1e22'))-mp.mpf('1.5'))
        return (f(ss+dd+hh)+f(ss+dd-hh)-2*f(ss))/(32*mp.pi**2)

    direct=[]; high={}; repeated=0.; direct_errors={'CW_eV4':0.,'CW_chi_eV4':0.,'s_eV2':0.,'d_eV2':0.,'B_abs2_eV4':0.}
    for phase,(cp,st) in PHASES.items():
        for m in MASSES:
            group=[]
            for ix in INDICES:
                values=[]
                for precision in (180,240):
                    with mp.workdps(precision):
                        xx=mp.mpf(str(chi[ix])); mm=mp.mpf(str(m))
                        cv=exact_cw(xx,mm,cp,st)
                        cvp=mp.diff(lambda v:exact_cw(v,mm,cp,st),xx,addprec=40)
                        ss,dd,hh=exact_spectrum(xx,mm,cp,st)
                        values.append([cv,cvp,ss,dd,hh*hh])
                with mp.workdps(250):
                    repeated=max(repeated,max(float(abs(a-bb)/max(abs(bb),mp.mpf('1e-200'))) for a,bb in zip(*values)))
                row=selected[(phase,m,ix)]
                rec=dict(phase=phase,mG_eV=m,index=ix)
                for key,value in zip(direct_errors,values[-1]):
                    rec[key]=float(value)
                direct.append(rec);group.append((row,rec))
            for key in direct_errors:
                scale=max(abs(a[key]) for _,a in group)
                err=max(abs(float(old[key])-new[key])/max(scale,1e-300) for old,new in group)
                direct_errors[key]=max(direct_errors[key],err)
                check('direct_'+phase+'_'+str(m)+'_'+key,err,1e-9)
    check('direct_180_240_precision',repeated,1e-50)
    dump_csv(OUT/'independent_supertrace.csv',direct)

    # Full two-field potential from W,D_Z W; no compact-potential import.
    def direct_v(x,y,beta,cp,st):
        ii=constants(); ee=ii['epsilon']
        q=mp.sqrt(ii['W_i'])*mp.exp(-(x-ii['chi_i']))
        l=3*ii['omega_lambda']/ee
        phase=mp.mpc(cp,st)
        w=-q*mp.exp(-1j*y)/mp.sqrt(2)+beta/ee*phase
        dz=q*(1+1j*ee*y)*mp.exp(-1j*y)-1j*mp.sqrt(2)*beta*y*phase
        return mp.exp(ee*y*y)*(abs(dz)**2+3*beta*beta/ee+l-3*ee*abs(w)**2)

    def compact_v(x,y,beta,cp,st):
        ii=constants(); ee=ii['epsilon'];uu=ii['W_i']*mp.exp(-2*(x-ii['chi_i']))
        cy=mp.cos(y)*cp-mp.sin(y)*st; sy=mp.sin(y)*cp+mp.cos(y)*st
        return mp.exp(ee*y*y)*(3*ii['omega_lambda']/ee+uu*(1-3*ee/2+ee*ee*y*y)+2*beta*beta*y*y+mp.sqrt(2*uu)*beta*((3-2*ee*y*y)*cy+2*y*sy))

    potential_rows=[]; potential_error=0.
    with mp.workdps(240):
        for phase,(cp,st) in PHASES.items():
            for m in MASSES:
                for ix,yy in zip(INDICES,(0.,.013,-.017)):
                    xx=mp.mpf(str(chi[ix]));y=mp.mpf(str(yy)); beta=mp.mpf(str(m))/mp.mpf(str(H))
                    vals=[direct_v(xx,y,beta,cp,st),mp.diff(lambda t:direct_v(t,y,beta,cp,st),xx),mp.diff(lambda t:direct_v(xx,t,beta,cp,st),y)]
                    pred=[compact_v(xx,y,beta,cp,st),mp.diff(lambda t:compact_v(t,y,beta,cp,st),xx),mp.diff(lambda t:compact_v(xx,t,beta,cp,st),y)]
                    errs=[float(abs(a-bb)/max(abs(bb),mp.mpf(1))) for a,bb in zip(vals,pred)]
                    potential_error=max(potential_error,*errs)
                    potential_rows.append(dict(phase=phase,mG_eV=m,index=ix,y=yy,V_error=errs[0],Vx_error=errs[1],Vy_error=errs[2]))
    check('54_direct_KW_value_gradients',potential_error,1e-100)
    dump_csv(OUT/'independent_potential.csv',potential_rows)

    # Independently locate all nine local valleys using numerical derivatives of W.
    valley_rows=[]; valley_errors={}; max_valley_gradient=0.
    with mp.workdps(100):
        for row in read_csv(OUT/'valleys.csv'):
            xx=mp.mpf(row['chi']); beta=mp.mpf(row['beta']); ix=int(row['index'])
            uu=mp.mpf(str(inputs['W_i']))*mp.exp(-2*(xx-mp.mpf(str(ci))))
            lam=3*mp.mpf(str(inputs['omega_lambda']))/mp.mpf(str(eps))
            vasym=mp.sqrt(uu)/(2*mp.sqrt(2))
            fn=lambda vv:mp.diff(lambda yy:direct_v(xx,yy,beta,0,1),vv/beta)/beta
            vv=mp.findroot(fn,(vasym*mp.mpf('.9'),vasym*mp.mpf('1.1')),tol=mp.mpf('1e-70'))
            yy=vv/beta
            Vv=direct_v(xx,yy,beta,0,1)
            Vx=mp.diff(lambda t:direct_v(t,yy,beta,0,1),xx)
            Vyy=mp.diff(lambda t:direct_v(xx,t,beta,0,1),yy,2)
            residual=abs(fn(vv))/max(1,mp.sqrt(2*uu))
            max_valley_gradient=max(max_valley_gradient,float(residual))
            rec=dict(index=ix,beta=float(beta),v=float(vv),y=float(yy),v_asymptotic=float(vasym),
                     relative_y_asymptotic_error=float(abs(vv/vasym-1)),VminusLambda_over_U=float((Vv-lam)/uu),
                     Vx_over_U=float(Vx/uu),transverse_m2_over_mG2=float(Vyy/beta**2))
            for key,value in rec.items():
                if key in ('index','beta'):continue
                expected=float(row[key]);scale=max(abs(value),abs(expected),1e-10)
                err=abs(value-expected)/scale
                # Tiny asymptotic deviations are better measured absolutely.
                if key=='relative_y_asymptotic_error':err=abs(value-expected)
                valley_errors[key]=max(valley_errors.get(key,0.),err)
            valley_rows.append(rec)
    for key,err in valley_errors.items():check('valley_'+key,err,1e-8)
    check('valley_direct_gradient_residual',max_valley_gradient,1e-60)
    dump_csv(OUT/'independent_valleys.csv',valley_rows)

    summary=json.loads((OUT/'summary.json').read_text())
    roots=[]
    for rec in summary['limits']['phases']:
        phase=rec['phase'];cp,st=PHASES[phase]
        rootlog=brentq(lambda lm:np.log10(max(stable(10**lm,cp,st)['loop'])/.1),-20,-14,xtol=1e-12)
        looplimit=10**rootlog
        axislimit=float(min(.1*oldforce/(np.sqrt(2*fs*U)*np.sqrt(9*cp*cp+st*st))))
        check(phase+'_root_CW',abs(looplimit/rec['charged_CW_limit_eV']-1))
        check(phase+'_root_axis',abs(axislimit/rec['initial_axis_force_limit_eV']-1))
        roots.append(dict(phase=phase,CW_limit_eV=looplimit,axis_limit_eV=axislimit))
    paths=[Path(__file__),inputs_path,trajectory_path,ROOT/'protocol.md',OUT/'scan.csv',OUT/'representative_curves.csv',OUT/'valleys.csv',OUT/'summary.json']
    result=dict(created_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=time.time()-started,
                python=platform.python_version(),mpmath=mp.__version__,
                source_hashes={str(p.relative_to(REPO)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                scope='Frozen-background reconstruction only; no ODE or new observations. Low-VA heavy-loop results are formal extrapolations.',
                scan_cases=len(scan),background_points=len(chi),direct_supertrace_points=len(direct),
                direct_supertrace_precisions=[180,240],direct_KW_points=len(potential_rows),valley_points=len(valley_rows),
                max_scan_errors=scan_errors,max_direct_errors=direct_errors,precision_repeat_error=repeated,
                max_KW_error=potential_error,max_valley_errors=valley_errors,max_valley_gradient_residual=max_valley_gradient,
                roots=roots,checks=checks,checks_passed=sum(c['passed'] for c in checks),check_count=len(checks))
    (OUT/'independent_summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    report=f'''# 隐藏部门扫描：独立数值核查

## Material Passport

- 日期：{result['created_utc']}。
- 主实现 `run_transfer.py` 未读取、未导入或调用；根据协议与先前独立解析推导重建。
- 使用父实验的固定背景及本轮归档 CSV；未求解新的宇宙 ODE，未重新拟合观测。
- 命令：`timeout 180 python experiments/hidden_transfer_v01/code/independent_check.py`。
- 复现所需哈希、版本和逐项结果在 [独立摘要](../results/independent_summary.json)。本次耗时 {result['elapsed_seconds']:.3f} 秒。

## 核查范围与结果

1. 独立四阶分裂展开重建 {len(scan)} 组、每组 {len(chi)} 点的带电 CW、二维树力、质量比和域标记；同时检查约定的规范共同软项换算。摘要数字最大归一差 {max(scan_errors.values()):.3e}。
2. 在固定三相位、六质量和三背景点的 {len(direct)} 点上，以 180 和 240 位计算原始超迹差并直接数值微分，保留 d 的完整场依赖；两精度最大差 {repeated:.3e}。主数值与直接值的组峰值归一差：CW {direct_errors['CW_eV4']:.3e}，CW 导数 {direct_errors['CW_chi_eV4']:.3e}。接近零时按同相位、同质量三个采样点的峰值归一，没有删除点。
3. 直接由 W、D_Z W 构成完整两实场势，覆盖 {len(potential_rows)} 个轴上及轴外点；势与两个梯度最大归一差 {potential_error:.3e}。带电谱也通过未展开的 W 表达式重建。
4. 以完整 K/W 势的数值导数独立寻找全部 {len(valley_rows)} 个重场谷；最大归一梯度残差 {max_valley_gradient:.3e}。渐近式的偏离作为结果保留，未用渐近式替换数值根。
5. 三相位的 10% 力预算界均独立求根，与主结果一致。

合计 **{result['checks_passed']}/{result['check_count']}** 项数值核验通过。它们是确定性实现核查，数量不代表物理证据强度或外部同行评审。

## 适用范围

这轮只有旧背景上的预算与局部重场谷，不构成新的宇宙演化预测。低 nilpotent 强耦合标记区域的 100 GeV 阈值数字仅作形式延拓，不能称为该低能理论内受控的重粒子计算。独立数值一致性不补足 UV 完成、真实物质部门、完整高圈匹配或重整化条件。
'''
    (REPORT/'independent_results_cn.md').write_text(report)
    print(f"{result['checks_passed']}/{result['check_count']} checks; elapsed {result['elapsed_seconds']:.3f}s")
    print(json.dumps({k:result[k] for k in ('max_direct_errors','precision_repeat_error','max_KW_error','max_valley_errors','max_valley_gradient_residual')},indent=2))
    if result['checks_passed']!=result['check_count']:raise SystemExit(1)


if __name__=='__main__':main()
