#!/usr/bin/env python3
"""Minimal nilpotent hidden-sector transfer: frozen-background diagnostics only."""
from pathlib import Path
from datetime import datetime, timezone
import csv, hashlib, json, platform, time
import numpy as np
import scipy
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
PARENT=REPO/'experiments/sugra_shift_v01'
OUT=ROOT/'results'
PHASES={'zero':(1.,0.),'quadrature':(0.,1.),'pi':(-1.,0.)}
MASSES=np.r_[0.,np.logspace(-36,3,157)]
REPRESENTATIVES=(0.,1e-33,1e-17,1e-6,1e-3,1.,1e3)
MI=1e11


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def write_csv(p,rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n')
        w.writeheader();w.writerows(rows)


def read_input():
    i=json.loads((PARENT/'results/inputs.json').read_text())
    with (PARENT/'results/trajectories.csv').open() as f:
        rows=[r for r in csv.DictReader(f) if r['case']=='shift_axis']
    r={k:np.array([float(v[k]) for v in rows]) for k in rows[0] if k not in ('case','model')}
    i['rho_lambda_eV4']=3*i['omega_lambda']*i['Mpl_eV']**2*i['H_ref_eV']**2
    return i,r


def spectrum(chi,mG,phase,i):
    eps=i['epsilon'];fs=i['F_squared_eV2'];mp2=i['Mpl_eV']**2
    xi=i['xi'];ci=i['chi_i'];co,si=PHASES[phase]
    U=i['potential_unit_eV4']*i['W_i']*np.exp(-2*(chi-ci))
    g=.5*(1+xi/chi**2);s=MI**2*g
    r=-xi/(chi*(chi**2+xi));rp=xi*(3*chi**2+xi)/(chi**2*(chi**2+xi)**2)
    sp=2*r*s
    a=np.sqrt(2*U/fs)*(r+eps/2)
    ap=np.sqrt(2*U/fs)*(rp-r-eps/2)
    C=np.sqrt(2*fs*U)*mG
    d=mG**2+(1-eps)*U/mp2+i['rho_lambda_eV4']/mp2+2*C*co/mp2
    dp=-2*(1-eps)*U/mp2-2*C*co/mp2
    h2=s*((a-mG*co)**2+(mG*si)**2)
    h2p=sp*((a-mG*co)**2+(mG*si)**2)+2*s*(a-mG*co)*ap
    return dict(U=U,s=s,sp=sp,d=d,dp=dp,h2=h2,h2p=h2p,C=C,
                extra_chi_force=-3*C*co,extra_y_force=-C*si,
                tree_force=2*(1-1.5*eps)*U)


def loop_budget(r):
    s,sp,d,dp,h2,h2p=(r[k] for k in ('s','sp','d','dp','h2','h2p'))
    L=np.log(s/MI**2);t=s+d;Lt=L+np.log1p(d/s)
    common=4*s*d*(L-1)+2*d*d*L+(2/3)*d**3/s
    commonp=4*(s*dp*(L-1)+d*sp*L)+4*d*dp*L+2*d*d*sp/s+2*d*d*dp/s-(2/3)*d**3*sp/s**2
    split=2*h2*Lt-h2*h2/(6*t*t)
    splitp=2*h2p*Lt+2*h2*(sp+dp)/t-h2*h2p/(3*t*t)+h2*h2*(sp+dp)/(3*t**3)
    return (common+split)/(32*np.pi**2),(commonp+splitp)/(32*np.pi**2)


def potential(chi,y,beta,phase,i):
    """Dimensionless V/(Fchi² Href²), with Lambda included exactly once."""
    eps=i['epsilon'];b=1-1.5*eps;U=i['W_i']*np.exp(-2*(chi-i['chi_i']))
    lam=3*i['omega_lambda']/eps
    co,si=PHASES[phase]
    ct=np.cos(y)*co-np.sin(y)*si;st=np.sin(y)*co+np.cos(y)*si
    C=np.sqrt(2*U)*beta
    A=(3-2*eps*y*y)*ct+2*y*st
    Ay=2*(1-2*eps)*y*ct+(-1+2*eps*y*y)*st
    Ayy=(1-4*eps+2*eps*y*y)*ct+(-2+8*eps)*y*st
    rest=U*(b+eps**2*y*y)+2*beta**2*y*y+C*A
    P=lam+rest;Py=2*eps**2*U*y+4*beta**2*y+C*Ay
    Pyy=2*eps**2*U+4*beta**2+C*Ayy
    ex=np.exp(eps*y*y)
    return dict(V=ex*P,VminusLambda=np.expm1(eps*y*y)*lam+ex*rest,
                Vx=ex*(-2*U*(b+eps**2*y*y)-C*A),
                Vy=ex*(Py+2*eps*y*P),
                Vyy=ex*(Pyy+4*eps*y*Py+(2*eps+4*eps**2*y*y)*P),U=U)


def main():
    start=time.time();OUT.mkdir(parents=True,exist_ok=True)
    i,old=read_input();chi=old['chi'];eps=i['epsilon'];b=1-1.5*eps
    checks=[];scan=[];curves=[];valleys=[]
    def check(name,value,tol=1e-8):
        checks.append(dict(name=name,value=float(value),tolerance=tol,
                           passed=bool(np.isfinite(value) and value<tol)))
    threshold=MI**2/(np.sqrt(3)*i['Mpl_eV'])
    kappa=8*(4*np.pi*i['alpha_reference'])*np.log(1e13/MI)/(16*np.pi**2)
    old_budgets=json.loads((REPO/'experiments/protected_threshold_v01/results/budgets.json').read_text())
    db=[x for x in old_budgets['records'] if x['mass_eV']==MI and x['mu_over_mi']==1.][0]
    d_limit=db['msoft_for_force_point1_eV']**2

    def evaluate(mG,phase):
        r=spectrum(chi,mG,phase,i);u1,u1p=loop_budget(r)
        force_tree=np.hypot(r['extra_chi_force'],r['extra_y_force'])/r['tree_force']
        force_loop=abs(u1p)/r['tree_force']
        cutoff=(3*mG*mG*i['Mpl_eV']**2+i['rho_lambda_eV4'])**.25
        fx=np.sqrt(3*mG*mG*i['Mpl_eV']**2+i['rho_lambda_eV4'])
        mlambda=(4*np.pi*i['alpha_reference']/2)*fx/i['Mpl_eV']
        rec=dict(phase=phase,mG_eV=float(mG),beta=mG/i['H_ref_eV'],
                 max_initial_axis_tree_force_ratio=float(np.max(force_tree)),
                 max_tree_index=int(np.argmax(force_tree)),max_charged_CW_force_ratio=float(np.max(force_loop)),
                 max_CW_index=int(np.argmax(force_loop)),
                 max_d_over_s=float(np.max(abs(r['d']/r['s']))),
                 max_h_over_s=float(np.max(np.sqrt(r['h2'])/r['s'])),
                 min_scalar_ratio=float(np.min(1+r['d']/r['s']-np.sqrt(r['h2'])/r['s'])),
                 VA_marker_eV=float(cutoff),VA_over_MI=float(cutoff/MI),
                 heavy_threshold_below_marker=bool(cutoff>=MI),
                 heavy_threshold_below_marker_by_ten=bool(cutoff>=10*MI),
                 tree_axis_budget_pass=bool(np.max(force_tree)<=.1),
                 CW_budget_pass=bool(np.max(force_loop)<=.1),
                 direct_gaugino_cX1_eV=float(mlambda),
                 gauge_common_msoft_cX1_eV=float(np.sqrt(kappa)*mlambda),
                 gauge_common_old_force_ratio_cX1=float(.1*kappa*mlambda**2/d_limit),
                 phase_offset_for_axis_longitudinal_budget=float(.1/np.max(3*r['C']/r['tree_force'])) if mG>0 else None)
        return rec,r,u1,u1p,force_tree,force_loop

    limits=[]
    for phase in PHASES:
        for mG in MASSES:scan.append(evaluate(float(mG),phase)[0])
        for mG in REPRESENTATIVES:
            rec,r,u1,u1p,ft,fl=evaluate(mG,phase)
            for n in range(len(chi)):
                curves.append(dict(phase=phase,mG_eV=mG,index=n,z=float(old['z'][n]),chi=float(chi[n]),
                    U_eV4=float(r['U'][n]),s_eV2=float(r['s'][n]),d_eV2=float(r['d'][n]),
                    B_abs2_eV4=float(r['h2'][n]),CW_eV4=float(u1[n]),CW_chi_eV4=float(u1p[n]),
                    tree_axis_force_ratio=float(ft[n]),charged_CW_force_ratio=float(fl[n])))
        coeff=evaluate(1.,phase)[0]['max_initial_axis_tree_force_ratio']
        axis_limit=.1/coeff
        rootlog=brentq(lambda v:np.log10(evaluate(10**v,phase)[0]['max_charged_CW_force_ratio'])+1,-30,-5,xtol=1e-12)
        loop_limit=10**rootlog
        check(phase+'_loop_limit_equation',abs(evaluate(loop_limit,phase)[0]['max_charged_CW_force_ratio']/.1-1))
        check(phase+'_axis_limit_equation',abs(evaluate(axis_limit,phase)[0]['max_initial_axis_tree_force_ratio']/.1-1))
        phase_rows=[r for r in scan if r['phase']==phase]
        check(phase+'_scalar_positive',0 if all(r['min_scalar_ratio']>0 for r in phase_rows) else 1,.5)
        check(phase+'_small_split_domain',max(max(r['max_d_over_s'],r['max_h_over_s']) for r in phase_rows),.01)
        # Large-mG formula is a diagnostic, independently evaluated with the exact finite split.
        r=spectrum(chi,1.,phase,i);_,p=loop_budget(r)
        leading=r['sp']*(3*np.log(r['s']/MI**2)+1)/(16*np.pi**2)
        check(phase+'_large_hidden_CW_limit',np.max(abs(p-leading))/np.max(abs(leading)))
        limits.append(dict(phase=phase,initial_axis_force_limit_eV=axis_limit,
                           charged_CW_limit_eV=loop_limit,
                           valley_leading_CW_limit_eV=loop_limit*np.sqrt((b-.25)/b) if phase=='quadrature' else None,
                           heavy_EFT_marker_over_loop_limit=threshold/loop_limit))

    for n in (0,512,1024):
        for beta in (1e2,1e4,1e6):
            U=i['W_i']*np.exp(-2*(chi[n]-i['chi_i']))
            va=np.sqrt(2*U)/4
            fun=lambda v:potential(chi[n],v/beta,beta,'quadrature',i)['Vy']/beta
            v=brentq(fun,0,2*va,xtol=1e-14);y=v/beta
            p=potential(chi[n],y,beta,'quadrature',i)
            residual=abs(fun(v))/max(np.sqrt(U),1.)
            check(f'valley_{n}_{beta:g}_gradient',residual,1e-9)
            check(f'valley_{n}_{beta:g}_positive_curvature',0 if p['Vyy']>0 else 1,.5)
            valleys.append(dict(index=n,z=float(old['z'][n]),chi=float(chi[n]),beta=beta,
                mG_eV=beta*i['H_ref_eV'],v=float(v),y=float(y),v_asymptotic=float(va),
                relative_y_asymptotic_error=float(abs(v/va-1)),
                VminusLambda_over_U=float(p['VminusLambda']/U),expected_coefficient=b-.25,
                Vx_over_U=float(p['Vx']/U),expected_Vx_coefficient=-2*(b-.25),
                gradient_residual=float(residual),transverse_m2_over_mG2=float(p['Vyy']/beta**2),
                VA_over_MI=float((3*(beta*i['H_ref_eV'])**2*i['Mpl_eV']**2+i['rho_lambda_eV4'])**.25/MI)))
    # rhoLambda moved into hidden F once: old on-axis potential and background unchanged at mG=0.
    baseline=potential(chi,0.,0.,'zero',i)
    check('zero_mG_axis_potential_matches_old_with_lambda',
          np.max(abs(baseline['V']-(old['V']+3*i['omega_lambda']/eps)))/np.max(abs(baseline['V'])))
    check('zero_mG_axis_gradient_matches_old',np.max(abs(baseline['Vx']-old['Vx']))/np.max(abs(old['Vx'])))
    limits_out=dict(rho_lambda_eV4=i['rho_lambda_eV4'],VA_marker_floor_eV=i['rho_lambda_eV4']**.25,
                   mG_for_VA_MI_eV=threshold,mG_for_VA_10MI_eV=100*threshold,
                   marker_scope='Order-one power-counting scale, necessary diagnostic only; unknown UV completion.',
                   phases=limits,gauge_kappa=kappa,gauge_msoft_over_gaugino=np.sqrt(kappa),
                   old_common_msoft_limit_eV=np.sqrt(d_limit),
                   gauge_gaugino_limit_from_old_common_eV=np.sqrt(d_limit/kappa),
                   direct_cX1_mG_limit_from_old_common_eV=np.sqrt(d_limit/kappa)/(2*np.pi*i['alpha_reference']*np.sqrt(3)),
                   no_grid_overlap_of_heavy_marker_and_CW_budget=not any(r['heavy_threshold_below_marker'] and r['CW_budget_pass'] for r in scan))
    write_csv(OUT/'scan.csv',scan);write_csv(OUT/'representative_curves.csv',curves);write_csv(OUT/'valleys.csv',valleys)
    dump(OUT/'inputs.json',dict(parent_inputs=i,masses_eV=MASSES.tolist(),phases=PHASES,representatives_eV=REPRESENTATIVES,
         diagnostics={'force_limit':.1,'matching_mu_eV':MI,'heavy_mass_initial_eV':MI,'UV_gauge_scale_eV':1e13}))
    result=dict(created_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=time.time()-start,
        python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
        source_hashes={str(p.relative_to(REPO)):sha(p) for p in
          [PARENT/'results/inputs.json',PARENT/'results/trajectories.csv',REPO/'experiments/protected_threshold_v01/results/budgets.json',ROOT/'protocol.md',Path(__file__)]},
        scan_rows=len(scan),representative_rows=len(curves),valley_rows=len(valleys),limits=limits_out,
        checks=checks,checks_passed=sum(x['passed'] for x in checks),check_count=len(checks),
        scope='Frozen-background diagnostic; low-VA heavy loops are formal extrapolations, not a controlled UV calculation.')
    dump(OUT/'summary.json',result)
    print(json.dumps(limits_out,indent=2),flush=True)
    print('Checks',result['checks_passed'],'/',result['check_count'],flush=True)
    if not all(x['passed'] for x in checks):raise SystemExit(1)


if __name__=='__main__':main()
