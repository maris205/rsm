#!/usr/bin/env python3
"""Independent analytic-mode and high-precision output audit.
Does not read/import run_thresholds.py. Input physics comes from protocol.md and
parent input/old trajectory; main outputs are read only for comparison.
"""
from pathlib import Path
import csv, json, hashlib, math, sys
from datetime import datetime, timezone
import numpy as np
import mpmath as mp
from scipy.optimize import brentq

BASE=Path(__file__).resolve().parents[1]
ROOT=BASE.parents[1]
RES=BASE/'results'
PARENT=ROOT/'experiments/mediator_chain_v01/results/inputs.json'
TRAJ=ROOT/'experiments/sugra_shift_v01/results/trajectories.csv'
PARAM=json.loads(PARENT.read_text())['parameters']
M=lambda x:mp.mpf(str(x))
P=float(PARAM['P']); FC=math.sqrt(float(PARAM['F_squared_eV2']))
UI=float(PARAM['Ui']); CHII=float(PARAM['chi_i'])
with TRAJ.open() as f: OLD=[r for r in csv.DictReader(f) if r['case']=='shift_axis']
CHI=np.array([float(r['chi']) for r in OLD]); U=UI*np.exp(-2*(CHI-CHII))
CHECKS=[]

def check(name,ok,**detail):
    CHECKS.append(dict(name=name,passed=bool(ok),**detail))

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def chain_hp(n,dps):
    """Sine eigenmodes, no inversion or diagonalizer from the primary code."""
    with mp.workdps(dps):
        q=mp.mpf(3)
        angle=[mp.pi*j/(n+1) for j in range(1,n+1)]
        eig=[q*q+1-2*q*mp.cos(a) for a in angle]
        endpoint=[2*mp.sin(a)**2/(n+1) for a in angle]
        weight=[u/e for u,e in zip(endpoint,eig)]
        g=mp.fsum(weight); h2=mp.fsum(u/e**2 for u,e in zip(endpoint,eig))
        ell=[g*e for e in eig]
        # Hermite divided difference of F'. Off diagonal is evaluated in mp.
        lh=mp.mpf(0)
        for i in range(n):
            lh+=weight[i]**2*2*mp.log(ell[i])
            for j in range(i):
                dd=2*((ell[i]*mp.log(ell[i])-ell[j]*mp.log(ell[j]))/(ell[i]-ell[j])-1)
                lh+=2*weight[i]*weight[j]*dd
        out={'G00':g,'h2':h2,'a_hat':2*h2/g**2,'k_hat':1/(2*mp.pi**2),
             'h_hat':3*lh/(8*mp.pi**2*g**2),'LH_at_Lambda':lh,
             'tau':(q-q**-1)/(q**n-q**-n),
             'min_mass_over_Lambda':mp.sqrt(g*eig[0]),
             'max_mass_over_Lambda':mp.sqrt(g*eig[-1]),'endpoint_weight_sum':mp.fsum(weight)}
        return {k:mp.nstr(v,dps) for k,v in out.items()}

CHAINS={}
HP={}
for n in (128,127):
    lo=chain_hp(n,100); hi=chain_hp(n,150)
    with mp.workdps(150):
        err=max(abs(M(lo[k])/M(hi[k])-1) for k in lo)
    check(f'analytic_modes_precision_n{n}',err<mp.mpf('1e-80'),relative_error=str(err),precisions=[100,150],tolerance='1e-80')
    HP[str(n)]=hi;CHAINS[n]={k:float(v) for k,v in hi.items()}
    print('mode coefficients',n,hi['h_hat'],flush=True)

# Derived invariants checked with numerical differentiation, not hand-cancelled
# comparison expressions. One check groups 10 masses/chains x 6 scale routes.
rgerr=mp.mpf(0);deriverr=mp.mpf(0); chainruleerr=mp.mpf(0); wronggap=[]
with mp.workdps(100):
    for n,lam in ((128,mp.mpf('1e12')),(127,mp.mpf('2e12'))):
        a=M(HP[str(n)]['a_hat']);k=M(HP[str(n)]['k_hat']);h=M(HP[str(n)]['h_hat'])
        for mg in map(mp.mpf,('1e-17','1e-15','1e-14','1e-12','1e-6')):
            s=3*M(PARAM['P'])*mg**2; b=4/lam**2; mx=mp.sqrt(b*s)
            def fixed(z):return z*z*(-a+h+k*(mp.log(b*z/lam**2)-mp.mpf('1.5')))/lam**4
            def running(z):
                mu=mp.sqrt(b*z)
                cl=-a+h+2*k*mp.log(mu/lam)
                return z*z*(cl-mp.mpf('1.5')*k)/lam**4
            expected=2*s*(-a+h+k*(mp.log(b*s/lam**2)-1))/lam**4
            deriverr=max(deriverr,abs(mp.diff(fixed,s)/expected-1))
            chainruleerr=max(chainruleerr,abs(mp.diff(running,s)/expected-1))
            wronggap.append(float(abs((-a+h-k)-expected*lam**4/(2*s))))
            for mh,ml in ((lam/2,mx/10),(lam,mx),(2*lam,10*mx),(lam/2,lam),(2*lam,lam/2),(mp.mpf('1.25')*lam,3*mx)):
                cuv=-a+5*k*mp.log(mh/lam)
                heavy=h-3*k*mp.log(mh/lam)
                cl=cuv+heavy+2*k*mp.log(ml/mh)
                vl=s*s*(cl+k*(mp.log(b*s/ml**2)-mp.mpf('1.5')))/lam**4
                force=2*s*(cl+k*(mp.log(b*s/ml**2)-1))/lam**4
                rgerr=max(rgerr,abs(vl/fixed(s)-1),abs(force/expected-1))
    check('RG_path_invariance_high_precision',rgerr<mp.mpf('1e-80'),max_relative_error=str(rgerr),routes=60,tolerance='1e-80')
    check('fixed_scale_source_derivative',deriverr<mp.mpf('1e-80'),max_relative_error=str(deriverr),points=10,tolerance='1e-80')
    check('field_scale_total_derivative',chainruleerr<mp.mpf('1e-80'),max_relative_error=str(chainruleerr),points=10,tolerance='1e-80')
    check('reset_control_changes_physics',min(wronggap)>0.1,min_force_coefficient_difference=min(wronggap))

# Low-energy couplings from exact K_{X bar X}=1-4|X|²/L².
# X=x/sqrt2, V=S/(1-2x²/L²). d^4V/dx^4|0=96S/L⁴.
with mp.workdps(100):
    s=mp.mpf('1.7');L=mp.mpf('2.3')
    v=lambda x:s/(1-2*x*x/L**2)
    gfun=lambda x:4*mp.sqrt(s)*x/(L**2-4*x*x)
    err1=abs(mp.diff(v,0,4)/(96*s/L**4)-1)
    err2=abs(mp.diff(gfun,0)/(4*mp.sqrt(s)/L**2)-1)
    check('low_energy_scalar_and_Goldstino_couplings',max(err1,err2)<mp.mpf('1e-80'),scalar_quartic='96S/Lambda^4',Goldstino_yukawa_squared='16S/Lambda^4')

def profiles(n,lam,mg,phase='zero'):
    c=CHAINS[n]; a=c['a_hat'];k=c['k_hat'];h=c['h_hat']
    rw=FC*np.sqrt(U/2)/(P*mg)
    co,si=(1.,0.) if phase=='zero' else (0.,1.)
    s=3*P*mg*mg*((co-rw)**2+si*si)
    dx=6*P*mg*mg*rw*(co-rw)
    dy=2*P*mg*mg*rw*si
    grad=np.hypot(dx,dy)
    ell=np.log(4*s/lam**4)
    unit=s*grad/(U*lam**4)
    tree=a*unit;loop=np.abs(h+k*(ell-1))*unit;total=np.abs(-a+h+k*(ell-1))*unit
    scalar=6*s*np.abs(ell)/(math.pi**2*lam**4)
    return {'tree_ratio':tree,'loop_ratio':loop,'total_ratio':total,
            'scalar_log_parameter':scalar,'goldstino_log_parameter':scalar/6}

def benchmark(n,lam,mg):
    c=CHAINS[n];a=c['a_hat'];k=c['k_hat'];h=c['h_hat'];s=3*P*mg*mg
    mx=math.sqrt(4*s/lam**2);ell=math.log(4*s/lam**4)
    c0=-a+h;cmx=c0+k*ell
    pr=profiles(n,lam,mg)
    tr=float(np.max(pr['tree_ratio']));lr=float(np.max(pr['loop_ratio']));tot=float(np.max(pr['total_ratio']))
    return dict(mX_eV=mx,log_mX2_over_Lambda2=ell,cL_hat_at_Lambda=c0,cL_hat_at_mX=cmx,
       force_hat_at_fixed_Lambda=c0+k*(ell-1),force_hat_at_mX_with_chainrule=cmx-k,
       wrong_force_hat_with_reset=c0-k,wrong_force_hat_no_chainrule=cmx-math.fsum((k,k/2)),
       max_tree_force_ratio=tr,max_selected_loop_force_ratio=lr,max_selected_total_force_ratio=tot,
       selected_loop_over_new_tree=lr/tr,low_scalar_log_parameter=max(pr['scalar_log_parameter']),
       low_goldstino_log_parameter=max(pr['goldstino_log_parameter']))

# Decimal-preserving main-output comparison: relative tolerance 3e-11,
# coefficients near zero additionally allow 2e-14 absolute discrepancy.
# Report one grouped pass for each table, not a pass per numerical cell.
TABLE_SUMMARIES={}
def compare_table(filename,expect):
    with (RES/filename).open() as f:rows=list(csv.DictReader(f))
    errs=[];fail=[];cells=0
    for i,row in enumerate(rows):
        ex=expect(row,i)
        for key,value in ex.items():
            actual=row[key];cells+=1
            if isinstance(value,(bool,np.bool_)):
                ok=actual==str(bool(value));err=0 if ok else float('inf')
            else:
                actual=float(actual);value=float(value)
                scale=max(abs(actual),abs(value),1e-300)
                err=abs(actual-value)/scale
                abs_tol=2e-14 if key.endswith('_hat') or key=='ell' else 0.
                ok=math.isfinite(actual) and abs(actual-value)<=3e-11*scale+abs_tol
            errs.append(err)
            if not ok:fail.append(dict(row=i,column=key,actual=actual,expected=value,relative_error=err))
    table=dict(rows=len(rows),numeric_and_bool_cells=cells,max_relative_error=max(errs,default=0),failures=fail)
    TABLE_SUMMARIES[filename]=table
    check('table_'+filename,len(fail)==0,**table,tolerance='3e-11 relative, 2e-14 absolute for dimensionless *_hat / ell coefficients only')
    return rows

current_inputs=json.loads((RES/'inputs.json').read_text())
coef_err=[]
for row in current_inputs['chains']:
    n=int(row['n']);lam=float(row['Lambda_eV']);c=CHAINS[n]
    for key,value in c.items():
        if key=='min_mass_over_Lambda':key='min_mass_eV';value*=lam
        elif key=='max_mass_over_Lambda':key='max_mass_eV';value*=lam
        coef_err.append(abs(float(row[key])/value-1))
check('main_chain_coefficients',max(coef_err)<3e-11,max_relative_error=max(coef_err),fields=len(coef_err))

def scale_row(r,i):
    n=int(r['n']);lam=float(r['Lambda_eV']);mg=float(r['mG_eV']);mh=float(r['muH_eV']);ml=float(r['muL_eV'])
    c=CHAINS[n];a=c['a_hat'];k=c['k_hat'];h=c['h_hat'];s=3*P*mg*mg
    ell=math.log(4*s/(lam**2*ml**2));cu=-a+5*k*math.log(mh/lam)
    heavy=h-3*k*math.log(mh/lam);cl=cu+heavy+2*k*math.log(ml/mh)
    return dict(cUV_hat=cu,heavy_threshold_hat=heavy,cL_hat=cl,ell=ell,potential_hat=cl+k*(ell-1.5),derivative_hat=cl+k*(ell-1))
scale=compare_table('scale_audit.csv',scale_row)
check('main_reported_path_errors_small',max(float(r['normalized_error']) for r in scale)<1e-12)
benches=compare_table('benchmarks.csv',lambda r,i:benchmark(int(r['n']),float(r['Lambda_eV']),float(r['mG_eV'])))

def massrow(r,i):
    b=benchmark(127,2e12,float(r['mG_eV']))
    out={k:b[k] for k in ('max_selected_total_force_ratio','max_tree_force_ratio','max_selected_loop_force_ratio','selected_loop_over_new_tree','mX_eV')}
    out.update(max_scalar_log_parameter=b['low_scalar_log_parameter'],max_goldstino_log_parameter=b['low_goldstino_log_parameter'],selected_force_under_10percent=b['max_selected_total_force_ratio']<.1)
    return out
mass=compare_table('mass_scan.csv',massrow)
prcache={phase:profiles(127,2e12,1e-15,phase) for phase in ('zero','quadrature')}
def profilerow(r,i):
    ix=i%len(CHI);out={k:v[ix] for k,v in prcache[r['phase']].items() if k!='goldstino_log_parameter'}
    out.update(chi=CHI[ix],U_eV4=U[ix]);return out
profiles_main=compare_table('clock_profiles.csv',profilerow)
c=CHAINS[127];s0=3*P*1e-30

def curverow(r,i):
    f=float(r['mu_over_Lambda']);cl=-c['a_hat']+c['h_hat']+2*c['k_hat']*math.log(f)
    explicit=c['k_hat']*(math.log(4*s0/(2e12)**4)-2*math.log(f)-1)
    return dict(cL_hat=cl,explicit_loop_force_hat=explicit,total_force_hat=cl+explicit)
curve=compare_table('scale_curve.csv',curverow)

check('all_declared_output_row_counts',[len(scale),len(benches),len(mass),len(profiles_main),len(curve)]==[350,10,65,2050,141])
check('quadrature_y_factor_one_third',np.max(np.abs(prcache['quadrature']['total_ratio']/prcache['zero']['total_ratio']-1/3))<1e-13,
      scope='Local prescribed auxiliary-field gradient only; not a new rolling trajectory.')
logroot=brentq(lambda lm:benchmark(127,2e12,10**lm)['max_selected_total_force_ratio']-.1,-15,-13,xtol=1e-13)
root=10**logroot
main_summary=json.loads((RES/'summary.json').read_text())
check('conditional_mass_limit',abs(root/main_summary['selected_total_10percent_mG_eV']-1)<3e-11,
      independent_mG_eV=root,main_mG_eV=main_summary['selected_total_10percent_mG_eV'],tolerance='3e-11 relative')
# Independent mp root at final old coordinate. Every old-coordinate force was
# separately evaluated above; monotonicity over the finite input grid is checked.
check('finite_input_maximum_at_last_coordinate',all(np.argmax(profiles(127,2e12,float(r['mG_eV']))['total_ratio'])==len(CHI)-1 for r in mass))
with mp.workdps(100):
    lam=mp.mpf('2e12');pp=M(PARAM['P']);fc=mp.sqrt(M(PARAM['F_squared_eV2']));ui=M(PARAM['Ui'])
    uv=ui*mp.exp(-2*(M(OLD[-1]['chi'])-M(PARAM['chi_i'])))
    a=M(HP['127']['a_hat']);k=M(HP['127']['k_hat']);h=M(HP['127']['h_hat'])
    def last_force(lmg):
        mg=mp.power(10,lmg);rw=fc*mp.sqrt(uv/2)/(pp*mg);s=3*pp*mg**2*(1-rw)**2
        sx=6*pp*mg**2*rw*(1-rw)
        return s*sx*abs(-a+h+k*(mp.log(4*s/lam**4)-1))/(uv*lam**4)
    highroot=mp.power(10,mp.findroot(lambda z:last_force(z)-mp.mpf('.1'),(-14,-13)))
    check('high_precision_mass_limit',abs(M(root)/highroot-1)<mp.mpf('3e-11'),mG_eV=mp.nstr(highroot,90),tolerance='3e-11 relative for float comparison')

primary=benchmark(127,2e12,1e-15)
check('primary_interpretation',primary['max_selected_total_force_ratio']<.1 and primary['selected_loop_over_new_tree']>1 and primary['low_scalar_log_parameter']<1e-20 and primary['low_goldstino_log_parameter']<1e-20,
      detail='Loop/new-tree >1 is not the low-energy successive-loop parameter; unknown high-scale terms remain open.')
files=[PARENT,TRAJ,BASE/'protocol.md',Path(__file__),RES/'inputs.json',RES/'summary.json']+[RES/x for x in TABLE_SUMMARIES]
result=dict(created_utc=datetime.now(timezone.utc).isoformat(),checks_passed=sum(c['passed'] for c in CHECKS),check_count=len(CHECKS),
            checks=CHECKS,source_hashes={str(f.relative_to(ROOT)):sha(f) for f in files},
            main_source_read_or_imported=False,high_precision_chain_coefficients=HP,conditional_mass_limit_eV=root,
            primary=primary,table_comparisons=TABLE_SUMMARIES,
            scope=['Independent implementation of declared Gaussian one-loop leading S^2 sector only.',
                   'Does not independently establish the parent full-supergravity Sy formula; adopts its archived protocol and checks implementation of factor 1/3.',
                   'No higher loops, all-sector UV completion, cosmic trajectory, observational inference or Riemann origin established.'],
            failure_policy='Failures remain in this JSON; tolerance changes require preserving previous outputs.',
            grouped_check_policy='All CSV rows and declared numerical fields are compared; one pass per table rather than a pass per cell.')
output=RES/'independent_comparison.json'
if output.exists():
    old=json.loads(output.read_text())
    if old.get('checks_passed')!=old.get('check_count'):
        stamp=old['created_utc'].replace(':','').replace('+','_')
        (RES/f'independent_failure_{stamp}.json').write_text(json.dumps(old,ensure_ascii=False,indent=2)+'\n')
output.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
print(json.dumps({'checks_passed':result['checks_passed'],'check_count':len(CHECKS),'mass_limit_eV':root,'failed':[c['name'] for c in CHECKS if not c['passed']]},indent=2),flush=True)
sys.exit(0 if all(c['passed'] for c in CHECKS) else 1)
