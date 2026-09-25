#!/usr/bin/env python3
"""Independent exact no-scale audit. No stabilization model or new FLRW solution."""
from pathlib import Path
import csv, hashlib, json, platform, time
import sympy as sp
import numpy as np
import mpmath as mp

HERE=Path(__file__).resolve().parent
REPO=HERE.parent
records=[]
started=time.time()

def symbolic(name, expr):
    residual=expr.applyfunc(sp.simplify) if isinstance(expr,sp.MatrixBase) else sp.simplify(expr)
    passed=all(x==0 for x in residual) if isinstance(residual,sp.MatrixBase) else residual==0
    records.append(dict(name=name, method='symbolic', passed=passed,
                        residual=str(residual)))

def numeric(name, residual, tol):
    records.append(dict(name=name, method='numeric', passed=bool(abs(residual)<tol),
                        residual=str(residual), tolerance=str(tol)))

P,Y,a,b,W,Wb,w,wb=sp.symbols('P Y a b W Wb w wb', nonzero=True)
g=sp.Matrix([[3*P/Y**2,-b/Y**2],[-a/Y**2,1/Y+a*b/(3*P*Y**2)]])
ki=sp.Matrix([-3*P/Y,a/Y]); kib=sp.Matrix([-3*P/Y,b/Y])
D=sp.Matrix([-3*W/Y,w+a*W/(P*Y)])
Db=sp.Matrix([-3*Wb/Y,wb+b*Wb/(P*Y)])
# Metric rows are holomorphic, columns antiholomorphic; contraction uses inverse transpose.
symbolic('no_scale_norm', (ki.T*g.inv().T*kib)[0]-3*P)
symbolic('full_W_independent_T_potential', ((D.T*g.inv().T*Db)[0]-3*W*Wb/P)/Y**3-w*wb/Y**2)
symbolic('F_matter', -(g.inv().T*Db)[1]/Y**sp.Rational(3,2)+wb/sp.sqrt(Y))
symbolic('F_modulus', -(g.inv().T*Db)[0]/Y**sp.Rational(3,2)-(Wb-a*wb/3)/(P*sp.sqrt(Y)))

# Three arbitrary flat-Hessian matter coordinates (Z,Q+,Q-), independent k_i and k_bar_i.
avec=sp.Matrix(sp.symbols('a0:3')); bvec=sp.Matrix(sp.symbols('b0:3'))
g4=sp.zeros(4); g4[0,0]=3*P/Y**2
g4[0,1:4]=-bvec.T/Y**2; g4[1:4,0]=-avec/Y**2
g4[1:4,1:4]=sp.eye(3)/Y+avec*bvec.T/(3*P*Y**2)
gi=sp.zeros(4); gi[0,0]=Y**2/(3*P)+Y*(bvec.T*avec)[0]/(9*P**2)
gi[0,1:4]=Y*bvec.T/(3*P); gi[1:4,0]=Y*avec/(3*P)
gi[1:4,1:4]=Y*sp.eye(3)
symbolic('four_field_inverse', g4*gi-sp.eye(4))
wvec=sp.Matrix(sp.symbols('w0:3')); wbvec=sp.Matrix(sp.symbols('wb0:3'))
D4=sp.Matrix([-3*W/Y,*list(wvec+avec*W/(P*Y))])
Db4=sp.Matrix([-3*Wb/Y,*list(wbvec+bvec*Wb/(P*Y))])
symbolic('four_field_full_potential', (D4.T*gi.T*Db4)[0]-3*W*Wb/P-Y*(wvec.T*wbvec)[0])

U,M,Mb,Mz,Mzb,t,qp,qpb,qm,qmb=sp.symbols('U M Mb Mz Mzb t qp qpb qm qmb', nonzero=True)
scale=sp.symbols('scale')
sq=qp*qpb+qm*qmb
VV=(U+wb*Mz*qp*qm+w*Mzb*qpb*qmb+M*Mb*sq)/(t-sq/(3*P))**2
expanded=sp.diff(VV.subs({qp:scale*qp,qpb:scale*qpb,qm:scale*qm,qmb:scale*qmb}),scale,2).subs(scale,0)/2
expected=(M*Mb/t**2+2*U/(3*P*t**3))*sq+(wb*Mz*qp*qm+w*Mzb*qpb*qmb)/t**2
symbolic('charged_quadratic_full_denominator', expanded-expected)
symbolic('canonical_common_d', t*(M*Mb/t**2+2*U/(3*P*t**3))-M*Mb/t-2*U/(3*P*t**2))
symbolic('canonical_B', t*wb*Mz/t**2-wb*Mz/t)
symbolic('canonical_fermion_mass_squared', (t*t/t**3)*M*Mb-M*Mb/t)

F,chi,y,r,ui,ci=sp.symbols('F chi y r ui ci', real=True, positive=True)
eps=F**2/P
YY=t-eps*y**2/3
V0=ui*sp.exp(-2*(chi-ci))/YY**2
symbolic('axis_transverse_gradient', sp.diff(V0,y).subs(y,0))
symbolic('axis_transverse_mass', sp.diff(V0,y,2).subs(y,0)*t/F**2-4*ui*sp.exp(-2*(chi-ci))/(3*P*t**2))
aa=sp.sqrt(sp.Rational(2,3))/sp.sqrt(P)
tt=sp.exp(aa*r)
Vr=ui*sp.exp(-2*(chi-ci))/tt**2
symbolic('canonical_modulus_slope', sp.diff(Vr,r)+2*aa*Vr)
symbolic('canonical_modulus_curvature', sp.diff(Vr,r,2)-8*Vr/(3*P))
symbolic('modulus_metric_canonical', 3*P/(2*tt**2)*sp.diff(tt,r)**2-1)
G=F**2/tt
Gamma=aa*G/2
symbolic('modulus_fixed_dynamic_condition', Gamma*(4*Vr/G)+sp.diff(Vr,r))
# Fixed-r radiation tracker: dot chi=2H, dot H=-2H^2, V=F^2 H^2/t.
H=sp.symbols('H',positive=True)
symbolic('radiation_tracker_clock_equation', -4*H**2+3*H*(2*H)-2*t/F**2*(F**2*H**2/t))
symbolic('radiation_tracker_modulus_equation', aa*(F**2/t)*(2*H)**2/2-2*aa*(F**2*H**2/t))
symbolic('radiation_tracker_energy_fraction', (sp.Rational(1,2)*F**2/t*(2*H)**2+F**2*H**2/t)/(3*P*H**2)-F**2/(P*t))

# At w_clock=0 the surviving clock Weyl fermion is massless after the goldstino is eaten.
zz,zzb,Ts,Tsb=sp.symbols('zz zzb Ts Tsb')
KK=-3*P*sp.log(Ts+Tsb+(zz-zzb)**2/(6*P))
fields=(Ts,zz); bars=(Tsb,zzb); axis={Ts:t/2,Tsb:t/2,zz:0,zzb:0}
gg=sp.Matrix([[sp.diff(KK,h,bb) for bb in bars] for h in fields])
gax=gg.subs(axis)
gamTzz=sum(gax.inv().T[0,j]*sp.diff(gg[1,j],zz).subs(axis) for j in range(2))
symbolic('clock_fermion_Kzz',sp.diff(KK,zz,zz).subs(axis)+1/t)
symbolic('clock_fermion_connection_TZZ',gamTzz-1/(3*P))
symbolic('clock_fermion_covariant_mass_Wconstant',sp.diff(KK,zz,zz).subs(axis)*W/P-gamTzz*(-3*W/t))
symbolic('clock_goldstino_projection_zero',sp.diff(KK,zz).subs(axis)*W/P)
gamTTT=sum(gax.inv().T[0,j]*sp.diff(gg[0,j],Ts).subs(axis) for j in range(2))
DT=-3*W/t
unprojectedTT=sp.diff(sp.diff(KK,Ts)*W/P,Ts).subs(axis)+sp.diff(KK,Ts).subs(axis)*DT/P-gamTTT*DT
symbolic('modulus_unprojected_fermion_mass',unprojectedTT-6*W/t**2)
symbolic('modulus_goldstino_projected_mass',unprojectedTT-sp.Rational(2,3)*DT**2/W)

# Fully complex direct K/W potential, 80-digit arithmetic, three Wc and off-axis Q points.
mp.mp.dps=90
complex_tests=[]
for pn,(zr,zi,qpr,qpi,qmr,qmi,tr) in enumerate([
    (.3,.1,.07,.03,-.02,.04,1.2),(.8,-.12,.02,-.05,.03,.01,.9),(.5,.0,.0,.0,.0,.0,1.)]):
    z=mp.mpc(zr,zi); q1=mp.mpc(qpr,qpi); q2=mp.mpc(qmr,qmi)
    pp=mp.mpf('1.7'); ff=mp.mpf('.6'); uu=mp.mpf('.08')
    kval=-(z-z.conjugate())**2/2+abs(q1)**2+abs(q2)**2
    yy=mp.mpf(tr)-kval/(3*pp)
    ka=[-(z-z.conjugate()),q1.conjugate(),q2.conjugate()]
    kb=[x.conjugate() for x in ka]
    gm=mp.matrix(4); gm[0,0]=3*pp/yy**2
    for j in range(3):
        gm[0,j+1]=-kb[j]/yy**2;gm[j+1,0]=-ka[j]/yy**2
        for k in range(3):gm[j+1,k+1]=(int(j==k)/yy+ka[j]*kb[k]/(3*pp*yy**2))
    ww=-ff*mp.sqrt(uu)/mp.sqrt(2)*mp.exp(-mp.sqrt(2)*z/ff)
    wz=mp.sqrt(uu)*mp.exp(-mp.sqrt(2)*z/ff)
    mm=mp.mpf('.4')*mp.exp(mp.mpf('.12')*z); mz=mp.mpf('.12')*mm
    wa=[wz+mz*q1*q2,mm*q2,mm*q1]
    target=sum(abs(v)**2 for v in wa)/yy**2
    for wc in (mp.mpc(0),mp.mpc('.7','.3'),mp.mpc('1e12','-3e11')):
        WT=ww+wc+mm*q1*q2
        dv=mp.matrix([-3*WT/yy,*[wa[j]+ka[j]*WT/(pp*yy) for j in range(3)]])
        db=dv.conjugate()
        direct=((dv.T*(gm**-1).T*db)[0]-3*abs(WT)**2/pp)/yy**3
        err=abs(direct-target)/abs(target)
        numeric(f'direct_complex_KW_{pn}_{str(wc)}',err,mp.mpf('1e-55'))
        complex_tests.append(dict(point=pn,Wc=str(wc),relative_error=str(err)))

# Independent finite-CW evaluation, using previous chi samples solely as coordinates.
inp=REPO/'experiments/sugra_shift_v01/results/inputs.json'
traj=REPO/'experiments/sugra_shift_v01/results/trajectories.csv'
iv=json.loads(inp.read_text())
with traj.open() as f: old=[x for x in csv.DictReader(f) if x['case']=='shift_axis']
cc=np.array([float(v['chi']) for v in old]); xi=iv['xi']; f2=iv['F_squared_eV2']; p2=iv['Mpl_eV']**2
uu=iv['potential_unit_eV4']*iv['W_i']*np.exp(-2*(cc-iv['chi_i']))
mass2=.5*1e22*(1+xi/cc**2); logr=-xi/(cc*(cc**2+xi))
logrp=xi*(3*cc**2+xi)/(cc**2*(cc**2+xi)**2)
d=2*uu/(3*p2);dp=-2*d;s=mass2;ss=2*logr*s
h2=2*uu*s*logr**2/f2
h2p=h2*(-2+2*logr+2*logrp/logr)
L=np.log(s/1e22)
lead=(4*(s*dp*(L-1)+d*ss*L)+2*h2p*L+2*h2*ss/s)/(32*np.pi**2)
ratio=np.abs(lead)/(2*uu)
fixed_grid=dict(count=len(old),max_force_ratio=float(ratio.max()),index=int(ratio.argmax()),
               chi_at_max=float(cc[ratio.argmax()]),z_at_max=float(old[ratio.argmax()]['z']),
               max_d_over_s=float(np.max(d/s)),max_absB_over_s=float(np.max(np.sqrt(h2)/s)),
               meaning='Fixed-t=1 local charged Coleman-Weinberg diagnostic, no cosmic trajectory.')
mp.mp.dps=190
mpiv={k:mp.mpf(str(iv[k])) for k in ('xi','F_squared_eV2','Mpl_eV','potential_unit_eV4','W_i','chi_i')}
cw_samples=[]
for n in (0,512,1024):
    x0=mp.mpf(old[n]['chi'])
    for mt in (mp.mpf('.5'),mp.mpf('1'),mp.mpf('2')):
        def terms(x):
            xi=mpiv['xi']; U=mpiv['potential_unit_eV4']*mpiv['W_i']*mp.exp(-2*(x-mpiv['chi_i']))
            S=mp.mpf('5e21')*(1+xi/x**2)/mt
            RR=-xi/(x*(x*x+xi)); DD=2*U/(3*mpiv['Mpl_eV']**2*mt**2)
            HH=2*U*(S*mt)*RR**2/(mpiv['F_squared_eV2']*mt**2)
            return S,DD,HH,U/mt**2
        def directcw(x):
            S,DD,HH,_=terms(x);h=mp.sqrt(HH)
            def fn(v):return v*v*(mp.log(v/mp.mpf('1e22'))-mp.mpf('1.5'))
            return (fn(S+DD+h)+fn(S+DD-h)-2*fn(S))/(32*mp.pi**2)
        def leadingcw(x):
            S,DD,HH,_=terms(x);LL=mp.log(S/mp.mpf('1e22'))
            return (4*S*DD*(LL-1)+2*HH*LL)/(32*mp.pi**2)
        val=directcw(x0); leadval=leadingcw(x0)
        grad=mp.diff(directcw,x0); leadgrad=mp.diff(leadingcw,x0)
        er=abs(grad-leadgrad)/max(abs(leadgrad),mp.mpf('1e-100'))
        ev=abs(val-leadval)/max(abs(leadval),mp.mpf('1e-100'))
        numeric(f'CW_direct_derivative_{n}_t{mt}',er,mp.mpf('1e-70'))
        numeric(f'CW_direct_value_{n}_t{mt}',ev,mp.mpf('1e-70'))
        cw_samples.append(dict(index=n,t=str(mt),value_eV4=str(val),gradient_eV4=str(grad),
                               relative_value_error=str(ev),relative_gradient_error=str(er)))

out=dict(created_utc='2026-09-25',scope='Specified W_T=0 no-scale model only; no stabilization and no new cosmological trajectory.',
         python=platform.python_version(),sympy=sp.__version__,mpmath=mp.__version__,numpy=np.__version__,
         elapsed_seconds=time.time()-started,
         checks=records,checks_passed=sum(v['passed'] for v in records),check_count=len(records),
         complex_KW=complex_tests,fixed_grid_CW=fixed_grid,CW_samples=cw_samples,
         hashes={str(p.relative_to(REPO)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),inp,traj)},
         primary_sources=[
             dict(url='https://arxiv.org/abs/1503.08867',title='Phenomenological Aspects of No-Scale Inflation Models',
                  scope='Abstract and metadata only, not a full-paper read; pure no-scale soft-term examples.',http_status=200,
                  sha256='dc06a0388ebde83c8c6cb7fcfa4fbd89879b2b7e0c8a818c31736cc967abdbb5',bytes=44528),
             dict(url='https://link.springer.com/article/10.1140/epjc/s10052-017-4805-x',title='No-scale SU(5) super-GUTs',
                  scope='Publisher HTML section 2.2, equations 4-8, and adjacent introduction; not all phenomenology.',http_status=200,
                  sha256='7c6688b481fe8429e5932b16ab9a8aa5d8e0099896128572e5414a792a6b6f00',bytes=705568,
                  retrieval_note='urllib response redirected to publisher URL with cookies_not_supported query; web reader independently exposed equations 4-8. Hash fingerprints delivered HTML only.')],
         tool_errors=[{'type':'search_expression','message':'First rg used unrecognized regex escape \\c; replaced with bounded file reads. No scientific result affected.'},
                      {'type':'source_fetch','message':'web open arxiv.org/html/1503.08867 returned 406; abstract and publisher full text used.'},
                      {'type':'initial_symbolic_check','message':'First run 45/46: expected general F^T contained k_bar instead of k_hol in its matter contraction; direct inverse established k_hol. Corrected expected identity; axis k_Z=0, scalar potential, spectra and CW unaffected.'}])
(HERE/'noscale_exact_20260925_checks.json').write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
print(json.dumps({'passed':out['checks_passed'],'count':out['check_count'],'fixed_grid_CW':fixed_grid},indent=2))
if not all(v['passed'] for v in records):raise SystemExit(1)
