#!/usr/bin/env python3
"""Independent finite-K/W valley and threshold verification; never imports main code."""
from pathlib import Path
import argparse,csv,json,hashlib,platform,time
import numpy as np
from scipy.optimize import brentq
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]; REPO=ROOT.parents[1]
OUT=ROOT/'results'; OUT.mkdir(exist_ok=True)
IP=REPO/'experiments/sugra_shift_v01/results/inputs.json'
TP=REPO/'experiments/sugra_shift_v01/results/trajectories.csv'
raw=json.loads(IP.read_text())
with TP.open() as f: old=[r for r in csv.DictReader(f) if r['case']=='shift_axis']
checks=[]; phase={'zero':(1,0),'quadrature':(0,1),'pi':(-1,0)}

def m(v):return mp.mpf(str(v))
def ck(name,err,tol='1e-50'):
    checks.append(dict(name=name,residual=str(abs(err)),tolerance=tol,passed=bool(abs(err)<m(tol))))
def params():
    return dict(P=m(raw['Mpl_eV'])**2,F=mp.sqrt(m(raw['F_squared_eV2'])),
                eps=m(raw['epsilon']),Ui=m(raw['potential_unit_eV4'])*m(raw['W_i']),
                ci=m(raw['chi_i']),xi=m(raw['xi']),
                rho=3*m(raw['omega_lambda'])*m(raw['Mpl_eV'])**2*m(raw['H_ref_eV'])**2)

def geometry(delta,v,y,eta,a,b,I):
    f=1+delta+eta*delta**2+a*delta**4+b*v**4
    p=1+2*eta*delta+4*a*delta**3-4j*b*v**3
    q=2*eta+12*a*delta**2+12*b*v**2
    om=f-I['F']**2*y*y/(3*I['P'])
    kz=-1j*mp.sqrt(2)*I['F']*y
    dt=abs(p)**2-q*(om+abs(kz)**2/(3*I['P']))
    gtt=3*I['P']*(abs(p)**2-om*q)/om**2
    return f,p,q,om,kz,dt,gtt

def clock(x,y,I):
    uu=I['Ui']*mp.exp(-2*(x-I['ci']))
    ww=-I['F']*mp.sqrt(uu/2)*mp.exp(-1j*y)
    wz=mp.sqrt(uu)*mp.exp(-1j*y)
    zeta=x+1j*y
    mass_pref=m('1e22')/(1+I['xi']/I['ci']**2)
    mm=mp.sqrt(mass_pref*(1+I['xi']/zeta**2))
    mmz=-mp.sqrt(2)/I['F']*mass_pref*I['xi']/(zeta**3*mm)
    return uu,ww,wz,mm,mmz

def potential(delta,v,x,y,eta,a,b,wc,I,clock_on=True):
    f,p,q,om,kz,dt,gtt=geometry(delta,v,y,eta,a,b,I)
    uu,ww,wz,_,_=clock(x,y,I)
    if not clock_on:ww=wz=mp.mpc(0);uu=mp.mpf(0)
    S=3*(wc+ww)-kz.conjugate()*wz
    return (abs(wz)**2+q*abs(S)**2/(3*I['P']*dt))/om**2

def direct_kw(delta,v,x,y,eta,a,b,wc,I):
    f,p,q,om,kz,dt,gtt=geometry(delta,v,y,eta,a,b,I)
    uu,ww,wz,_,_=clock(x,y,I);wt=wc+ww
    g=mp.matrix([[gtt,-p*kz.conjugate()/om**2],
                 [-p.conjugate()*kz/om**2,1/om+abs(kz)**2/(3*I['P']*om**2)]])
    dw=mp.matrix([-3*p*wt/om,wz+kz*wt/(I['P']*om)])
    return mp.re(((dw.T*(g**-1).T*dw.conjugate())[0]-3*abs(wt)**2/I['P'])/om**3)

def tune(mg,a,b,mode,I):
    target=I['rho'];s0=48*a*mg**2
    vacuum_scalar=-3*s0**2/(64*mp.pi**2)
    if mode=='scalar_retuned_eta':target-=vacuum_scalar
    e0=target/(6*mg**2*I['P']);d0=e0/(6*a);wc=mg*I['P']
    def off(delta,eta):return potential(delta,0,I['ci'],0,eta,a,b,wc,I,False)
    def fn(dn,en):
        de=dn*d0;ee=en*e0
        return ((off(de,ee)-target)/target,
                mp.diff(lambda dd:off(dd,ee),de)/(mg**2*I['P']*e0))
    dn,en=mp.findroot(fn,(1,1),tol=mp.mpf(10)**(-mp.mp.dps+30),maxsteps=25)
    eta=en*e0;dc=dn*d0
    return eta,dc,target,vacuum_scalar,e0

def mass_squares(delta,x,y,eta,a,b,wc,I):
    gg=geometry(delta,0,y,eta,a,b,I)[6]/2
    mt=mp.diff(lambda dd:potential(dd,0,x,y,eta,a,b,wc,I),delta,2)/gg
    mv=mp.diff(lambda vv:potential(delta,vv,x,y,eta,a,b,wc,I),0,2)/gg
    return mt,mv

def spectrum(delta,x,y,eta,a,b,wc,I):
    f,p,q,om,kz,dt,gtt=geometry(delta,0,y,eta,a,b,I)
    uu,ww,wz,mm,mmz=clock(x,y,I)
    S=3*(wc+ww)-kz.conjugate()*wz
    V=potential(delta,0,x,y,eta,a,b,wc,I)
    B=(wz.conjugate()*mmz+q*S.conjugate()*(mm-kz.conjugate()*mmz)/(3*I['P']*dt))/om
    return abs(mm)**2/om,2*V/(3*I['P']),abs(B)**2

def leading_cw(s,d,h2):
    L=mp.log(s/m('1e22'))
    return (4*s*d*(L-1)+2*h2*L)/(32*mp.pi**2)

def one_point(mg,a,ph,n,mode,I,tuned):
    b=a; co,si=phase[ph];wc=mg*I['P']*mp.mpc(co,si);x=m(old[n]['chi'])
    eta,dc,target,vacuum_scalar,e0=tuned
    U=clock(x,0,I)[0];W=clock(x,0,I)[1]+wc
    ds=eta/(6*a)+U*I['P']/(36*a*abs(W)**2)
    fun=lambda de:potential(de,0,x,0,eta,a,b,wc,I)
    norm=U+abs(eta)*mg**2*I['P']
    dn=mp.findroot(lambda dd:mp.diff(fun,dd*ds)/norm,1,tol=mp.mpf(10)**(-mp.mp.dps+35),maxsteps=25)
    de=dn*ds; V=fun(de)
    f,p,q,om,kz,dt,gtt=geometry(de,0,0,eta,a,b,I)
    Vdd=mp.diff(fun,de,2)
    dx=-mp.diff(lambda xx:mp.diff(lambda dd:potential(dd,0,xx,0,eta,a,b,wc,I),de),x)/Vdd
    dy=-mp.diff(lambda yy:mp.diff(lambda dd:potential(dd,0,x,yy,eta,a,b,wc,I),de),0)/Vdd
    Vx=mp.diff(lambda xx:potential(de,0,xx,0,eta,a,b,wc,I),x)
    Vy=mp.diff(lambda yy:potential(de,0,x,yy,eta,a,b,wc,I),0)
    mt,mv=mass_squares(de,x,0,eta,a,b,wc,I)
    mts=[];mvs=[]
    for kind,coord,dd in [('chi',x,dx),('y',mp.mpf(0),dy)]:
        def vals(c):
            return mass_squares(de,c,0,eta,a,b,wc,I) if kind=='chi' else mass_squares(de,x,c,eta,a,b,wc,I)
        mm=[mp.diff(lambda c:vals(c)[j],coord)+dd*mp.diff(lambda de_:mass_squares(de_,x,0,eta,a,b,wc,I)[j],de) for j in (0,1)]
        mts.append(mm[0]);mvs.append(mm[1])
    s,d,h2=spectrum(de,x,0,eta,a,b,wc,I)
    derivs=[]
    for kind,coord,dd in [('chi',x,dx),('y',mp.mpf(0),dy)]:
        def vals(c):
            return spectrum(de,c,0,eta,a,b,wc,I) if kind=='chi' else spectrum(de,x,c,eta,a,b,wc,I)
        derivs.append([mp.diff(lambda c:vals(c)[j],coord)+dd*mp.diff(lambda de_:spectrum(de_,x,0,eta,a,b,wc,I)[j],de) for j in range(3)])
    def loop_grad(j):
        sp,dp,hp=derivs[j];L=mp.log(s/m('1e22'))
        common=4*(s*dp*(L-1)+d*sp*L)+4*d*dp*L+2*d*d*sp/s
        common+=2*d*d*dp/s-mp.mpf(2)/3*d**3*sp/s**2
        st=s+d;Lt=mp.log(st/m('1e22'))
        split=2*hp*Lt+2*h2*(sp+dp)/st-h2*hp/(3*st**2)+h2*h2*(sp+dp)/(3*st**3)
        return (common+split)/(32*mp.pi**2)
    mu2=48*a*mg**2
    Tgrad=[sum(ms*msd*(mp.log(ms/mu2)-1)/(32*mp.pi**2) for ms,msd in ((mt,mts[j]),(mv,mvs[j]))) for j in range(2)]
    prefix=f'{mode}_a{a}_m{mg}_{ph}_{n}'
    ck(prefix+'_vacuum_tune',(potential(dc,0,x,0,eta,a,b,mg*I['P'],I,False)-target)/target)
    ck(prefix+'_valley_gradient',mp.diff(fun,de)/norm)
    ck(prefix+'_KW_potential',(direct_kw(de,0,x,0,eta,a,b,wc,I)-V)/V)
    # Check nonzero y too: this probes the S(y) term and Kähler mixing.
    vytest=potential(de,0,x,m('.07'),eta,a,b,wc,I)
    ck(prefix+'_KW_offaxis',(direct_kw(de,0,x,m('.07'),eta,a,b,wc,I)-vytest)/vytest)
    ck(prefix+'_positive_geometry_and_masses',0 if min(f,dt,mt,mv,s+d-mp.sqrt(h2))>0 else 1,'0.5')
    rec=dict(mode=mode,a=str(a),mG_eV=str(mg),phase=ph,index=n,precision=mp.mp.dps,
             eta=str(eta),eta_leading=str(e0),delta_off=str(dc),delta=str(de),delta_leading=str(ds),
             delta_chi=str(dx),delta_y=str(dy),V_tree=str(V),V_chi=str(Vx),V_y=str(Vy),
             target_vacuum=str(target),scalar_vacuum_fixed=str(vacuum_scalar),
             f=str(f),D=str(dt),mT_squared=str(mt),mV_squared=str(mv),
             mT_chi=str(mts[0]),mT_y=str(mts[1]),mV_chi=str(mvs[0]),mV_y=str(mvs[1]),
             s=str(s),d=str(d),B_abs2=str(h2),
             s_chi=str(derivs[0][0]),d_chi=str(derivs[0][1]),B_abs2_chi=str(derivs[0][2]),
             s_y=str(derivs[1][0]),d_y=str(derivs[1][1]),B_abs2_y=str(derivs[1][2]),
             charged_CW_chi=str(loop_grad(0)),charged_CW_y=str(loop_grad(1)),
             scalar_CW_chi=str(Tgrad[0]),scalar_CW_y=str(Tgrad[1]),
             charged_force_ratio=str(mp.sqrt(loop_grad(0)**2+loop_grad(1)**2)/(2*U)),
             scalar_force_ratio=str(mp.sqrt(Tgrad[0]**2+Tgrad[1]**2)/(2*U)),
             relative_delta_asymptotic_error=str(abs(de/ds-1)))
    return rec

def save(name,rows,start):
    with (OUT/(name+'.csv')).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0],lineterminator='\n');w.writeheader();w.writerows(rows)
    data=dict(created_utc='2026-09-25',scope='Independent full K/W valley and local threshold audit; no ODE or observations.',
              python=platform.python_version(),mpmath=mp.__version__,rows=len(rows),checks=checks,
              passed=sum(c['passed'] for c in checks),count=len(checks),elapsed_seconds=time.time()-start,
              input_hashes={str(p.relative_to(REPO)):hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in (IP,TP,ROOT/'protocol.md',Path(__file__))},
              corrections=[dict(stage='high_precision_mass_normalization',
                  description='Initial independent run used 0.5*mi^2 prefactor, which misses M(chi_i)=100GeV by about 7e-18 after decimal promotion of rounded xi and chi_i. Normalized by 1+xi/chi_i^2 as mandated by the fixed initial mass; all point results rerun. Initial attempt logs retained.'),
                  dict(stage='stable_CW_expansion',
                  description='The extra retuned high-mass initial matching points resolved d^2 terms: first leading-only audit 523/526, maximal local relative discrepancy 4.83e-9. Added stable d^3 and B^4 terms, keeping the original strict 1e-40 comparison threshold; initial audit JSON retained.')])
    (OUT/(name+'.json')).write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:data[k] for k in ('rows','passed','count','elapsed_seconds')}),flush=True)
    if data['passed']!=data['count']:raise SystemExit(1)

def vector_rebuild(mg,a,ph,muf,mode):
    """Independent NumPy form of the protocol's retained adiabatic order."""
    P=float(raw['Mpl_eV'])**2;F2=float(raw['F_squared_eV2']);F=np.sqrt(F2)
    ci=float(raw['chi_i']);xi=float(raw['xi'])
    x=np.array([float(r['chi']) for r in old]);co,si=phase[ph]
    U=float(raw['potential_unit_eV4'])*float(raw['W_i'])*np.exp(-2*(x-ci))
    w=-np.sqrt(F2*U/2);rw=w/(mg*P)
    wn=1+2*rw*co+rw*rw;W=mg*P*(co+1j*si)+w
    rho=3*float(raw['omega_lambda'])*P*float(raw['H_ref_eV'])**2
    mt0=48*a*mg**2;L0=-2*np.log(muf)
    vac=mt0*mt0*(L0-1.5)/(32*np.pi**2)
    eta=(rho-(vac if mode=='scalar_retuned_eta' else 0))/(6*P*mg**2)
    delta=eta/(6*a)+U/(36*a*mg**2*P*wn)
    small=delta+eta*delta**2+a*delta**4
    f=1+small;p=1+2*eta*delta+4*a*delta**3;q=2*eta+12*a*delta**2
    D=p*p-f*q
    V=U/f**2+3*q*mg**2*P*wn/(f*f*D)
    ex=-2*U*np.expm1(-2*np.log1p(small))-6*q*mg**2*P*(rw*co+rw*rw)/(f*f*D)
    ey=-2*q*mg**2*P*rw*si/(f*f*D)
    Vx=-2*U+ex;Vy=ey
    ss=1e22*(1+xi/x**2)/(1+xi/ci**2)
    M=np.sqrt(ss);rr=-xi/(x*(x*x+xi))
    rrp=xi*(3*x*x+xi)/(x*x*(x*x+xi)**2)
    Mx=M*rr;Mxx=M*(rr*rr+rrp)
    Mz=np.sqrt(2)/F*Mx;Mzx=np.sqrt(2)/F*Mxx
    B=(np.sqrt(U)*Mz+q*np.conj(W)*M/(P*D))/f
    Bx=(np.sqrt(U)*(Mzx-Mz)+q*(np.conj(W)*Mx-w*M)/(P*D))/f
    By=1j*(np.sqrt(U)*(Mzx+Mz)+q*(w*M/3-np.conj(W)*Mx)/(P*D))/f
    s=ss/f;sx=2*rr*s
    d=2*V/(3*P);dx=2*Vx/(3*P);dy=2*Vy/(3*P)
    h2=np.abs(B)**2;h2x=2*np.real(np.conj(B)*Bx);h2y=2*np.real(np.conj(B)*By)
    L=np.log(s/1e22)
    Qx=(4*(s*dx*(L-1)+d*sx*L)+2*h2x*L+2*h2*sx/s)/(32*np.pi**2)
    Qy=(4*s*dy*(L-1)+2*h2y*L)/(32*np.pi**2)
    mt=mt0*wn;mtx=-2*mt0*(rw*co+rw*rw);mty=-(2.0/3.0)*mt0*rw*si
    LT=np.log1p(2*rw*co+rw*rw)+L0
    Tx=mt*mtx*(LT-1)/(16*np.pi**2);Ty=mt*mty*(LT-1)/(16*np.pi**2)
    tree=np.hypot(ex,ey)/(2*U);charged=np.hypot(Qx,Qy)/(2*U);mod=np.hypot(Tx,Ty)/(2*U)
    selected=np.hypot(ex+Qx+Tx,ey+Qy+Ty)/(2*U);components=tree+charged+mod
    hh=np.array([float(r['E']) for r in old])*float(raw['H_ref_eV'])
    qx=np.array([float(r['qx']) for r in old])
    dkin=-F2*(hh*qx)**2/(144*a*mg**2*P*wn)
    marker=float(raw['Mpl_eV'])/np.sqrt(a);fcomp=np.abs(q*np.conj(W)/(P*np.sqrt(f)*D))
    g2=4*np.pi*float(raw['alpha_reference'])
    summary=dict(eta=eta,scalar_CW_clockoff_eV4=vac,max_abs_delta=np.max(abs(delta)),
        max_abs_delta_kinetic=np.max(abs(dkin)),min_f=np.min(f),min_D=np.min(D),
        max_tree_force_ratio=np.max(tree),max_charged_force_ratio=np.max(charged),
        max_modulus_force_ratio=np.max(mod),max_selected_force_ratio=np.max(selected),
        max_component_sum_ratio=np.max(components),
        max_selected_x_ratio=np.max(abs(ex+Qx+Tx)/(2*U)),max_selected_y_ratio=np.max(abs(ey+Qy+Ty)/(2*U)),
        min_modulus_mass_eV=np.min(np.sqrt(mt)),max_H_over_modulus_mass=np.max(hh/np.sqrt(mt)),
        max_mass_over_geometric_marker=max(np.max(np.sqrt(s)),np.max(np.sqrt(mt)))/marker,
        geometric_marker_eV=marker,max_d_over_s=np.max(abs(d/s)),max_B_over_s=np.max(np.sqrt(h2)/s),
        min_scalar_over_fermion_mass2=np.min(1+d/s-np.sqrt(h2)/s),
        max_Fcomp_abs_eV=np.max(fcomp),max_above_threshold_AMSB_abs_eV=np.max(g2*2*fcomp/(16*np.pi**2)))
    return summary

def audit():
    start=time.time();mp.mp.dps=240
    def read(n):
        with (OUT/n).open() as f:return list(csv.DictReader(f))
    t220=read('independent_tree_eta_220.csv');r220=read('independent_scalar_retuned_eta_220.csv')
    t160=sum([read(f'independent_tree_eta_160_a{a}.csv') for a in ('0.1','1','10')],[])
    r160=read('independent_scalar_retuned_eta_160.csv')
    key=lambda r:(r['mode'],float(r['a']),float(r['mG_eV']),r['phase'],int(r['index']))
    low={key(r):r for r in t160+r160};precision_errors={}
    fields=['eta','delta','V_tree','V_chi','V_y','mT_squared','mV_squared','mT_chi','mT_y','mV_chi','mV_y',
            's','d','B_abs2','charged_CW_chi','charged_CW_y','scalar_CW_chi','scalar_CW_y']
    for field in fields:
        err=max(abs(m(r[field])-m(low[key(r)][field]))/max(abs(m(r[field])),m('1e-200')) for r in t220+r220)
        precision_errors[field]=str(err);ck('precision_160_220_'+field,err,'1e-65')
    # Direct, unexpanded charged supertrace at 180 and 240 digits for both eta choices.
    direct=[]
    for r in [r for r in t220+r220 if float(r['a'])==1]:
        vals={}
        for prec in (180,240):
            with mp.workdps(prec):
                s,d,h2=(m(r[k]) for k in ('s','d','B_abs2'));h=mp.sqrt(h2)
                def f(v):return v*v*(mp.log(v/m('1e22'))-m('1.5'))
                def fp(v):return 2*v*(mp.log(v/m('1e22'))-1)
                V=(f(s+d+h)+f(s+d-h)-2*f(s))/(32*mp.pi**2)
                gg=[]
                for direction in ('chi','y'):
                    sp,dp,hp=(m(r[k+'_'+direction]) for k in ('s','d','B_abs2'))
                    dh=hp/(2*h)
                    gg.append((fp(s+d+h)*(sp+dp+dh)+fp(s+d-h)*(sp+dp-dh)-2*fp(s)*sp)/(32*mp.pi**2))
                vals[prec]=(V,*gg)
        norm=max(abs(x) for x in vals[240])
        er=max(abs(x-y) for x,y in zip(vals[180],vals[240]))/norm
        ck('direct_CW_precision_'+str(key(r)),er,'1e-65')
        # The leading split expansion omits terms of order tiny d/s and h^2/s^2.
        eg=max(abs(vals[240][j+1]-m(r['charged_CW_'+dr])) for j,dr in enumerate(('chi','y')))/max(abs(vals[240][1]),abs(vals[240][2]),m('1e-200'))
        ck('direct_CW_vs_expansion_'+str(key(r)),eg,'1e-40')
        direct.append(dict(mode=r['mode'],a=r['a'],mG_eV=r['mG_eV'],phase=r['phase'],index=r['index'],
                           CW_value=str(vals[240][0]),CW_chi=str(vals[240][1]),CW_y=str(vals[240][2]),
                           precision_error=str(er),expansion_error=str(eg)))
    with (OUT/'independent_direct_CW.csv').open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=direct[0],lineterminator='\n');wr.writeheader();wr.writerows(direct)
    # Compare exact local results to production after normalizing each vector by its nonzero norm.
    reps={ (r['phase'],float(r['mG_eV']),int(r['index'])):r for r in read('representative_curves.csv') }
    comparison=[]
    for r in [r for r in t220+r220 if float(r['a'])==1]:
        ref=reps[(r['phase'],float(r['mG_eV']),int(r['index']))]
        pre='tree_eta_' if r['mode']=='tree_eta' else 'retuned_'
        for label,own,names in [
            ('charged',('charged_CW_chi','charged_CW_y'),('charged_CW_x','charged_CW_y')),
            ('modulus',('scalar_CW_chi','scalar_CW_y'),('modulus_CW_x','modulus_CW_y')),
            ('tree',('V_chi','V_y'),('Vtree_x','Vtree_y'))]:
            ev=[m(r[k]) for k in own];pv=[m(ref[pre+k]) for k in names]
            # Pointwise tiny matching-point components can depend on retained O(delta) terms;
            # normalize by the corresponding whole three-chi vector group.
            group=[rr for rr in t220+r220 if float(rr['a'])==1 and rr['mode']==r['mode']
                   and rr['phase']==r['phase'] and rr['mG_eV']==r['mG_eV']]
            scale=max(abs(m(rr[k])) for rr in group for k in own)
            er=max(abs(ev[j]-pv[j]) for j in (0,1))/max(scale,m('1e-200'))
            ck('exact_vs_main_'+label+'_'+str(key(r)),er,'1e-8')
            comparison.append(dict(mode=r['mode'],mG_eV=r['mG_eV'],phase=r['phase'],index=r['index'],
                                   component=label,group_normalized_error=str(er),
                                   exact_chi=str(ev[0]),main_chi=str(pv[0]),exact_y=str(ev[1]),main_y=str(pv[1])))
        for field in ('s','d','B_abs2'):
            ck('exact_spectrum_vs_main_'+field+'_'+str(key(r)),
               abs(m(r[field])-m(ref[pre+field]))/max(abs(m(r[field])),m('1e-200')),'1e-8')
    with (OUT/'independent_exact_main_comparison.csv').open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=comparison[0],lineterminator='\n');wr.writeheader();wr.writerows(comparison)
    scan_errors={};scan_rows=read('scheme_scan.csv')
    primary_scan=read('scan.csv')
    scheme_index={(r['mG_eV'],r['a'],r['phase'],r['mu_factor'],r['scheme']):r for r in scan_rows}
    embedded=all(r==scheme_index.get((r['mG_eV'],r['a'],r['phase'],r['mu_factor'],r['scheme'])) for r in primary_scan)
    ck('primary_scan_873_rows_is_exact_subset',0 if len(primary_scan)==873 and embedded else 1,'0.5')
    for j,r in enumerate(scan_rows):
        out=vector_rebuild(float(r['mG_eV']),float(r['a']),r['phase'],float(r['mu_factor']),r['scheme'])
        for k,v in out.items():
            ref=float(r[k]);er=abs(v-ref)/max(abs(v),abs(ref),1e-300)
            if k not in scan_errors or er>scan_errors[k]['error']:
                scan_errors[k]=dict(error=float(er),row=j,main=ref,independent=float(v))
    for k,v in scan_errors.items():ck('full_scheme_scan_'+k,v['error'],'1e-8')
    bound_results=[]
    for r in read('budget_limits.csv'):
        args=(float(r['a']),r['phase'],float(r['mu_factor']),r['scheme'])
        fn=lambda z:np.log10(vector_rebuild(10**z,*args)['max_component_sum_ratio'])+1
        root=10**brentq(fn,0,12,xtol=1e-12)
        er=abs(root/float(r['component_sum_10percent_mG_eV'])-1)
        ck('budget_root_'+str(args),er,'1e-8')
        bound_results.append(dict(**r,independent_bound=root,relative_error=er))
    summary=dict(created_utc='2026-09-25',scope='Independent high-precision and full-summary audit; not a full SUGRA loop.',
        elapsed_seconds=time.time()-start,precision_point_pairs=108,tree_points=81,retuned_points=27,
        direct_CW_points=len(direct),primary_scan_rows=len(primary_scan),scheme_rows=len(scan_rows),reconstructed_background_points=len(scan_rows)*len(old),
        limits=len(bound_results),precision_errors=precision_errors,full_scan_errors=scan_errors,bounds=bound_results,
        checks=checks,passed=sum(c['passed'] for c in checks),count=len(checks),
        max_main_vector_error=max(float(r['group_normalized_error']) for r in comparison),
        input_hashes={str(p.relative_to(REPO)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [Path(__file__),ROOT/'protocol.md',IP,TP,OUT/'representative_curves.csv',OUT/'scan.csv',OUT/'scheme_scan.csv',OUT/'budget_limits.csv']},
        corrections=[
            'Initial decimal promotion violated the fixed M(chi_i)=100 GeV normalization by roughly 7e-18; normalization enforced and all exact samples rerun.',
            'Initial leading-only split audit 523/526: three additional retuned matching-point relative residuals reached 4.83e-9. Stable d^3/B^4 series then used without weakening the 1e-40 direct-CW threshold; initial JSON preserved.'],
        interpretation_notes=[
            'Quadrature subleading chi modulus-mass derivatives include U/Mpl^2 terms; leading y force is independently recovered. Comparisons use the predeclared nonzero vector/group norm.',
            'The moduli masses here are canonical local T-direction potential curvatures at a clock-conditioned valley. They do not establish a full cosmological perturbation spectrum or global stability.',
            'The vector audit independently implements the declared adiabatic expansion; exact high-precision points test its domain. No primary code was read or imported.'],
        execution_errors=['One audit launch had a mismatched dictionary-comprehension bracket after adding metadata. Fixed before execution; SyntaxError log retained.'])
    (OUT/'independent_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({k:summary[k] for k in ('passed','count','scheme_rows','direct_CW_points','elapsed_seconds','max_main_vector_error')}),flush=True)
    if summary['passed']!=summary['count']:raise SystemExit(1)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--mode',choices=['tree_eta','scalar_retuned_eta'])
    ap.add_argument('--dps',type=int,default=160);ap.add_argument('--a',default=None)
    ap.add_argument('--audit',action='store_true')
    ar=ap.parse_args()
    if ar.audit:return audit()
    if not ar.mode:ap.error('--mode is required unless --audit')
    mp.mp.dps=ar.dps;I=params();rows=[];start=time.time()
    avec=[m(ar.a)] if ar.a else ([m('.1'),m('1'),m('10')] if ar.mode=='tree_eta' else [m('1')])
    for a in avec:
        for mg in (m('1e-6'),m('1'),m('1e7')):
            tuned=tune(mg,a,a,ar.mode,I)
            for ph in phase:
                for n in (0,512,1024):
                    rows.append(one_point(mg,a,ph,n,ar.mode,I,tuned))
            print(f'done {ar.mode} a={a} mG={mg} precision={mp.mp.dps}',flush=True)
    suffix=f'_a{ar.a}' if ar.a else ''
    save(f'independent_{ar.mode}_{ar.dps}{suffix}',rows,start)

if __name__=='__main__':main()
