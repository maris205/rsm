#!/usr/bin/env python3
"""Adversarial physical-unit fixed-T valley check, independent of the main run.

Rebuild the complete K/W matrix at 140 and 220 digits; no experiment imports.
"""
from pathlib import Path
import csv
import hashlib
import json
import mpmath as mp

HERE=Path(__file__).resolve()
ROOT=HERE.parents[1]
PAR=ROOT/'experiments/sugra_shift_v01/results/inputs.json'
TRAJ=ROOT/'experiments/sugra_shift_v01/results/trajectories.csv'
raw=json.loads(PAR.read_text(),parse_float=str)
coords=[row['chi'] for row in csv.DictReader(TRAJ.open()) if row['case']=='shift_axis']
checks=[]


def ck(name,actual,expected,scale=None,tol='1e-70'):
    den=max(abs(expected),mp.mpf('1e-130')) if scale is None else mp.mpf(scale)
    error=abs(actual-expected)/den
    checks.append(dict(name=name,error=mp.nstr(error,18),tolerance=tol,passed=bool(error<mp.mpf(tol))))


def compute(precision):
    mp.mp.dps=precision
    P=mp.mpf(raw['Mpl_eV'])**2
    F2=mp.mpf(raw['F_squared_eV2']);F=mp.sqrt(F2)
    H=mp.mpf(raw['H_ref_eV'])
    rho=3*mp.mpf(raw['omega_lambda'])*P*H*H
    mg=mp.mpf('1e-6');KK=mp.mpf('1e14');Lambda=mp.mpf('1e12')
    A=mp.zeta(3)/(48*mp.pi**4)*KK*KK/P
    rr=rho/(3*P*mg*mg)
    def geo(delta,eta):
        t=1+delta
        return (t+eta*delta**2+delta**4-A/t**2,
                1+2*eta*delta+4*delta**3+2*A/t**3,
                2*eta+12*delta**2-6*A/t**4,
                24*delta+24*A/t**5)
    def eq(delta,eta):
        g,p,q,qp=geo(delta,eta)
        return ((q-rr*p*p/(g*(1+rr)))/A,
                (qp-2*p*q/g+q*q/p)/A)
    delta,eta=mp.findroot(eq,(-A,3*A+rr/2),tol=mp.power(10,-precision+25))
    g,p,q,qp=geo(delta,eta);t=1+delta
    h=1+A/t**3;hp=-3*A/t**4;hpp=12*A/t**5
    Fhol=F*mp.sqrt(g/h);Lhol=Lambda*mp.sqrt(g/h)
    Wc=1j*P*g**mp.mpf('1.5')*mg
    D=p*p-g*q
    ck(f'{precision}_vacuum_energy',3*q*abs(Wc)**2/(P*g*g*D),rho,tol='1e-65')
    ck(f'{precision}_vacuum_stationary',qp-2*p*q/g+q*q/p,0,scale=A,tol='1e-90')
    rows=[]
    for which,chi in (('initial',mp.mpf(coords[0])),('final',mp.mpf(coords[-1]))):
        U=mp.mpf(raw['potential_unit_eV4'])*mp.mpf(raw['W_i'])*mp.exp(-2*(chi-mp.mpf(raw['chi_i'])))
        w=-g**mp.mpf('1.5')*F*mp.sqrt(U/2)
        def fields(y,x=mp.mpf(0)):
            k=Fhol**2*y*y-Fhol**4*y**4/Lhol**2
            ky=2*Fhol**2*y-4*Fhol**4*y**3/Lhol**2
            kyy=2*Fhol**2-12*Fhol**4*y*y/Lhol**2
            Om=g-h*k/(3*P);Ot=p-hp*k/(3*P);Ott=q-hpp*k/(3*P)
            Oy=-h*ky/(3*P);Oyy=-h*kyy/(3*P);Oty=-hp*ky/(3*P)
            Kt=-3*P*Ot/Om;Ky=-3*P*Oy/Om
            Ktt=-3*P*(Ott/Om-Ot*Ot/Om**2)
            Kyy=-3*P*(Oyy/Om-Oy*Oy/Om**2)
            Kty=-3*P*(Oty/Om-Ot*Oy/Om**2)
            root=mp.sqrt(2)*Fhol
            mat=mp.matrix([[Ktt,1j*Kty/root],[-1j*Kty/root,Kyy/root**2]])
            wr=w*mp.exp(-x);W=Wc+wr*mp.exp(-1j*y)
            Wz=-mp.sqrt(2)*wr/Fhol*mp.exp(-1j*y)
            dw=mp.matrix([Kt*W/P,Wz-1j*Ky*W/(root*P)])
            aux=-Om**mp.mpf('-1.5')*(mat**-1).T*dw.conjugate()
            V=mp.re(Om**-3*((dw.conjugate().T*(mat**-1)*dw)[0]-3*abs(W)**2/P))
            return dict(V=V,Omega=Om,KTT=Ktt,det=mp.re(mp.det(mat)),FT=aux[0],FZ=aux[1],W=W)
        V=lambda yy,xx=mp.mpf(0):fields(yy,xx)['V']
        J=mp.diff(V,0);Hyy=mp.diff(V,0,2)
        yq=-J/Hyy
        # Solve for an order-one scale factor so the tiny y coordinate does
        # not hide convergence; residual is normalized by the actual force.
        scale=mp.findroot(lambda z:mp.diff(V,yq*z)/J,(mp.mpf('.99'),mp.mpf('1.01')),
                          tol=mp.power(10,-precision+35))
        ys=yq*scale
        result=fields(ys)
        stationarity=mp.diff(V,ys)/J
        mass=Hyy/F2
        leadingJ=-6*mp.sqrt(2)*A*F*mp.sqrt(U)*mg
        wm=F*mp.sqrt(U/2);rw=wm/(mg*P)
        baseline=4*(1+rr)*U/(3*P)+2*rr*mg*mg*(2+rr)*(1+rw*rw)
        approxmass=baseline-24*A*mg*mg+12*U/Lambda**2
        approxmass_x=-8*(1+rr)*U/(3*P)-4*rr*mg*mg*(2+rr)*rw*rw-24*U/Lambda**2
        approx_y=-leadingJ/(F2*approxmass)
        approx_dv=-leadingJ**2/(2*F2*approxmass)
        approx_dvx=approx_dv*(-2-approxmass_x/approxmass)
        dv=V(ys)-V(0)
        dvx=mp.diff(lambda xx:V(ys,xx)-V(0,xx),0)
        exact_quadratic=-J*J/(2*Hyy)
        label=f'{precision}_{which}'
        ck(label+'_stationarity',stationarity,0,scale=1,tol='1e-65')
        ck(label+'_positive_metric',int(result['Omega']>0 and result['KTT']>0 and result['det']>0),1)
        ck(label+'_positive_valley_curvature',int(mp.diff(V,ys,2)>0),1)
        ck(label+'_quadratic_location',ys,yq,tol='1e-30')
        ck(label+'_quadratic_energy',dv,exact_quadratic,tol='1e-30')
        ck(label+'_main_leading_J',J,leadingJ,tol='1e-20')
        ck(label+'_main_mass',mass,approxmass,tol='1e-20')
        ck(label+'_main_valley_location',ys,approx_y,tol='1e-20')
        ck(label+'_main_valley_energy',dv,approx_dv,tol='1e-20')
        ck(label+'_main_valley_force',dvx,approx_dvx,tol='1e-15')
        # Full auxiliary-field B at Q=y=0 versus h_Q derivative disabled,
        # with identical physical charged mass. This isolates the beta sign.
        axis=fields(0);ellp=hp/h;Mphys=mp.mpf('1e11')
        xi=mp.mpf(raw['xi']);r=-xi/(chi*(chi*chi+xi));Mz_over_M=mp.sqrt(2)*r/Fhol
        common=-mp.conj(axis['W'])/(P*g**mp.mpf('1.5'))+axis['FT']*p/g-axis['FZ']*Mz_over_M
        Bactual=Mphys*(common-2*ellp*axis['FT']);Bcontrol=Mphys*common
        dB=Bactual-Bcontrol;beta=2*ellp
        ck(label+'_B5d_signed_auxiliary',dB,-beta*Mphys*axis['FT'])
        ck(label+'_B5d_opposite_sign_rejected',int(abs(dB-beta*Mphys*axis['FT'])>abs(dB)),1)
        row=dict(precision=precision,endpoint=which,chi=mp.nstr(chi,40),U=mp.nstr(U,50),
                 J=mp.nstr(J,65),J_leading=mp.nstr(leadingJ,65),
                 physical_mass_squared=mp.nstr(mass,65),main_mass_squared=mp.nstr(approxmass,65),
                 y_stationary=mp.nstr(ys,65),y_quadratic=mp.nstr(yq,65),y_main=mp.nstr(approx_y,65),
                 stationary_residual_over_J=mp.nstr(stationarity,30),
                 delta_V=mp.nstr(dv,65),delta_V_quadratic=mp.nstr(exact_quadratic,65),
                 delta_V_main=mp.nstr(approx_dv,65),delta_V_chi=mp.nstr(dvx,65),
                 delta_V_chi_main=mp.nstr(approx_dvx,65),
                 metric_boundary_fraction=mp.nstr(mp.sqrt(6)*F*abs(ys)/Lambda,45),
                 Omega=mp.nstr(result['Omega'],65),metric_TT=mp.nstr(result['KTT'],45),
                 metric_determinant=mp.nstr(result['det'],45),
                 delta_B5d_real=mp.nstr(mp.re(dB),45),delta_B5d_imag=mp.nstr(mp.im(dB),45))
        rows.append(row)
    return rows


rows=[]
for digits in (140,220):
    rows.extend(compute(digits))
for endpoint in ('initial','final'):
    lo=next(r for r in rows if r['precision']==140 and r['endpoint']==endpoint)
    hi=next(r for r in rows if r['precision']==220 and r['endpoint']==endpoint)
    for key in ('J','physical_mass_squared','y_stationary','delta_V','delta_V_chi'):
        ck(f'precision_{endpoint}_{key}',mp.mpf(lo[key]),mp.mpf(hi[key]),tol='1e-40')
result=dict(scope='Fixed-T complete local K/W potential; physical units, finite-rho reference, no cosmic or full heavy-spectrum solution',
            check_count=len(checks),passed=sum(c['passed'] for c in checks),all_pass=all(c['passed'] for c in checks),
            checks=checks,rows=rows,code_sha256=hashlib.sha256(HERE.read_bytes()).hexdigest(),
            source_hashes={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in (PAR,TRAJ)},
            reviewed_main='experiments/clock_stabilization_v01/code/run_stabilization.py',
            execution_notes='First physical-unit adversarial run; existing 843-check script was not rerun or altered.')
HERE.with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ('check_count','passed','all_pass')},indent=2))
if not result['all_pass']:
    print(json.dumps([c for c in checks if not c['passed']],indent=2))
    raise SystemExit(1)
