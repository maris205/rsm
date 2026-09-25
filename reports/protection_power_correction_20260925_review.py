#!/usr/bin/env python3
"""Adversarial audit: independent full auxiliary fields and h_Z correction to B_y."""
from pathlib import Path
import csv
import hashlib
import json
import mpmath as mp
import numpy as np

mp.mp.dps=100
ROOT=Path(__file__).resolve().parents[1]
checks=[]
def ck(name,actual,expected,scale=None):
    scale=max(abs(expected),mp.mpf('1e-30')) if scale is None else scale
    error=abs(actual-expected)/scale
    checks.append(dict(name=name,error=mp.nstr(error,12),tolerance='1e-70',passed=bool(error<mp.mpf('1e-70'))))
def jets(t,eta,A):
    return t+eta*(t-1)**2+(t-1)**4-A/t**2,1+2*eta*(t-1)+4*(t-1)**3+2*A/t**3,2*eta+12*(t-1)**2-6*A/t**4

def Bfull(t,y,eta,A,hA,Wc,w,F,chi=mp.mpf(140),xi=mp.mpf(19000)):
    """P=1, h_Q=1 control isolates the new clock metric's contribution."""
    g,p,q=jets(t,eta,A)
    h=1+hA/t**3;hp=-3*hA/t**4;hpp=12*hA/t**5
    Om=g-h*F**2*y*y/3;Ot=p-hp*F**2*y*y/3;Oy=-2*h*F**2*y/3
    Ott=q-hpp*F**2*y*y/3;Oyy=-2*h*F**2/3;Oty=-2*hp*F**2*y/3
    Kt=-3*Ot/Om;Ky=-3*Oy/Om;root=mp.sqrt(2)*F
    mat=mp.matrix([[-3*(Ott/Om-Ot*Ot/Om**2),-3j*(Oty/Om-Ot*Oy/Om**2)/root],
                   [3j*(Oty/Om-Ot*Oy/Om**2)/root,-3*(Oyy/Om-Oy*Oy/Om**2)/root**2]])
    W=Wc+w*mp.exp(-1j*y);wz=-mp.sqrt(2)*w/F*mp.exp(-1j*y)
    DW=mp.matrix([Kt*W,wz-1j*Ky/root*W])
    aux=-Om**mp.mpf('-1.5')*(mat**-1).T*DW.conjugate()
    Oz=-1j*Oy/root
    r=-xi/((chi+1j*y)*((chi+1j*y)**2+xi))
    M=mp.sqrt(1+xi/(chi+1j*y)**2)
    return M/mp.sqrt(Om)*(-mp.conj(W)/Om**mp.mpf('1.5')+aux[0]*Ot/Om+aux[1]*(Oz/Om-mp.sqrt(2)*r/F))

for A in (mp.mpf('1e-6'),mp.mpf('1e-10')):
    t=mp.findroot(lambda xx:(xx-1)*xx**5+A,1,tol=mp.mpf('1e-95'))
    eta=-6*(t-1)**2+3*A/t**4
    g,p,q=jets(t,eta,A)
    F0=mp.mpf('.07');w=mp.mpf('-1e-5');chi=mp.mpf(140);xi=mp.mpf(19000)
    for kap in (-1,0,1):
        h=1+kap*A/t**3;hp=-3*kap*A/t**4
        Fhol=F0/mp.sqrt(h) # Both controls use identical g; fixes the physical metric.
        for phase in (0,mp.pi/2):
            Wc=mp.mpf('.5')*mp.exp(1j*phase)
            f=lambda yy:Bfull(t,yy,eta,A,kap*A,Wc,w,Fhol)
            base=lambda yy:Bfull(t,yy,eta,A,0,Wc,w,F0)
            actual=mp.diff(lambda yy:f(yy)-base(yy),0)
            r=-xi/(chi*(chi**2+xi));M=mp.sqrt(1+xi/chi**2)/mp.sqrt(g)
            FT=mp.conj(Wc+w)/(p*mp.sqrt(g))
            expected=2j*M*hp/h*(r*FT+w/(3*p*mp.sqrt(g)))
            tag=f'A{A}_kap{kap}_phase{mp.nstr(phase,5)}'
            ck(tag+'_B_axis_normalization',f(0),base(0))
            ck(tag+'_B_y_clock_metric',actual,expected)

# Numerical size of the selected extra two-loop interference; old chi coordinates only.
parpath=ROOT/'experiments/sugra_shift_v01/results/inputs.json'
trajpath=ROOT/'experiments/sugra_shift_v01/results/trajectories.csv'
par=json.loads(parpath.read_text())
chi=np.array([float(row['chi']) for row in csv.DictReader(trajpath.open()) if row['case']=='shift_axis'])
P=par['Mpl_eV']**2;F=np.sqrt(par['F_squared_eV2']);rho=3*par['omega_lambda']*P*par['H_ref_eV']**2
U=par['potential_unit_eV4']*par['W_i']*np.exp(-2*(chi-par['chi_i']));xi=par['xi']
r=-xi/(chi*(chi**2+xi));rel=xi*(chi**-2-par['chi_i']**-2)/(1+xi/par['chi_i']**2)
s=1e22*(1+rel);M=np.sqrt(s);L=np.log1p(rel);C=float(mp.zeta(3)/(48*mp.pi**4))
kk=1e14;A=C*kk**2/P;c=12*A
rows=[]
for mg in np.logspace(-24,2,105):
    rw=F*np.sqrt(U/2)/(mg*P);ft=mg*(-rw-1j) # Quadrature gives maximal leading phase projection.
    B0=M*(np.sqrt(2*U)/F*r+rho/(3*mg**2*P)*ft)
    deltaBy=-6j*A*M*(r*ft-mg*rw/3)
    yforce=L*2*np.real(np.conj(B0)*deltaBy)/(16*np.pi**2)
    R=float(np.max(np.abs(yforce)/(2*U)))
    X=mg**2*(rw**2+1);dX=np.array([-2*mg**2*rw**2,2*mg**2*rw/3]);ds=np.array([2*r*s,np.zeros_like(s)])
    kc=-(X*L*ds+s*(L-1)*dX)/(8*np.pi**2)
    cR=float(np.max(np.hypot(*(c*kc))/(2*U)))
    rows.append(dict(mG_eV=float(mg),mKK_eV=kk,additional_hZ_B_interference_ratio=R,ratio_to_selected_c=cR and R/cR))
result=dict(scope='Independent auxiliary-field audit; isolated h_Z B_y term, not complete two-loop matching',
            check_count=len(checks),passed=sum(c['passed'] for c in checks),all_pass=all(c['passed'] for c in checks),checks=checks,
            max_additional_hZ_interference_ratio=max(r['additional_hZ_B_interference_ratio'] for r in rows),
            max_ratio_to_selected_c=max(r['ratio_to_selected_c'] for r in rows),scan=rows,
            source_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),parpath,trajpath)},
            known_omissions='General h_Q+h_Z cross terms beyond first order A, finite-rho corrections to this small extra jet, full mixed two-loop matching and readout evolution.')
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ('check_count','passed','all_pass','max_additional_hZ_interference_ratio','max_ratio_to_selected_c')},indent=2))
if not result['all_pass']:
    raise SystemExit(1)
