#!/usr/bin/env python3
"""Bounded algebraic audit of main interference jets; no production imports."""
from pathlib import Path
import hashlib,json
import sympy as s
root=Path(__file__).resolve().parents[1]
x,y,u,v,P,F,H,M,r,rp,eta,beta=s.symbols('x y u v P F H M r rp eta beta',real=True)
wr=-F*H/s.sqrt(2)
zeta=x+s.I*y
w=wr*s.exp(-zeta)
wz=H*s.exp(-zeta)
Wh=P*(u+s.I*v)+w
S=3*Wh-s.I*s.sqrt(2)*F*y*wz
FT=s.conjugate(S)/(3*P)
Mh=M*(1+r*zeta+(r*r+rp)*zeta*zeta/2)
Mz=M*(r+(r*r+rp)*zeta)
B0=s.conjugate(wz)*s.sqrt(2)/F*Mz+2*eta*FT*(Mh-2*s.I*y*Mz)
Delta=beta*Mh*FT
at0=lambda e:s.simplify(s.expand_complex(e.subs({x:0,y:0})))
ft=u-s.I*v+wr/P
wm=-wr
hc=s.sqrt(2)*H/F
expected=[('FT_center',at0(FT)-ft),('FT_x',at0(s.diff(FT,x))-wm/P),('FT_y',at0(s.diff(FT,y))+s.I*wm/(3*P)),('X_chi',at0(s.diff(s.conjugate(FT)*FT,x))-2*wm/P*(u-wm/P)),('X_y',at0(s.diff(s.conjugate(FT)*FT,y))-2*wm*v/(3*P)),('B0_center',at0(B0)-M*(hc*r+2*eta*ft)),('B0_chi',at0(s.diff(B0,x))-M*(hc*(r*r+rp-r)+2*eta*(r*ft+wm/P))),('B0_y',at0(s.diff(B0,y))-s.I*M*(hc*(r*r+rp+r)+2*eta*(-wm/(3*P)-r*ft))),('DeltaB_chi',at0(s.diff(Delta,x))-beta*M*(r*ft+wm/P)),('DeltaB_y',at0(s.diff(Delta,y))-s.I*beta*M*(r*ft-wm/(3*P)))]
checks=[]
for name,e in expected:
 residual=s.simplify(e);checks.append({'name':name,'passed':residual==0,'residual':str(residual)})
main=root/'experiments/eft_matching_v01/code/run_matching.py'
protocol=root/'experiments/eft_matching_v01/protocol.md'
res={'scope':'Only local-center first jets entering selected X and baseline B interference; no full off-axis valley, loops or new cosmic solution. Symbolic W/M jets are independent of production implementation.','source_hashes':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (main,protocol,Path(__file__))},'checks':checks,'passed':sum(c['passed'] for c in checks),'total':len(checks),'sympy':s.__version__}
Path(__file__).with_suffix('.json').write_text(json.dumps(res,ensure_ascii=False,indent=2)+'\n')
print(res['passed'], '/',res['total'])
if res['passed']!=res['total']:raise SystemExit(1)
