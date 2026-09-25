#!/usr/bin/env python3
"""Specified geometric protection: local asymptotic budgets, not new cosmic solutions."""
from pathlib import Path
from datetime import datetime,timezone
import csv,hashlib,json,platform,time
import numpy as np
from scipy.optimize import brentq
import scipy
ROOT=Path(__file__).resolve().parents[1]; REPO=ROOT.parents[1];OUT=ROOT/'results'
PHASES={'zero':(1.,0.),'quadrature':(0.,1.),'pi':(-1.,0.)}
MASSES=np.logspace(-12,12,97);COEFFS=(.1,1.,10.);SCALES=(.5,1.,2.)
REP=(1e-6,1.,1e5,1e7,1e12);MI=1e11

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def csvout(p,rows):
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
def inputs():
 d=json.loads((REPO/'experiments/sugra_shift_v01/results/inputs.json').read_text())
 with (REPO/'experiments/sugra_shift_v01/results/trajectories.csv').open() as f:
  rows=[r for r in csv.DictReader(f) if r['case']=='shift_axis']
 arr={k:np.array([float(r[k]) for r in rows]) for k in ('chi','z','E','qx')}
 d['rho_lambda_eV4']=3*d['omega_lambda']*d['Mpl_eV']**2*d['H_ref_eV']**2
 return d,arr

def cw(s,der_s,d,der_d,h2,der_h2):
 """Small-split finite supertrace, retaining d^3,h^4, derivative in either direction."""
 L=np.log(s/MI**2);t=s+d;Lt=L+np.log1p(d/s)
 val=4*s*d*(L-1)+2*d*d*L+2*d**3/(3*s)+2*h2*Lt-h2*h2/(6*t*t)
 der=4*(s*der_d*(L-1)+d*der_s*L)+4*d*der_d*L+2*d*d*der_s/s+2*d*d*der_d/s-2*d**3*der_s/(3*s*s)
 der+=2*der_h2*Lt+2*h2*(der_s+der_d)/t-h2*der_h2/(3*t*t)+h2*h2*(der_s+der_d)/(3*t**3)
 return val/(32*np.pi**2),der/(32*np.pi**2)

def calculate(mg,a,phase,mu_factor,scheme,par,old):
 P=par['Mpl_eV']**2;F=np.sqrt(par['F_squared_eV2']);rho=par['rho_lambda_eV4']
 chi=old['chi'];co,si=PHASES[phase]
 U=par['potential_unit_eV4']*par['W_i']*np.exp(-2*(chi-par['chi_i']))
 wm=F*np.sqrt(U)/np.sqrt(2);rw=wm/(mg*P)
 h2=(co-rw)**2+si**2
 # Avoid subtracting 1 from h2 for the infinitesimal change entering log.
 logh=np.log1p(-2*co*rw+rw*rw)
 hx=2*rw*(co-rw);hy=2*rw*si/3
 X=mg*mg*h2;Xx=mg*mg*hx;Xy=mg*mg*hy
 c=48*a;L0=-2*np.log(mu_factor)
 vmod0=2*c*c*mg**4*(L0-1.5)/(64*np.pi**2)
 eta=rho/(6*mg*mg*P) if scheme=='tree_eta' else (rho-vmod0)/(6*mg*mg*P)
 # Leading expansion around the exact geometric vacuum; the independent calculation solves full K/W.
 delta=eta/(6*a)+U/(36*a*mg*mg*P*h2)
 delta_x=U/(36*a*mg*mg*P*h2)*(-2-hx/h2)
 delta_y=-U*hy/(36*a*mg*mg*P*h2*h2)
 q=2*eta+12*a*delta**2;qx=24*a*delta*delta_x;qy=24*a*delta*delta_y
 f=1+delta+eta*delta**2+a*delta**4;p=1+2*eta*delta+4*a*delta**3;D=p*p-f*q
 fx=p*delta_x;fy=p*delta_y
 Dx=(p*q-24*a*f*delta)*delta_x;Dy=(p*q-24*a*f*delta)*delta_y
 Wbar=mg*P*((co-rw)-1j*si)
 Vtree=U/f**2+3*q*mg*mg*P*h2/(f*f*D)
 # Envelope derivatives; explicit off-axis S gives the 1/3 factor in y.
 Vx=-2*U/f**2+3*q*P*Xx/(f*f*D)
 Vy=3*q*P*Xy/(f*f*D)
 extra_x=Vx+2*U;extra_y=Vy
 # The small correction is below double precision; reconstruct it rather than subtract nearly equal forces.
 extra_x=2*U*(-np.expm1(-2*np.log1p(delta+eta*delta**2+a*delta**4)))+3*q*P*Xx/(f*f*D)
 s0=.5*MI**2*(1+par['xi']/chi**2);M=np.sqrt(s0)
 r=-par['xi']/(chi*(chi**2+par['xi']));rp=par['xi']*(3*chi**2+par['xi'])/(chi**2*(chi**2+par['xi'])**2)
 s=s0/f;sx=s*(2*r-fx/f);sy=-s*fy/f
 d=2*Vtree/(3*P);dx=2*Vx/(3*P);dy=2*Vy/(3*P)
 hclock=np.sqrt(2*U)/F
 bh=q*Wbar/(P*D)
 B=M*(hclock*r+bh)/f
 Bx=M/f*(hclock*(r*r+rp-r)+(r-fx/f)*bh+((qx-q*Dx/D)*Wbar+q*wm)/(P*D)) - M*hclock*r*fx/f**2
 # y derivative of the full off-axis soft expression, plus the tiny valley displacement chain rule.
 By=1j*M/f*(hclock*(r*r+rp+r)+q*(-wm/3-r*Wbar)/(P*D))
 By+=M/f*((qy-q*Dy/D)*Wbar/(P*D)-(hclock*r+bh)*fy/f)
 b2=abs(B)**2;b2x=2*np.real(np.conj(B)*Bx);b2y=2*np.real(np.conj(B)*By)
 vq,vqx=cw(s,sx,d,dx,b2,b2x);_,vqy=cw(s,sy,d,dy,b2,b2y)
 # Local heavy-scalar CW only; use the full-K/W leading heavy eigenvalue's y jet.
 L=L0+logh
 common=2*c*c*X*(L-1)/(32*np.pi**2)
 vtx=common*Xx;vty=common*Xy
 scalar_vac=2*c*c*X*X*(L-1.5)/(64*np.pi**2)
 fcomp_abs=mg*np.sqrt(h2)*q/(np.sqrt(f)*D)
 g2=4*np.pi*par['alpha_reference'];amsb_abs=2*g2*fcomp_abs/(16*np.pi**2)
 denom=2*U
 norm=lambda x,y:np.hypot(x,y)/denom
 rt=norm(extra_x,extra_y);rq=norm(vqx,vqy);rm=norm(vtx,vty)
 selected_x=extra_x+vqx+vtx;selected_y=extra_y+vqy+vty
 selected=norm(selected_x,selected_y);component_sum=rt+rq+rm
 mass=np.sqrt(c*X);cutoff=par['Mpl_eV']/np.sqrt(a)
 kinetic=.5*F*F*(old['E']*par['H_ref_eV']*old['qx'])**2
 delta_kinetic=-kinetic/(72*a*mg*mg*P*h2)
 arrays=dict(U=U,eta=np.full_like(U,eta),delta=delta,delta_kinetic=delta_kinetic,
  s=s,d=d,B_abs2=b2,B_real=B.real,B_imag=B.imag,Vtree=Vtree,Vtree_x=Vx,Vtree_y=Vy,
  extra_tree_x=extra_x,extra_tree_y=extra_y,charged_CW=vq,charged_CW_x=vqx,charged_CW_y=vqy,
  modulus_CW_x=vtx,modulus_CW_y=vty,modulus_CW_vac=scalar_vac,
  selected_force_ratio=selected,component_sum_ratio=component_sum,
  tree_force_ratio=rt,charged_force_ratio=rq,modulus_force_ratio=rm,
  modulus_mass_eV=mass,Fcomp_abs_eV=fcomp_abs,above_threshold_AMSB_abs_eV=amsb_abs)
 rec=dict(mG_eV=float(mg),a=float(a),phase=phase,mu_factor=float(mu_factor),scheme=scheme,
  eta=float(eta),scalar_CW_clockoff_eV4=float(vmod0),max_abs_delta=float(max(abs(delta))),
  max_abs_delta_kinetic=float(max(abs(delta_kinetic))),min_f=float(min(f)),min_D=float(min(D)),
  max_tree_force_ratio=float(max(rt)),max_charged_force_ratio=float(max(rq)),max_modulus_force_ratio=float(max(rm)),
  max_selected_force_ratio=float(max(selected)),max_selected_index=int(np.argmax(selected)),
  max_component_sum_ratio=float(max(component_sum)),max_component_sum_index=int(np.argmax(component_sum)),
  max_selected_x_ratio=float(max(abs(selected_x)/denom)),max_selected_y_ratio=float(max(abs(selected_y)/denom)),
  min_modulus_mass_eV=float(min(mass)),max_H_over_modulus_mass=float(max(old['E']*par['H_ref_eV']/mass)),
  max_mass_over_geometric_marker=float(max(MI,max(mass))/cutoff),geometric_marker_eV=float(cutoff),
  max_d_over_s=float(max(abs(d/s))),max_B_over_s=float(max(abs(B)/s)),
  min_scalar_over_fermion_mass2=float(min(1+d/s-abs(B)/s)),
  max_Fcomp_abs_eV=float(max(fcomp_abs)),max_above_threshold_AMSB_abs_eV=float(max(amsb_abs)),
  selected_budget_pass=bool(max(selected)<=.1),component_sum_budget_pass=bool(max(component_sum)<=.1))
 return rec,arrays

def main():
 start=time.time();OUT.mkdir(exist_ok=True,parents=True);par,old=inputs();rows=[];reps=[];checks=[];limits=[]
 def check(name,val,tol):checks.append(dict(name=name,value=float(val),tolerance=float(tol),passed=bool(np.isfinite(val) and val<tol)))
 for a in COEFFS:
  for phase in PHASES:
   for mg in MASSES:
    for mu in SCALES:
     for scheme in ('tree_eta','scalar_retuned_eta'):
      rec,_=calculate(float(mg),a,phase,mu,scheme,par,old);rows.append(rec)
   for mu in SCALES:
    for scheme in ('tree_eta','scalar_retuned_eta'):
     def objective(lm):return np.log10(calculate(10**lm,a,phase,mu,scheme,par,old)[0]['max_component_sum_ratio'])+1
     lm=brentq(objective,-12,12,xtol=1e-11)
     rec,_=calculate(10**lm,a,phase,mu,scheme,par,old)
     limits.append(dict(a=a,phase=phase,mu_factor=mu,scheme=scheme,component_sum_10percent_mG_eV=10**lm,
       selected_force_at_component_bound=rec['max_selected_force_ratio']))
     check(f'budget_root_{a}_{phase}_{mu}_{scheme}',abs(rec['max_component_sum_ratio']/.1-1),1e-8)
 for phase in PHASES:
  for mg in REP:
   _,t=calculate(mg,1.,phase,1.,'tree_eta',par,old)
   _,r=calculate(mg,1.,phase,1.,'scalar_retuned_eta',par,old)
   for n in range(len(old['chi'])):
    item=dict(phase=phase,mG_eV=mg,index=n,chi=float(old['chi'][n]),z=float(old['z'][n]))
    for prefix,arr in [('tree_eta',t),('retuned',r)]:
     for k,v in arr.items():item[prefix+'_'+k]=float(v[n])
    reps.append(item)
 check('positive_geometry',0 if all(r['min_D']>0 and r['min_f']>0 for r in rows) else 1,.5)
 check('charged_spectra_positive',0 if all(r['min_scalar_over_fermion_mass2']>0 for r in rows) else 1,.5)
 check('small_charged_splitting',max(max(r['max_d_over_s'],r['max_B_over_s']) for r in rows),.01)
 check('small_modulus_displacement',max(r['max_abs_delta'] for r in rows),1e-8)
 check('heavy_modulus_adiabatic',max(r['max_H_over_modulus_mass'] for r in rows),1e-3)
 check('geometric_marker_margin',max(r['max_mass_over_geometric_marker'] for r in rows),.1)
 # Generic phases have y/chi gradients in the 1/3 ratio, not an artificial quadrature exemption.
 for mg in (1e-6,1.,1e7):
  _,z=calculate(mg,1.,'zero',1.,'tree_eta',par,old);_,q=calculate(mg,1.,'quadrature',1.,'tree_eta',par,old)
  check(f'heavy_y_jet_factor_{mg}',max(abs(3*q['modulus_CW_y']/z['modulus_CW_x']-1)),1e-8)
  check(f'charged_1ev_small_{mg}',max(z['charged_force_ratio']),1e-20)
 csvout(OUT/'scheme_scan.csv',rows);csvout(OUT/'representative_curves.csv',reps);csvout(OUT/'budget_limits.csv',limits)
 primary=[r for r in rows if r['mu_factor']==1 and r['scheme']=='tree_eta'];csvout(OUT/'scan.csv',primary)
 example=[r for r in rows if r['a']==1 and r['mu_factor']==1 and r['mG_eV'] in (1e-6,1.,1e5,1e7,1e12)]
 dump(OUT/'inputs.json',dict(parent_inputs=par,masses_eV=MASSES.tolist(),coefficients=COEFFS,phases=PHASES,
  scales=SCALES,representative_masses=REP,heavy_mass_initial_eV=MI,force_budget=.1,
  scientific_scope='Leading adiabatic geometry; tree_eta and scalar-vacuum-retuned diagnostic, not full SUGRA loops or new trajectories'))
 summary=dict(created_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=time.time()-start,
  python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
  source_hashes={str(p.relative_to(REPO)):sha(p) for p in [REPO/'experiments/sugra_shift_v01/results/inputs.json',REPO/'experiments/sugra_shift_v01/results/trajectories.csv',ROOT/'protocol.md',Path(__file__)]},
  scan_rows=len(primary),scheme_rows=len(rows),representative_rows=len(reps),limit_rows=len(limits),
  checks=checks,checks_passed=sum(c['passed'] for c in checks),check_count=len(checks),examples=example,
  limits=limits,scope='Selected local finite contributions only; no complete loop protection or cosmological fit is claimed.')
 dump(OUT/'summary.json',summary)
 print('Checks',summary['checks_passed'],'/',summary['check_count'],flush=True)
 print('Rows',len(primary),len(rows),len(reps),flush=True)
 for r in limits:
  if r['a']==1 and r['mu_factor']==1: print(r,flush=True)
 if not all(c['passed'] for c in checks):raise SystemExit(1)

if __name__=='__main__':main()
