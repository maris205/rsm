#!/usr/bin/env python3
"""Independent full K/W/VD block check; no main-stage implementation imports."""
from pathlib import Path
import json, hashlib, argparse
import mpmath as mp
HERE=Path(__file__).resolve()
parser=argparse.ArgumentParser();parser.add_argument('--dps',type=int,default=220);parser.add_argument('--output');args=parser.parse_args()
mp.mp.dps=args.dps
checks=[]
def ck(name,a,b,tol='1e-100',scale=None):
    den=max(abs(b),mp.mpf('1e-120')) if scale is None else mp.mpf(scale)
    err=abs(a-b)/den
    checks.append(dict(name=name,actual=mp.nstr(a,40),expected=mp.nstr(b,40),relative_error=mp.nstr(err,12),tolerance=tol,passed=bool(err<mp.mpf(tol))))
def H(v):
    n=len(v)
    return mp.matrix([10*v[i]-(3*v[i-1] if i else 0)-(3*v[i+1] if i<n-1 else 0) for i in range(n)])
def Hi(v):
    n=len(v);d=[mp.mpf(10)]*n;r=list(v)
    for i in range(1,n):
        fac=-3/d[i-1];d[i]-=fac*(-3);r[i]-=fac*r[i-1]
    x=[mp.mpf(0)]*n;x[-1]=r[-1]/d[-1]
    for i in range(n-2,-1,-1):x[i]=(r[i]+3*x[i+1])/d[i]
    return mp.matrix(x)
def dot(a,b):return sum(a[i]*b[i] for i in range(len(a)))
def geom(t,eta,A):
    z=t-1
    return t+eta*z*z+z**4-A/t**2,1+2*eta*z+4*z**3+2*A/t**3,2*eta+12*z*z-6*A/t**4,24*z+24*A/t**5

def ref(P,mg,A,rho):
    t0=mp.findroot(lambda t:(t-1)*t**5+A,(mp.mpf('.995'),mp.mpf('1.005')),tol=mp.mpf('1e-140'))
    eta0=-6*(t0-1)**2+3*A/t0**4
    rr=rho/(3*P*mg*mg)
    def eq(t,e):
        g,p,q,qq=geom(t,e,A)
        return q-rr*p*p/(g*(1+rr)),qq-2*p*q/g+q*q/p
    t,e=mp.findroot(eq,(t0,eta0),tol=mp.mpf('1e-140'))
    g,p,q,_=geom(t,e,A)
    return dict(t=t,eta=e,g=g,p=p,q=q,h=1+A/t**3,hp=-3*A/t**4,hpp=12*A/t**5,GT=3*P*(p*p-g*q)/g**2)

def full(y,j,P,mg,A,Fp,U,theta,M,r,shift=0,return_matrix=False,tshift=0,vshift=0):
    """Full chiral metric Schur complement at T=T*, any y and bottom j."""
    n=len(j);gs,p,q,hs,hp,hpp,GT=(r[k] for k in ('g','p','q','h','hp','hpp','GT'))
    c=gs/(3*P);a=mp.sqrt(2)/M;F=Fp*mp.sqrt(gs/hs)
    t=r['t']+tshift;gv,gt,gtt,_=geom(t,r['eta'],A);gv+=vshift**4
    ht=1+A/t**3;hp=-3*A/t**4;hpp=12*A/t**5
    xc2=GT*(tshift*tshift+vshift*vshift)/4
    hj=H(j);R=dot(j,hj);b=a*j[len(j)-1];kc=Fp**2*y*y
    u=mp.matrix(hj);u[0]+=a*kc;u[n-1]+=a*xc2
    Om=gv-c*(R/2+a*j[0]*kc+b*xc2)-ht*F*F*y*y/(3*P)
    HH=ht+hs*a*j[0]
    oi=mp.matrix([gt-4j*vshift**3-hp*F*F*y*y/(3*P)-c*b*GT*(tshift-1j*vshift)/2,1j*mp.sqrt(2)*F*HH*y/(3*P)])
    E=mp.matrix([[gtt+12*vshift*vshift-hpp*F*F*y*y/(3*P)-c*b*GT,-1j*mp.sqrt(2)*hp*F*y/(3*P)],
                 [1j*mp.sqrt(2)*hp*F*y/(3*P),-HH/(3*P)]])
    KL=mp.matrix(2)
    for i in range(2):
        for k in range(2):KL[i,k]=-3*P*(E[i,k]/Om-oi[i]*mp.conj(oi[k])/Om**2)
    B=[-gs*oi[0]*u/Om**2,-gs*oi[1]*u/Om**2]
    B[0][n-1]+=gs*a*GT*(tshift-1j*vshift)/(2*Om)
    B[1][0]+=-1j*mp.sqrt(2)*gs*a*Fp*Fp*y/(F*Om)
    hiu=Hi(u);uu=dot(u,hiu)
    def Dinv(v):return Om/gs*(Hi(v)-c*hiu*dot(hiu,v)/(Om+c*uu))
    kd=-3*P*oi/Om;kh=gs*u/Om
    wr=-gs**mp.mpf('1.5')*Fp*mp.sqrt(U/2)*mp.exp(-shift-1j*y)
    W=P*gs**mp.mpf('1.5')*mg*mp.exp(1j*theta)+wr
    wz=-mp.sqrt(2)*wr/F
    DL=kd*W/P+mp.matrix([0,wz]);DH=kh*W/P
    diDH=Dinv(DH);s=mp.matrix(KL);de=mp.matrix(DL)
    for i in range(2):
        de[i]-=dot(B[i],diDH)
        for k in range(2):s[i,k]-=dot(B[i],Dinv(B[k].conjugate()))
    norm=(de.conjugate().T*(s**-1)*de)[0]+dot(DH.conjugate(),diDH)
    VF=mp.re(Om**-3*(norm-3*abs(W)**2/P))
    VD=M*M*gs*gs*dot(u,u)/(4*Om*Om)
    out=dict(V=VF+VD,VF=VF,VD=VD,Omega=Om,Schur=s,W=W)
    if return_matrix:
        K=mp.matrix(n+2);DW=mp.matrix(list(DL)+list(DH))
        for i in range(2):
            for k in range(2):K[i,k]=KL[i,k]
            for k in range(n):K[i,k+2]=B[i][k];K[k+2,i]=mp.conj(B[i][k])
        for i in range(n):
            for k in range(n):K[i+2,k+2]=gs/Om*((10 if i==k else -3 if abs(i-k)==1 else 0)+c*u[i]*u[k]/Om)
        out['directV']=mp.re(Om**-3*((DW.conjugate().T*(K**-1)*DW)[0]-3*abs(W)**2/P))+VD
    return out

def reduced(j,P,mg,Fp,U,theta,M,r,shift=0):
    g,p,q,hs,GT=(r[k] for k in ('g','p','q','h','GT'))
    a=mp.sqrt(2)/M;c=g/(3*P);R=dot(j,H(j));Om=g-c*R/2;b=a*j[len(j)-1]
    C=q-c*b*GT;W0=P*g**mp.mpf('1.5')*mg
    rw=Fp*mp.sqrt(U/2)/(P*mg)*mp.exp(-shift)
    w2=W0**2*(1-2*rw*mp.cos(theta)+rw*rw)
    den=p*p-C*(Om+c*R)
    vf=3*w2*C/(P*Om*Om*den)+U*mp.exp(-2*shift)*g*g/(Om*Om*(1+a*j[0]))
    return vf+M*M*g*g*dot(H(j),H(j))/(4*Om*Om)

# Dimensionless, non-no-scale/finite-y tests verify actual complete block inversion.
for n in (2,5):
    P,mg,A,Fp,U,M=map(mp.mpf,('3','.2','.001','.13','.0003','2'))
    rr=ref(P,mg,A,mp.mpf('.0001'))
    j=mp.matrix([mp.mpf('.003')*(i+1) for i in range(n)])
    for theta in (mp.mpf(0),mp.mpf('.6')):
        for y in (mp.mpf(0),mp.mpf('.02')):
            o=full(y,j,P,mg,A,Fp,U,theta,M,rr,return_matrix=True)
            ck(f'block_vs_direct_n{n}_phase{theta}_y{y}',o['V'],o['directV'])
            if y==0:ck(f'closed_axis_n{n}_phase{theta}',o['V'],reduced(j,P,mg,Fp,U,theta,M,rr))

n=128;e=mp.matrix([0]*(n-1)+[1]);u=Hi(e);v=Hi(u);e0=mp.matrix([1]+[0]*(n-1));u0=Hi(e0)
H00=u0[0];hee=dot(u,u)
ck('chain_H_inverse_residual',mp.norm(H(u)-e),0,scale=1)
ck('chain_H_squared_inverse_residual',mp.norm(H(H(v))-e),0,scale=1)
P=mp.mpf('2.435e27')**2;Fp=mp.sqrt(mp.mpf('5.929225e50'));rho=mp.mpf('2.5181378331717977e-11')
A=mp.zeta(3)/(48*mp.pi**4)*mp.mpf('1e14')**2/P
Lam=mp.mpf('1e12');M=Lam*mp.sqrt(H00);Ui=mp.mpf('6.161810271353063e-14')
chi0=mp.mpf('137.6422999547138');chi1=mp.mpf('139.79424666735443')
rows=[]
radial=[]
for n,Lam in ((128,mp.mpf('1e12')),(127,mp.mpf('2e12'))):
    e=mp.matrix([0]*(n-1)+[1]);e0=mp.matrix([1]+[0]*(n-1));u=Hi(e);v=Hi(u);u0=Hi(e0)
    H00=u0[0];hee=dot(u,u);M=Lam*mp.sqrt(H00)
    for mg in (mp.mpf('1e-6'),mp.mpf('1e-14'),mp.mpf('1e-15')):
        rr=ref(P,mg,A,rho)
        for chi in (chi0,chi1):
            U=Ui*mp.exp(-2*(chi-chi0))
            for phase,theta in (('aligned',mp.mpf(0)),('quadrature',mp.pi/2)):
                rw=Fp*mp.sqrt(U/2)/(P*mg);S=3*P*mg*mg*(1-2*rw*mp.cos(theta)+rw*rw)
                kappa=4*S*hee/M**4;x=mp.findroot(lambda x:x*(1+x)**2-kappa,(kappa/(1+kappa)**2,kappa),tol=mp.mpf('1e-140'))
                j=2*mp.sqrt(2)*S/(M**3*(1+x)**2)*v;jzero=mp.matrix(n,1)
                dV=lambda yy,ss:full(yy,j,P,mg,A,Fp,U,theta,M,rr,shift=ss)['V']-full(yy,jzero,P,mg,A,Fp,U,theta,M,rr,shift=ss)['V']
                delta=dV(0,0)
                axis=reduced(j,P,mg,Fp,U,theta,M,rr)-reduced(jzero,P,mg,Fp,U,theta,M,rr)
                tag=f'n{n}_L{Lam}_mg{mg}_chi{chi}_phase{phase}'
                ck(tag+'_full_vs_closed',delta,axis,tol='1e-85')
                SX=6*P*mg*mg*rw*(mp.cos(theta)-rw)
                SYnaive=6*P*mg*mg*rw*mp.sin(theta)
                predV=S*((1+mp.mpf('1.5')*x)/(1+x)**2-1)
                predX=-x/(1+x)*SX
                dx=mp.diff(lambda sh:dV(0,sh),0)
                dy=mp.diff(lambda yy:dV(yy,0),0)
                # y prediction is independently derived after full K metric reduction.
                row=dict(n=n,Lambda_eV=mp.nstr(Lam,25),mg_eV=mp.nstr(mg,25),chi=mp.nstr(chi,25),phase=phase,x=mp.nstr(x,30),kappa=mp.nstr(kappa,30),deltaV_eV4=mp.nstr(delta,40),rigid_deltaV_eV4=mp.nstr(predV,40),deltaV_chi_eV4=mp.nstr(dx,40),rigid_deltaV_chi_eV4=mp.nstr(predX,40),deltaV_y_eV4=mp.nstr(dy,40),naive_S_y=mp.nstr(SYnaive,40),force_chi_over_2U=mp.nstr(abs(dx)/(2*U),40),force_y_over_2U=mp.nstr(abs(dy)/(2*U),40),R_over_P=mp.nstr(dot(j,H(j))/P,30))
                ck(tag+'_rigid_axis_force',dx,predX,tol='1e-24')
                ck(tag+'_full_y_factor_one_third',dy,-x/(1+x)*SYnaive/3,tol='1e-24',scale=max(abs(x/(1+x)*SYnaive/3),mp.mpf('1e-80')))
                # Exact reduced full-SUGRA heavy force at the leading stationary bottom.
                gs,pp,qq,GT=(rr[k] for k in ('g','p','q','GT'));c=gs/(3*P);aa=mp.sqrt(2)/M
                hj=H(j);h2j=H(hj);R=dot(j,hj);E=dot(hj,hj);Om=gs-c*R/2;bb=aa*j[n-1];ll=aa*j[0]
                CC=qq-c*bb*GT;den=pp*pp-CC*(Om+c*R);w2=P*P*gs**3*mg*mg*(1-2*rw*mp.cos(theta)+rw*rw)
                vf=3*w2*CC/(P*Om*Om*den);vz=U*gs*gs/(Om*Om*(1+ll));vd=M*M*gs*gs*E/(4*Om*Om)
                vr=c/Om*(vf+vz+vd)+vf*CC*c/(2*den);ve=M*M*gs*gs/(4*Om*Om)
                vb=-c*GT*3*w2*pp*pp/(P*Om*Om*den*den);vl=-vz/(1+ll)
                grad=2*vr*hj+2*ve*h2j+aa*vb*e+aa*vl*e0
                ck(tag+'_full_heavy_stationarity_residual',mp.norm(grad),0,tol='1e-24',scale=mp.norm(2*ve*h2j))
                row['heavy_stationarity_residual_fraction']=mp.nstr(mp.norm(grad)/mp.norm(2*ve*h2j),30)
                ck(tag+'_rigid_potential',delta,predV,tol='1e-24')
                if chi==chi1 and phase=='aligned':
                    fn=lambda tt,vv:full(0,j,P,mg,A,Fp,U,theta,M,rr,tshift=tt,vshift=vv)['V']
                    kval=mp.re(full(0,j,P,mg,A,Fp,U,theta,M,rr)['Schur'][0,0])
                    mt=2*mp.diff(lambda tt:fn(tt,0),0,2)/kval
                    mv=2*mp.diff(lambda vv:fn(0,vv),0,2)/kval
                    C=2*S/(M*M*(1+x)**3);pred_mass=2*C*H00
                    ft=mp.sqrt(2/kval)*mp.diff(lambda tt:fn(tt,0),0)
                    fv=mp.sqrt(2/kval)*mp.diff(lambda vv:fn(0,vv),0)
                    ck(tag+'_bulk_modulus_radial_mass',mt,pred_mass,tol='1e-24')
                    ck(tag+'_bulk_modulus_imaginary_mass',mv,pred_mass,tol='1e-24')
                    ck(tag+'_bulk_modulus_imaginary_force',fv,0,tol='1e-100',scale=1)
                    radial.append(dict(n=n,Lambda_eV=str(Lam),mg_eV=str(mg),mass2_t_eV2=mp.nstr(mt,40),mass2_v_eV2=mp.nstr(mv,40),rigid_mass2_eV2=mp.nstr(pred_mass,40),canonical_force_t_eV3=mp.nstr(ft,40),linear_center_shift_eV=mp.nstr(-ft/mt,40),linear_shift_over_M=mp.nstr(abs(ft/mt)/M,40)))
                rows.append(row)
# Once-only local-curvature rematch; exact leading normal-coordinate SUGRA, not radial stationarity.
rematch=[]
for kap in (mp.mpf('1e-4'),mp.mpf('1e-20'),mp.mpf('.1')):
    alpha_x=lambda x:(1+mp.mpf('1.5')*x)/(1+x)**2
    x=mp.findroot(lambda x:alpha_x(x)**3*x*(1+x)**2-kap,(kap,kap*mp.mpf('1.1')),tol=mp.mpf('1e-140'))
    alpha=alpha_x(x);b=alpha*x;S=kap/2;d=mp.mpf(1)
    V=S*(1/(alpha+b)-1)+b*b/(4*d)
    coeff=1/(alpha+b)-1
    ck('rematch_vacuum_'+str(kap),V,0,scale=1)
    ck('rematch_stationary_'+str(kap),b*(alpha+b)**2,2*d*S)
    ck('rematch_residual_force_'+str(kap),coeff,-x/(2+3*x))
    rematch.append(dict(kappa=str(kap),x=mp.nstr(x,35),alpha=mp.nstr(alpha,35),force_coefficient=mp.nstr(coeff,35)))
precision_file=HERE.with_name('mediator_sugra_response_20260926_150_checks.json')
precision_checks=[]
if args.dps==220 and precision_file.exists():
    low=json.loads(precision_file.read_text())
    if low.get('script_sha256')!=hashlib.sha256(HERE.read_bytes()).hexdigest():
        raise RuntimeError('150-digit baseline must use the same final checker source')
    for table,high_rows,keys in (
        ('physical_rows',rows,('x','deltaV_eV4','deltaV_chi_eV4','deltaV_y_eV4','force_chi_over_2U','force_y_over_2U','R_over_P')),
        ('radial_rows',radial,('mass2_t_eV2','mass2_v_eV2','canonical_force_t_eV3','linear_center_shift_eV','linear_shift_over_M'))):
        for i,(hi,lo) in enumerate(zip(high_rows,low[table],strict=True)):
            for key in keys:
                ck(f'precision150_220_{table}_{i}_{key}',mp.mpf(hi[key]),mp.mpf(lo[key]),tol='1e-35',scale=max(abs(mp.mpf(lo[key])),mp.mpf('1e-100')))
    precision_checks=[c for c in checks if c['name'].startswith('precision150_220_')]
res=dict(passed=sum(c['passed'] for c in checks),total=len(checks),dps=mp.mp.dps,script_sha256=hashlib.sha256(HERE.read_bytes()).hexdigest(),chain_cases=[dict(n=128,q=3,Lambda_eV='1e12'),dict(n=127,q=3,Lambda_eV='2e12')],cross_precision_checks=len(precision_checks),physical_rows=rows,radial_rows=radial,curvature_rematch=rematch,checks=checks,scope='Full K/W and D-term plus local radial Hessian, fixed heavy bottom at leading stationary value. No solved full radial extremum or mediator loops.')
destination=Path(args.output) if args.output else HERE.with_suffix('.json')
destination.write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps(dict(passed=res['passed'],total=res['total'],rows=rows),indent=2))
if res['passed']!=res['total']:raise SystemExit(1)
