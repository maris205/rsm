#!/usr/bin/env python3
"""Independent full K/W check of a locally centered persistent shared current.

No experiment implementation or earlier checker is imported. All two-field
metric, auxiliaries, scalar potential and charged soft jets follow from Omega.
"""
from pathlib import Path
import hashlib
import json
import mpmath as mp

mp.mp.dps = 100
HERE = Path(__file__).resolve()
checks = []


def ck(name, actual, expected, scale=None, tol='1e-70'):
    den = max(abs(expected), mp.mpf('1e-65')) if scale is None else mp.mpf(scale)
    err = abs(actual-expected)/den
    checks.append(dict(name=name, actual=mp.nstr(actual, 30), expected=mp.nstr(expected, 30),
                       error=mp.nstr(err, 12), tolerance=tol, passed=bool(err < mp.mpf(tol))))


def geometry(t, eta, A):
    d = t-1
    return (t+eta*d*d+d**4-A/t**2,
            1+2*eta*d+4*d**3+2*A/t**3,
            2*eta+12*d*d-6*A/t**4,
            24*d+24*A/t**5)


def reference(P, mg, A, rho, hqfactor=2):
    t0 = mp.findroot(lambda t:(t-1)*t**5+A, (mp.mpf('.995'),mp.mpf('1.005')), tol=mp.mpf('1e-95'))
    eta0 = -6*(t0-1)**2+3*A/t0**4
    if rho:
        rr = rho/(3*P*mg*mg)
        def eq(t, eta):
            g,p,q,qp = geometry(t,eta,A)
            return q-rr*p*p/(g*(1+rr)), qp-2*p*q/g+q*q/p
        t,eta = mp.findroot(eq,(t0,eta0),tol=mp.mpf('1e-90'))
    else:
        t,eta=t0,eta0
    g,p,q,_=geometry(t,eta,A)
    D=p*p-g*q
    return dict(tstar=t,eta=eta,gs=g,ps=p,qs=q,Ds=D,hs=1+A/t**3,
                hqs=1+hqfactor*A/t**3,GT=3*P*D/g**2)


def full(t, v, y, *, P, A, Fp, mg, U, phase, L, tau, sigma, ref,
         strength=1, chishift=mp.mpf(0), hqfactor=2):
    """Exact finite-field two-chiral K/W and charged quadratic spectrum."""
    ts,eta,gs,hs,hqs,GT=(ref[k] for k in ('tstar','eta','gs','hs','hqs','GT'))
    F=Fp*mp.sqrt(gs/hs)
    root=mp.sqrt(2)*F
    g,gt,gtt,_=geometry(t,eta,A)
    g+=v**4
    gv,gvv=4*v**3,12*v*v
    h,hp,hpp=1+A/t**3,-3*A/t**4,12*A/t**5
    hq,hqp,hqpp=1+hqfactor*A/t**3,-3*hqfactor*A/t**4,12*hqfactor*A/t**5
    J=Fp**2*y*y+tau*GT*((t-ts)**2+v*v)/4
    Ji=mp.matrix([tau*GT*(t-ts)/2,tau*GT*v/2,2*Fp**2*y])
    Jij=mp.diag([tau*GT/2,tau*GT/2,2*Fp**2])
    b=strength*gs/(3*P*L*L)
    Om=g-h*F*F*y*y/(3*P)+b*J*J
    Oi=mp.matrix([gt-hp*F*F*y*y/(3*P),gv,-2*h*F*F*y/(3*P)])+2*b*J*Ji
    Oij=mp.matrix([[gtt-hpp*F*F*y*y/(3*P),0,-2*hp*F*F*y/(3*P)],
                   [0,gvv,0],[-2*hp*F*F*y/(3*P),0,-2*h*F*F/(3*P)]])
    for i in range(3):
        for j in range(3):
            Oij[i,j]+=2*b*(Ji[i]*Ji[j]+J*Jij[i,j])
    Ki=-3*P*Oi/Om
    Kij=mp.matrix(3)
    for i in range(3):
        for j in range(3):
            Kij[i,j]=-3*P*(Oij[i,j]/Om-Oi[i]*Oi[j]/Om**2)
    metric=mp.matrix([[Kij[0,0]+Kij[1,1],(1j*Kij[0,2]+Kij[1,2])/root],
                      [(-1j*Kij[0,2]+Kij[1,2])/root,Kij[2,2]/root**2]])
    wr=-gs**mp.mpf('1.5')*Fp*mp.sqrt(U/2)*mp.exp(-chishift)
    Wc=P*gs**mp.mpf('1.5')*mg*mp.exp(1j*phase)
    W=Wc+wr*mp.exp(-1j*y)
    Wz=-mp.sqrt(2)*wr/F*mp.exp(-1j*y)
    kt=Ki[0]-1j*Ki[1]
    kz=-1j*Ki[2]/root
    DW=mp.matrix([kt*W/P,Wz+kz*W/P])
    inv=metric**-1
    aux=-Om**mp.mpf('-1.5')*inv.T*DW.conjugate()
    V=mp.re(Om**-3*((DW.conjugate().T*inv*DW)[0]-3*abs(W)**2/P))
    # The shared current is defined with frozen matching coefficients.
    c=strength*2*sigma*hqs/(L*L)
    H=hq-c*J
    Hi=mp.matrix([hqp,0,0])-c*Ji
    Hij=mp.diag([hqpp,0,0])-c*Jij
    elli=Hi/H
    ellij=mp.matrix(3)
    for i in range(3):
        for j in range(3):
            ellij[i,j]=Hij[i,j]/H-Hi[i]*Hi[j]/H**2
    ell=mp.matrix([[ellij[0,0]+ellij[1,1],(1j*ellij[0,2]+ellij[1,2])/root],
                   [(-1j*ellij[0,2]+ellij[1,2])/root,ellij[2,2]/root**2]])
    d=2*V/(3*P)-mp.re(sum(aux[i]*mp.conj(aux[j])*ell[i,j] for i in range(2) for j in range(2)))
    zeta=mp.mpf(140)+chishift+1j*y
    xi=mp.mpf(19000)
    Mhol=mp.sqrt(gs)*hqs*mp.sqrt((1+xi/zeta**2)/(1+xi/mp.mpf(140)**2))
    M=Mhol/(mp.sqrt(Om)*H)
    r=-xi/(zeta*(zeta*zeta+xi))
    et=elli[0]-1j*elli[1]
    ez=-1j*elli[2]/root
    Ot=Oi[0]-1j*Oi[1]
    Oz=-1j*Oi[2]/root
    B=M*(-mp.conj(W)/(P*Om**mp.mpf('1.5'))+aux[0]*(Ot/Om-2*et)
         +aux[1]*(Oz/Om-2*ez-mp.sqrt(2)*r/F))
    return dict(V=V,FT=aux[0],FZ=aux[1],FX=mp.sqrt(GT)*aux[0],
                FZc=mp.sqrt(hs/gs)*aux[1],B=B,d=d,s=abs(M)**2,
                Omega=Om,KTT=metric[0,0],determinant=mp.re(mp.det(metric)),
                matter_numerator=H, current=J)


rows=[]
for A in (mp.mpf(0),mp.mpf('1e-3')):
    for rho in (mp.mpf(0),mp.mpf('1e-6')):
        for U in (mp.mpf(0),mp.mpf('1e-8')):
            for phase in (mp.mpf(0),mp.mpf('.37')):
                P,mg,Fp,L,tau,sigma=map(mp.mpf,('1.7','.3','.07','.2','.003','1'))
                ref=reference(P,mg,A,rho)
                ts,GT=ref['tstar'],ref['GT']
                pars=dict(P=P,A=A,Fp=Fp,mg=mg,U=U,phase=phase,L=L,tau=tau,sigma=sigma,ref=ref)
                new=lambda tt,vv,yy:full(tt,vv,yy,**pars)
                old=lambda tt,vv,yy:full(tt,vv,yy,**pars,strength=0)
                axis=new(ts,0,0)
                delta=lambda tt,vv,yy:new(tt,vv,yy)['V']-old(tt,vv,yy)['V']
                f,u=abs(axis['FX'])**2,abs(axis['FZc'])**2
                C=axis['FX']*mp.conj(axis['FZc'])
                tag=f'A{A}_rho{rho}_U{U}_phase{phase}'
                ck(tag+'_U_normalization',u,U,scale=1)
                if U==0:
                    ck(tag+'_physical_mg_FX',f,3*P*mg*mg*ref['ps']**2/ref['Ds'])
                    ck(tag+'_rho_axis',axis['V'],rho,scale=1)
                    ck(tag+'_reference_stationarity_t',mp.diff(lambda tt:new(tt,0,0)['V'],ts),0,scale=1)
                for key in ('V','FT','FZ','B','s'):
                    ck(tag+'_'+key+'_axis',new(ts,0,0)[key]-old(ts,0,0)[key],0,scale=1)
                for orders in ((1,0,0),(0,1,0),(0,0,1)):
                    ck(tag+'_V_first_'+str(orders),mp.diff(delta,(ts,mp.mpf(0),mp.mpf(0)),orders),0,scale=1)
                for key in ('FT','FZ'):
                    for orders in ((1,0,0),(0,1,0),(0,0,1)):
                        fn=lambda tt,vv,yy:new(tt,vv,yy)[key]-old(tt,vv,yy)[key]
                        ck(tag+'_'+key+'_first_'+str(orders),mp.diff(fn,(ts,mp.mpf(0),mp.mpf(0)),orders),0,scale=1)
                hh=mp.matrix(3)
                conv=(mp.sqrt(2/GT),mp.sqrt(2/GT),1/Fp)
                for i in range(3):
                    for j in range(i,3):
                        orders=[0,0,0];orders[i]+=1;orders[j]+=1
                        hh[i,j]=hh[j,i]=mp.diff(delta,(ts,mp.mpf(0),mp.mpf(0)),tuple(orders))*conv[i]*conv[j]
                expect=mp.matrix([[(4*tau*tau*f+2*tau*u)/L**2,0,-4*tau*mp.im(C)/L**2],
                                  [0,(4*tau*tau*f+2*tau*u)/L**2,4*tau*mp.re(C)/L**2],
                                  [-4*tau*mp.im(C)/L**2,4*tau*mp.re(C)/L**2,(4*tau*f+12*u)/L**2]])
                for i in range(3):
                    for j in range(i,3):
                        ck(tag+f'_H_{i}{j}',hh[i,j],expect[i,j],scale=1)
                ck(tag+'_added_H_positive',int(min(mp.eigsy(hh,eigvals_only=True))>0),1)
                dd=new(ts,0,0)['d']-old(ts,0,0)['d']
                ck(tag+'_d_Q',dd,2*sigma*(u+tau*f)/L**2)
                ck(tag+'_persistent_mass_contact_correlation',dd-2*sigma*u/L**2,
                   sigma/2*(hh[2,2]-12*u/L**2))
                # Frozen radial extension: unlike the previous h(t) prefactor,
                # its mixed matter Hessian contains C'/C = -d ln(h_Q)/dt.
                h,hq,g,p,q=(ref[k] for k in ('hs','hqs','gs','ps','qs'))
                hp=-3*A/ts**4
                lq=-6*A/(ts**4*hq)
                Fhol=Fp*mp.sqrt(g/h)
                fty=mp.diff(lambda yy:new(ts,0,yy)['FT'],0)
                fty_pred=1j*Fhol*h/(mp.sqrt(2)*g*GT)*(2*hp/h+p/g)*axis['FZ']
                ck(tag+'_FT_y_exact',fty,fty_pred,scale=1)
                dy=mp.diff(lambda yy:new(ts,0,yy)['d']-old(ts,0,yy)['d'],0)
                dy_pred=2*mp.sqrt(2)*Fhol*(2*sigma/L**2)*h/g*((1+tau)*hp/h+lq-q/p+tau*p/(2*g))*mp.im(axis['FT']*mp.conj(axis['FZ']))
                ck(tag+'_soft_y_frozen_extension',dy,dy_pred,scale=1)
                def ddchi(xx):
                    return full(ts,0,0,**pars,chishift=xx)['d']-full(ts,0,0,**pars,strength=0,chishift=xx)['d']
                rw=Fp*mp.sqrt(U/2)/(P*mg)
                fchi=6*P*(p*p/ref['Ds'])*mg*mg*rw*(mp.cos(phase)-rw)
                ck(tag+'_soft_chi',mp.diff(ddchi,0),-4*sigma*U/L**2+2*sigma*tau*fchi/L**2,scale=1)
                # First B jet: auxiliary jets are unchanged, only matter current contributes.
                bdiff=lambda tt,vv,yy:new(tt,vv,yy)['B']-old(tt,vv,yy)['B']
                by=mp.diff(lambda yy:bdiff(ts,0,yy),0)
                bx=mp.diff(lambda tt:bdiff(tt,0,0),ts)*mp.sqrt(2/GT)
                bv=mp.diff(lambda vv:bdiff(ts,vv,0),0)*mp.sqrt(2/GT)
                ck(tag+'_B_y',by,-4j*mp.sqrt(2)*sigma*Fp*axis['FZc']/L**2,scale=1)
                ck(tag+'_B_x',bx,2*mp.sqrt(2)*sigma*tau*axis['FX']/L**2,scale=1)
                ck(tag+'_B_v',bv,-2j*mp.sqrt(2)*sigma*tau*axis['FX']/L**2,scale=1)
                ck(tag+'_metric_origin',int(axis['Omega']>0 and axis['KTT']>0 and axis['determinant']>0 and axis['matter_numerator']>0),1)
                for xx,vv,aa in ((mp.mpf('.001'),0,mp.mpf('.001')),(mp.mpf('.002'),mp.mpf('.001'),mp.mpf('.002'))):
                    pt=new(ts+mp.sqrt(2/GT)*xx,mp.sqrt(2/GT)*vv,aa/Fp)
                    ck(tag+f'_metric_neighbor_{xx}',int(pt['Omega']>0 and pt['KTT']>0 and pt['determinant']>0 and pt['matter_numerator']>0),1)
                rows.append(dict(A=str(A),rho=str(rho),U=str(U),phase=str(phase),tstar=mp.nstr(ts,30),
                                 delta_clock_mass=mp.nstr(hh[2,2],30),delta_radion_mass=mp.nstr(hh[0,0],30),
                                 delta_soft_Q=mp.nstr(dd,30),min_added_hessian_eigenvalue=mp.nstr(min(mp.eigsy(hh,eigvals_only=True)),30)))

# A physically normalized, nonzero-A counterexample: positive restoration and
# finite matter feedback are obtained together, without freezing F^T.
P=mp.mpf('2.435e27')**2; mg=mp.mpf('1e-6'); L=mp.mpf('1e12')
mKK=mp.mpf('1e14'); A=mp.zeta(3)/(48*mp.pi**4)*mKK*mKK/P
rho=mp.mpf('2.5181378331717977e-11'); Fp=mp.sqrt(mp.mpf('5.929225e50'))
ref=reference(P,mg,A,rho,hqfactor=1)
target=24*A*mg*mg
tau=2*target*L*L/(12*P*mg*mg)
pars=dict(P=P,A=A,Fp=Fp,mg=mg,U=mp.mpf(0),phase=mp.mpf(0),L=L,tau=tau,sigma=mp.mpf(1),ref=ref,hqfactor=1)
ts=ref['tstar'];new=lambda yy:full(ts,0,yy,**pars);old=lambda yy:full(ts,0,yy,**pars,strength=0)
delta_mass=mp.diff(lambda yy:new(yy)['V']-old(yy)['V'],0,2)/Fp**2
clock_mass=mp.diff(lambda yy:new(yy)['V'],0,2)/Fp**2
ck('physical_nonzero_A_persistent_mass',delta_mass,4*tau*abs(new(0)['FX'])**2/L**2)
ck('physical_nonzero_A_contact_relation',new(0)['d']-old(0)['d'],delta_mass/2)
ck('physical_nonzero_A_positive_clock_mass',int(clock_mass>0),1)
physical_hessians=[]
for uphys in (mp.mpf(0),mp.mpf('8.328200199084861e-16')):
    for theta in (mp.mpf(0),mp.pi/2):
        ps=dict(pars,U=uphys,phase=theta)
        fn=lambda tt,vv,yy:full(tt,vv,yy,**ps)['V']
        hessian=mp.matrix(3)
        conv=(mp.sqrt(2/ref['GT']),mp.sqrt(2/ref['GT']),1/Fp)
        for i in range(3):
            for j in range(i,3):
                orders=[0,0,0];orders[i]+=1;orders[j]+=1
                hessian[i,j]=hessian[j,i]=mp.diff(fn,(ts,mp.mpf(0),mp.mpf(0)),tuple(orders))*conv[i]*conv[j]
        eigenvalues=mp.eigsy(hessian,eigvals_only=True)
        ck(f'physical_full_local_H_positive_U{uphys}_phase{mp.nstr(theta,5)}',int(min(eigenvalues)>0),1)
        physical_hessians.append(dict(U=mp.nstr(uphys,25),phase=mp.nstr(theta,25),
                                      eigenvalues_eV2=[mp.nstr(x,35) for x in eigenvalues],
                                      matrix_eV2=[[mp.nstr(x,25) for x in row] for row in hessian.tolist()],
                                      scope='Coordinate local potential Hessian with reference canonical normalization; at U=0 the reference is a stationary point. At U>0 this is not the full covariant rolling perturbation operator.'))
physical=dict(A=mp.nstr(A,35),mg_eV=str(mg),mKK_eV=str(mKK),Lambda_eV=str(L),
              tau=mp.nstr(tau,35),tstar_minus_one=mp.nstr(ts-1,35),eta=mp.nstr(ref['eta'],35),
              target_negative_mass_squared=mp.nstr(-target,35),persistent_mass_squared=mp.nstr(delta_mass,35),
              full_local_clock_mass_squared=mp.nstr(clock_mass,35),delta_soft_Q_eV2=mp.nstr(new(0)['d']-old(0)['d'],35),
              canonical_FX_squared=mp.nstr(abs(new(0)['FX'])**2,35),full_local_hessians=physical_hessians,
              generic_tau_one_added_mass_eV2=mp.nstr(4*abs(new(0)['FX'])**2/L**2,35),
              generic_tau_one_soft_Q_eV2=mp.nstr(2*abs(new(0)['FX'])**2/L**2,35),
              scope='Local U=0 slice only; tau chosen algebraically to exceed the specified negative curvature by factor two; not a mechanism fixing tau.')

out=dict(scope='Independent full local K/W and charged quadratic spectrum for an explicitly declared centered four-dimensional current; no UV or Riemann derivation, loops, or cosmological trajectory',
         precision_decimal_digits=mp.mp.dps,check_count=len(checks),passed=sum(c['passed'] for c in checks),all_pass=all(c['passed'] for c in checks),
         checks=checks,local_points=rows,physical_counterexample=physical,code_sha256=hashlib.sha256(HERE.read_bytes()).hexdigest(),
         execution_notes='Initial mass/current execution passed 507/507; extended with exact soft first jets and four full local Hessian examples. All fixed tolerances retained, no failed checks. Dimensionless amplified points use hQ=1+2A/t^3 to test separate metrics; the physical example uses hQ=hZ=1+A/t^3.',
         input_read_scope=['experiments/clock_stabilization_v01/protocol.md','reports/clock_quartic_geometry_20260926.md'],
         external_sources='No new external scientific claim used; all identities explicitly derived from the declared local K/W.')
HERE.with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ('check_count','passed','all_pass','physical_counterexample')},indent=2))
if not out['all_pass']:
    print(json.dumps([c for c in checks if not c['passed']],indent=2))
    raise SystemExit(1)
