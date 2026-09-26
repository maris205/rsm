#!/usr/bin/env python3
"""Independent full two-chiral K/W checks of a local quartic clock operator.

No project experiment code is imported. Dimensionless amplified points test
identities; all ordinary and charged quantities follow from Omega derivatives.
"""
from pathlib import Path
import hashlib
import json
import mpmath as mp

mp.mp.dps = 100
HERE = Path(__file__).resolve()
checks = []


def ck(name, actual, expected, scale=None, tol='1e-70'):
    den = max(abs(expected), mp.mpf('1e-45')) if scale is None else mp.mpf(scale)
    err = abs(actual-expected)/den
    checks.append(dict(name=name, error=mp.nstr(err, 15), tolerance=tol,
                       passed=bool(err < mp.mpf(tol))))


def model(y, *, P, g, p, q, h, hp, hpp, hq, lq1, lq2,
          F, Wc, w, lam, chi=mp.mpf(140), xi=mp.mpf(19000),
          sigma=mp.mpf(0), gamma=None):
    """Full exact K matrix at fixed T; g,p,q,h jets define this local slice."""
    a = mp.mpf(0) if lam is None else 1/lam**2
    k = F**2*y**2-a*F**4*y**4
    ky = 2*F**2*y-4*a*F**4*y**3
    kyy = 2*F**2-12*a*F**4*y**2
    Om = g-h*k/(3*P)
    Ot = p-hp*k/(3*P)
    Ott = q-hpp*k/(3*P)
    Oy = -h*ky/(3*P)
    Oyy = -h*kyy/(3*P)
    Oty = -hp*ky/(3*P)
    Kt, Ky = -3*P*Ot/Om, -3*P*Oy/Om
    Ktt = -3*P*(Ott/Om-Ot**2/Om**2)
    Kyy = -3*P*(Oyy/Om-Oy**2/Om**2)
    Kty = -3*P*(Oty/Om-Ot*Oy/Om**2)
    root = mp.sqrt(2)*F
    metric = mp.matrix([[Ktt, 1j*Kty/root], [-1j*Kty/root, Kyy/root**2]])
    W = Wc+w*mp.exp(-1j*y)
    Wz = -mp.sqrt(2)*w/F*mp.exp(-1j*y)
    DW = mp.matrix([Kt*W/P, Wz-1j*Ky*W/(root*P)])
    aux = -Om**mp.mpf('-1.5')*(metric**-1).T*DW.conjugate()
    V = mp.re(Om**-3*((DW.conjugate().T*(metric**-1)*DW)[0]-3*abs(W)**2/P))
    mbar = mp.conj(W)/(P*Om**mp.mpf('1.5'))
    Mhol = mp.sqrt(1+xi/(chi+1j*y)**2)
    r = -xi/((chi+1j*y)*((chi+1j*y)**2+xi))
    # For the square [k + sigma*gamma*S]^2, its Q-quadratic contact
    # uses the original quadratic k0 even though Q=0 Omega uses k_eff.
    gamma = hq/h if gamma is None else gamma
    C = 2*sigma*gamma*h/hq*a
    Cp = C*(hp/h-lq1)
    Cpp = C*((hpp/h-(hp/h)**2)-lq2+(hp/h-lq1)**2)
    k0 = F*F*y*y
    k0z, k0zb = -1j*mp.sqrt(2)*F*y, 1j*mp.sqrt(2)*F*y
    normal = 1-C*k0
    contact_tt = -Cpp*k0/normal-Cp**2*k0**2/normal**2
    contact_tzb = -Cp*k0zb/normal**2
    contact_ztb = -Cp*k0z/normal**2
    contact_zz = -C/normal-C*C*k0z*k0zb/normal**2
    ellt = lq1-Cp*k0/normal
    ellz = -C*k0z/normal
    M = Mhol/(mp.sqrt(Om)*hq*normal)
    Oz = -1j*Oy/root
    B = M*(-mbar+aux[0]*(Ot/Om-2*ellt)+aux[1]*(Oz/Om-2*ellz-mp.sqrt(2)*r/F))
    d = 2*V/(3*P)-lq2*abs(aux[0])**2-mp.re(abs(aux[0])**2*contact_tt+abs(aux[1])**2*contact_zz
                                         +aux[0]*mp.conj(aux[1])*contact_tzb
                                         +aux[1]*mp.conj(aux[0])*contact_ztb)
    det_formula = 9*P**2/(2*F**2*Om**3)*(Om*(Ott*Oyy-Oty**2)-Ott*Oy**2-Oyy*Ot**2+2*Oty*Ot*Oy)
    return dict(V=V, FT=aux[0], FZ=aux[1], B=B, d=d, s=abs(M)**2,
                Omega=Om, KTT=Ktt, determinant=mp.re(mp.det(metric)),
                determinant_formula=det_formula, matter_normal=normal)


for A in (mp.mpf('1e-3'), mp.mpf('1e-8')):
    t = mp.mpf('1.03')
    h = 1+A/t**3
    hp, hpp = -3*A/t**4, 12*A/t**5
    hq = 1+2*A/t**3
    hqp, hqpp = -6*A/t**4, 24*A/t**5
    for q in (mp.mpf(0), mp.mpf('.004')):
        for phase in (mp.mpf(0), mp.pi/2, mp.mpf('.37')):
            for lam in (mp.mpf('.09'), mp.mpf('.23')):
                pars = dict(P=mp.mpf('1.7'), g=mp.mpf('1.1'), p=mp.mpf('.97'), q=q,
                            h=h, hp=hp, hpp=hpp, hq=hq, lq1=hqp/hq,
                            lq2=hqpp/hq-(hqp/hq)**2, F=mp.mpf('.07'),
                            Wc=mp.mpf('.5')*mp.exp(1j*phase), w=mp.mpf('-1e-5'))
                new = lambda yy: model(yy, **pars, lam=lam)
                old = lambda yy: model(yy, **pars, lam=None)
                tag = f'A{A}_q{q}_phase{mp.nstr(phase,5)}_L{lam}'
                P, g, p, F = (pars[x] for x in ('P', 'g', 'p', 'F'))
                w, Wc = pars['w'], pars['Wc']
                D = p*p-g*q
                W = Wc+w
                Uhol = 2*w*w/F**2
                Uphys = Uhol/(g*g*h)
                Fphys2 = F*F*h/g
                Lphys2 = h/g*lam**2
                ck(tag+'_axis_V', new(0)['V'], Uphys+3*q*abs(W)**2/(P*g*g*D))
                ck(tag+'_axis_FT', new(0)['FT'], p*mp.conj(W)/(P*mp.sqrt(g)*D))
                ck(tag+'_axis_FZ', new(0)['FZ'], mp.sqrt(2)*w/(F*mp.sqrt(g)*h))
                for key in ('V', 'FT', 'FZ', 'B', 'd', 's'):
                    for n in (0, 1):
                        actual = mp.diff(lambda yy: new(yy)[key]-old(yy)[key], 0, n)
                        ck(tag+f'_{key}_unchanged_jet{n}', actual, 0, scale=1)
                dv2 = mp.diff(lambda yy: new(yy)['V']-old(yy)['V'], 0, 2)
                ck(tag+'_additional_physical_mass', dv2/Fphys2, 12*Uphys/Lphys2)
                dfz2 = mp.diff(lambda yy: new(yy)['FZ']-old(yy)['FZ'], 0, 2)
                ck(tag+'_FZ_second_jet', dfz2, 12*F*F/lam**2*old(0)['FZ'])
                dft2 = mp.diff(lambda yy: new(yy)['FT']-old(yy)['FT'], 0, 2)
                ck(tag+'_FT_second_jet', dft2, 0, scale=1)
                dd2 = mp.diff(lambda yy: new(yy)['d']-old(yy)['d'], 0, 2)
                ck(tag+'_charged_diagonal_second_jet', dd2, 2*dv2/(3*P))
                r = -mp.mpf(19000)/(mp.mpf(140)*(mp.mpf(140)**2+mp.mpf(19000)))
                Mphys = mp.sqrt(1+mp.mpf(19000)/mp.mpf(140)**2)/(mp.sqrt(g)*hq)
                Bclock = -Mphys*old(0)['FZ']*mp.sqrt(2)*r/F
                db2 = mp.diff(lambda yy: new(yy)['B']-old(yy)['B'], 0, 2)
                ck(tag+'_charged_B_second_jet', db2, 12*F*F/lam**2*Bclock)
                ck(tag+'_charged_mass_second_jet', mp.diff(lambda yy: new(yy)['s']-old(yy)['s'], 0, 2), 0, scale=1)
                for frac in (mp.mpf('.01'), mp.mpf('.10'), mp.mpf('.15')):
                    yy = mp.sqrt(frac)*lam/F
                    point = new(yy)
                    ck(tag+f'_determinant_frac{frac}', point['determinant'], point['determinant_formula'])
                    ck(tag+f'_metric_positive_frac{frac}', int(point['Omega']>0 and point['KTT']>0 and point['determinant']>0), 1)

# Exact uniform-geometry range: determinant changes sign at k/L^2 = 1/6.
uniform = dict(P=mp.mpf(1), g=mp.mpf('1.2'), p=mp.mpf('.9'), q=mp.mpf(0),
               h=mp.mpf('1.1'), hp=mp.mpf(0), hpp=mp.mpf(0), hq=mp.mpf(1),
               lq1=mp.mpf(0), lq2=mp.mpf(0), F=mp.mpf('.1'), Wc=mp.mpf('.5'),
               w=mp.mpf('-1e-5'), lam=mp.mpf('.2'))
for frac in (mp.mpf('.01'), mp.mpf('.1'), mp.mpf('.16'), mp.mpf('.17'), mp.mpf('.2')):
    yy = mp.sqrt(frac)*uniform['lam']/uniform['F']
    point = model(yy, **uniform)
    exact = 3*uniform['P']*uniform['h']*uniform['p']**2/point['Omega']**3*(1-6*frac)
    ck(f'uniform_frac{frac}_determinant', point['determinant'], exact)
    ck(f'uniform_frac{frac}_sign', int(point['determinant']>0), int(frac<mp.mpf(1)/6))

# Generic-phase local valleys: nonconstant finite correction is retained.
valleys = []
A = mp.mpf('1e-5'); g=mp.mpf(1); p=mp.mpf(1); h=1+A; hp=-3*A; hpp=12*A
Fphys = mp.mpf('.1'); mg=mp.mpf('.3'); Lphys=mp.mpf('.002'); theta=mp.pi/2
Fhol=Fphys/mp.sqrt(h); Lhol=Lphys/mp.sqrt(h)
for U in (mp.mpf('1e-8'), mp.mpf('3e-9'), mp.mpf('1e-9')):
    w=-Fphys*mp.sqrt(U/2)
    pars=dict(P=mp.mpf(1), g=g, p=p, q=mp.mpf(0), h=h, hp=hp, hpp=hpp,
              hq=mp.mpf(1), lq1=mp.mpf(0), lq2=mp.mpf(0), F=Fhol,
              Wc=mg*mp.exp(1j*theta), w=w, lam=Lhol)
    f=lambda yy:model(yy,**pars)['V']
    J=mp.diff(f,0); H=mp.diff(f,0,2)
    yquad=-J/H
    ys=mp.findroot(lambda yy:mp.diff(f,yy),(yquad*mp.mpf('.99'),yquad*mp.mpf('1.01')),tol=mp.mpf('1e-85'))
    point=model(ys,**pars)
    ck(f'valley_U{U}_stationarity',mp.diff(f,ys),0,scale=1,tol='1e-65')
    ck(f'valley_U{U}_positive_metric',int(point['Omega']>0 and point['KTT']>0 and point['determinant']>0),1)
    ck(f'valley_U{U}_positive_local_curvature',int(mp.diff(f,ys,2)>0),1)
    const_leading=-3*A*A*mg*mg*Lphys*Lphys
    vals=dict(U=mp.nstr(U,25), J=mp.nstr(J,30), H=mp.nstr(H,30), y_quadratic=mp.nstr(yquad,30),
              y_exact=mp.nstr(ys,30), delta_V_exact=mp.nstr(f(ys)-f(0),30),
              delta_V_quadratic=mp.nstr(-J*J/(2*H),30),
              constant_leading=mp.nstr(const_leading,30),
              relative_to_leading=mp.nstr((f(ys)-f(0))/const_leading-1,30),
              metric_fraction=mp.nstr(Fhol**2*ys**2/Lhol**2,30))
    valleys.append(vals)

# Shared-current contacts: full ell_Q Hessian, including mixed t/Z entries.
for A in (mp.mpf('1e-3'), mp.mpf('1e-8')):
    t=mp.mpf('1.03'); h=1+A/t**3; hp=-3*A/t**4; hpp=12*A/t**5
    hq=1+2*A/t**3; hqp=-6*A/t**4; hqpp=24*A/t**5
    for q in (mp.mpf(0),mp.mpf('.004')):
        for phase in (mp.pi/2,mp.mpf('.37')):
            for sigma in (-1,1):
                pars=dict(P=mp.mpf('1.7'),g=mp.mpf('1.1'),p=mp.mpf('.97'),q=q,h=h,hp=hp,hpp=hpp,
                          hq=hq,lq1=hqp/hq,lq2=hqpp/hq-(hqp/hq)**2,F=mp.mpf('.07'),
                          Wc=mp.mpf('.5')*mp.exp(1j*phase),w=mp.mpf('-1e-5'),lam=mp.mpf('.09'),gamma=hq/h)
                new=lambda yy:model(yy,**pars,sigma=sigma)
                old=lambda yy:model(yy,**pars,sigma=0)
                g,p,F,L=(pars[x] for x in ('g','p','F','lam'))
                C=2*sigma/L**2
                FT,FZ=old(0)['FT'],old(0)['FZ']
                Uphys=2*pars['w']**2/(F**2*g*g*h)
                Fphys=F*mp.sqrt(h/g);Lphys2=h/g*L*L
                tag=f'contact_A{A}_q{q}_phase{mp.nstr(phase,5)}_sig{sigma}'
                ck(tag+'_FZ_first_jet',mp.diff(lambda yy:old(yy)['FZ'],0),1j*FZ-1j*mp.sqrt(2)*F*(hp/h-q/p)*FT)
                ck(tag+'_d_axis',new(0)['d']-old(0)['d'],2*sigma*Uphys/Lphys2)
                dy=mp.diff(lambda yy:new(yy)['d']-old(yy)['d'],0)
                ck(tag+'_d_y_mixed_cancellation',dy,2*mp.sqrt(2)*F*C*(hqp/hq-q/p)*mp.im(FT*mp.conj(FZ)))
                def dchi(x):
                    pp=dict(pars,w=pars['w']*mp.exp(-x),chi=mp.mpf(140)+x)
                    return model(0,**pp,sigma=sigma)['d']-model(0,**pp,sigma=0)['d']
                ck(tag+'_d_chi',mp.diff(dchi,0),-4*sigma*Uphys/Lphys2)
                ck(tag+'_B_axis',new(0)['B']-old(0)['B'],0,scale=1)
                Mphys=mp.sqrt(1+mp.mpf(19000)/mp.mpf(140)**2)/(mp.sqrt(g)*hq)
                by=mp.diff(lambda yy:new(yy)['B']-old(yy)['B'],0)
                ck(tag+'_B_y',by,4j*mp.sqrt(2)*sigma*Mphys*Fphys*mp.sqrt(Uphys)/Lphys2)
                for key in ('s','V'):
                    for order in (0,1):
                        ck(tag+f'_{key}_unchanged_jet{order}',mp.diff(lambda yy:new(yy)[key]-old(yy)[key],0,order),0,scale=1)
                atedge=new(mp.sqrt(mp.mpf('.15'))*L/F)
                ck(tag+'_matter_metric_positive',int(atedge['matter_normal']>0),1)

out=dict(scope='Independent exact local K/W and charged-spectrum jets of an explicitly added 4D quartic clock operator; no UV completion, complete loop matching or cosmic trajectory',
         precision_decimal_digits=mp.mp.dps,check_count=len(checks),passed=sum(c['passed'] for c in checks),
         all_pass=all(c['passed'] for c in checks),checks=checks,generic_phase_valleys=valleys,
         code_sha256=hashlib.sha256(HERE.read_bytes()).hexdigest(),
         input_read_scope=['experiments/sequestering_test_v01/protocol.md','reports/protection_power_correction_20260925.md'],
         external_sources='No new external scientific claim used: full K/W identities derived explicitly from the declared local model.',
         execution_notes='Initial quartic-only run passed 667/667 checks; scope then extended to shared-current contact jets. Dimensionless amplified values intentionally expose corrections. Same fixed tolerances for all checks; no observational statistics.')
HERE.with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ('check_count','passed','all_pass','generic_phase_valleys')},indent=2))
if not out['all_pass']:
    print(json.dumps([c for c in checks if not c['passed']],indent=2))
    raise SystemExit(1)
