#!/usr/bin/env python3
"""Independent charged Gaussian determinant of the frozen hidden-current EFT.

Builds the full two-chiral Kähler metric and matter-metric jets directly.
No experiment or other review implementation is imported. This is the selected
charged determinant, not the heavy-vector, graviton or full SUGRA determinant.
"""
from pathlib import Path
import hashlib
import json
import mpmath as mp

HERE = Path(__file__).resolve()
REPO = HERE.parent.parent
SOURCE = REPO / 'experiments/clock_stabilization_v01/results/inputs.json'
INPUT = json.loads(SOURCE.read_text())
PARENT = INPUT['parent_parameters']
CHECKS = []


def check(name, actual, expected, tol='1e-70', scale=None):
    den = max(abs(expected), mp.mpf('1e-100')) if scale is None else mp.mpf(scale)
    err = abs(actual-expected)/den
    CHECKS.append(dict(name=name, actual=mp.nstr(actual, 35), expected=mp.nstr(expected, 35),
                       relative_error=mp.nstr(err, 15), tolerance=tol,
                       passed=bool(err < mp.mpf(tol))))


def reference(dps):
    mp.mp.dps = dps
    P = mp.mpf('2.435e27')**2
    mg, KK, lam = mp.mpf('1e-6'), mp.mpf('1e14'), mp.mpf('1e12')
    A = mp.zeta(3)/(48*mp.pi**4)*KK**2/P
    rho = mp.mpf(str(PARENT['rho']))
    rr = rho/(3*P*mg**2)
    def jets(t, eta):
        dt = t-1
        return (t+eta*dt**2+dt**4-A/t**2,
                1+2*eta*dt+4*dt**3+2*A/t**3,
                2*eta+12*dt**2-6*A/t**4,
                24*dt+24*A/t**5)
    def equations(t, eta):
        g, p, q, qt = jets(t, eta)
        return ((q-rr*p*p/(g*(1+rr)))/A,
                (qt-2*p*q/g+q*q/p)/A)
    t, eta = mp.findroot(equations, (1-A+rho/(36*P*mg**2),
                                    3*A+rho/(6*P*mg**2)),
                        tol=mp.mpf(10)**(-dps+25))
    g, p, q, qt = jets(t, eta)
    h, hp, hpp = 1+A/t**3, -3*A/t**4, 12*A/t**5
    F2 = mp.mpf(str(PARENT['F_squared_eV2']))
    F = mp.sqrt(F2*g/h)
    GTT = 3*P*(p*p-g*q)/(g*g)
    return dict(P=P, mg=mg, A=A, rho=rho, lam=lam, t=t, eta=eta,
                g=g, p=p, q=q, h=h, hp=hp, hpp=hpp, F=F, F2=F2,
                GTT=GTT, Xi=mp.mpf(str(PARENT['xi'])),
                chii=mp.mpf(str(INPUT['chi_first'])), Mi=mp.mpf('1e11'),
                Ui=mp.mpf(str(INPUT['Umax_eV4'])))


def spectrum(chi, y, rH, phase, clock_on, p, sigma=1, contacts=True):
    """Fixed T*=reference; frozen canonical current continued off the clock axis.

    Xc=sqrt(GTT*) (T-T*), kc=(h*/g*) k, Sc=(hQ*/g*) S.
    Ω=g-(h k+hQ S)/(3P)+g*/(3PΛ²)(kc+σSc+τ|Xc|²)^2.
    All starred coefficients are frozen under T/Z differentiation.
    """
    P, g, gt, gtt, h, hp, hpp, F, L = [p[z] for z in
        ('P','g','p','q','h','hp','hpp','F','lam')]
    tau = rH*p['A']*L**2/P
    # No-contact control removes the complete square, not just its Q couplings.
    a = 1/L**2 if contacts else mp.mpf(0)
    kcscale = h/g
    k, ky, kyy = F**2*y*y, 2*F**2*y, 2*F**2
    kc, kcy, kcyy = kcscale*k, kcscale*ky, kcscale*kyy
    O = g-h*k/(3*P)+g*a*kc*kc/(3*P)
    Ot = gt-hp*k/(3*P)
    Ott = gtt-hpp*k/(3*P)+2*g*a*tau*kc*p['GTT']/(3*P)
    Oy = -h*ky/(3*P)+2*g*a*kc*kcy/(3*P)
    Oyy = -h*kyy/(3*P)+2*g*a*(kcy*kcy+kc*kcyy)/(3*P)
    Oty = -hp*ky/(3*P)
    Kt, Ky = -3*P*Ot/O, -3*P*Oy/O
    Ktt = -3*P*(Ott/O-Ot*Ot/O**2)
    Kyy = -3*P*(Oyy/O-Oy*Oy/O**2)
    Kty = -3*P*(Oty/O-Ot*Oy/O**2)
    rt = mp.sqrt(2)*F
    G = mp.matrix([[Ktt,1j*Kty/rt],[-1j*Kty/rt,Kyy/rt**2]])
    U = p['Ui']*mp.exp(-2*(chi-p['chii'])) if clock_on else mp.mpf(0)
    wr = -g**mp.mpf('1.5')*mp.sqrt(p['F2']*U/2)
    Wc = P*g**mp.mpf('1.5')*p['mg']*mp.exp(1j*phase)
    W, Wz = Wc+wr*mp.exp(-1j*y), -mp.sqrt(2)*wr/F*mp.exp(-1j*y)
    DW = mp.matrix([Kt*W/P, Wz-1j*Ky*W/(rt*P)])
    inv = G**-1
    aux = -O**mp.mpf('-1.5')*inv.T*DW.conjugate()
    V = mp.re(O**-3*((DW.conjugate().T*inv*DW)[0]-3*abs(W)**2/P))
    # H_Q=Ω Z_Q, with derivatives of the complete neutral-current dependence.
    H = h-2*sigma*h*a*kc
    Ht = hp
    Htt = hpp-2*sigma*h*a*tau*p['GTT']
    kcz, kczb = -1j*kcscale*mp.sqrt(2)*F*y, 1j*kcscale*mp.sqrt(2)*F*y
    Hz, Hzb = -2*sigma*h*a*kcz, -2*sigma*h*a*kczb
    Hzz = -2*sigma*h*a*kcscale
    ellt, ellz = Ht/H, Hz/H
    elltt = Htt/H-Ht*Ht/H**2
    elltz = -Ht*Hzb/H**2
    ellzt = -Hz*Ht/H**2
    ellzz = Hzz/H-Hz*Hzb/H**2
    d = 2*V/(3*P)-mp.re(abs(aux[0])**2*elltt+abs(aux[1])**2*ellzz
            +aux[0]*mp.conj(aux[1])*elltz+aux[1]*mp.conj(aux[0])*ellzt)
    z = chi+1j*y
    raw = mp.sqrt((1+p['Xi']/z**2)/(1+p['Xi']/p['chii']**2))
    Mhol = mp.sqrt(g)*h*p['Mi']*raw
    mass = Mhol/(mp.sqrt(O)*H)
    r = -p['Xi']/(z*(z*z+p['Xi']))
    B = mass*(-mp.conj(W)/(P*O**mp.mpf('1.5'))
              +aux[0]*(Ot/O-2*ellt)
              +aux[1]*(-1j*Oy/(rt*O)-2*ellz-mp.sqrt(2)*r/F))
    return dict(s=abs(mass)**2, d=d, B=B, V=V, FT=aux[0], FZ=aux[1],
                FX2=p['GTT']*abs(aux[0])**2, U=U, Omega=O,
                GTT=G[0,0], determinant=mp.re(mp.det(G)), matter_H=H)


def cw(sp, p):
    s, d, b = sp['s'], sp['d'], abs(sp['B'])
    def f(x):
        return x*x*(mp.log(x/p['Mi']**2)-mp.mpf('1.5'))
    return (f(s+d+b)+f(s+d-b)-2*f(s))/(32*mp.pi**2)


def cw_chain(chi, rH, phase, on, p, contacts=True):
    axis = spectrum(chi, 0, rH, phase, on, p, contacts=contacts)
    values = {}
    def eig(c, y):
        z = spectrum(c,y,rH,phase,on,p,contacts=contacts)
        return [z['s']+z['d']+abs(z['B']), z['s']+z['d']-abs(z['B']), z['s']]
    lam = eig(chi, 0)
    def fp(x): return 2*x*(mp.log(x/p['Mi']**2)-1)
    def fpp(x): return 2*mp.log(x/p['Mi']**2)
    for var, order in (('chi',1),('y',1),('y',2)):
        js = []
        for k in range(3):
            fn = (lambda z, k=k:eig(z,0)[k]) if var=='chi' else (lambda z,k=k:eig(chi,z)[k])
            x0 = chi if var=='chi' else 0
            first = mp.diff(fn,x0)
            val = fp(lam[k])*first if order==1 else (
                fpp(lam[k])*first**2+fp(lam[k])*mp.diff(fn,x0,2))
            js.append(val)
        values[var+str(order)] = (js[0]+js[1]-2*js[2])/(32*mp.pi**2)
    values['value'] = cw(axis,p)
    values['my2'] = mp.diff(lambda y:spectrum(chi,y,rH,phase,on,p,contacts=contacts)['V'],0,2)/p['F2']
    values['loop_my2'] = values['y2']/p['F2']
    values['ddyy'] = mp.diff(lambda y:spectrum(chi,y,rH,phase,on,p,contacts=contacts)['d'],0,2)
    values['B_y'] = mp.diff(lambda y:spectrum(chi,y,rH,phase,on,p,contacts=contacts)['B'],0)
    values['B_yy'] = mp.diff(lambda y:spectrum(chi,y,rH,phase,on,p,contacts=contacts)['B'],0,2)
    values['s_yy'] = mp.diff(lambda y:spectrum(chi,y,rH,phase,on,p,contacts=contacts)['s'],0,2)
    values['linear_d_cw_chi'] = mp.diff(lambda c:(lambda z:z['s']*z['d']*(mp.log(z['s']/p['Mi']**2)-1)/(8*mp.pi**2))(
            spectrum(c,0,rH,phase,on,p,contacts=contacts)),chi)
    values['B2_cw_chi'] = mp.diff(lambda c:(lambda z:abs(z['B'])**2*mp.log(z['s']/p['Mi']**2)/(16*mp.pi**2))(
            spectrum(c,0,rH,phase,on,p,contacts=contacts)),chi)
    return axis, values


def run(dps):
    p = reference(dps)
    rows = []
    for label, chi, on, phase in [('first',p['chii'],True,mp.mpf(0)),
                           ('last',mp.mpf(str(INPUT['chi_last'])),True,mp.mpf(0)),
                           ('off_last',mp.mpf(str(INPUT['chi_last'])),False,mp.mpf(0)),
                           ('last_quadrature',mp.mpf(str(INPUT['chi_last'])),True,mp.pi/2)]:
        for rH in (0,2,4):
            axis, val = cw_chain(chi,rH,phase,on,p)
            old = spectrum(chi,0,rH,phase,on,p,contacts=False)
            tau = rH*p['A']*p['lam']**2/p['P']
            tag = f'dps{dps}_{label}_rH{rH}'
            expected_d = 2*(axis['U']+tau*axis['FX2'])/p['lam']**2
            check(tag+'_shared_axis_d',axis['d']-old['d'],expected_d)
            check(tag+'_B_axis_unchanged',axis['B'],old['B'])
            check(tag+'_axis_V_unchanged',axis['V'],old['V'])
            oldmy = mp.diff(lambda y:spectrum(chi,y,rH,phase,on,p,contacts=False)['V'],0,2)/p['F2']
            expected_my = 12*axis['U']/p['lam']**2+4*tau*axis['FX2']/p['lam']**2
            check(tag+'_full_mass_shift',val['my2']-oldmy,expected_my)
            delta_dy = mp.diff(lambda y:spectrum(chi,y,rH,phase,on,p)['d']
                            -spectrum(chi,y,rH,phase,on,p,contacts=False)['d'],0)
            dy_formula = (4*mp.sqrt(2)*p['F']*(p['h']/p['g'])/p['lam']**2
                *((2+tau)*p['hp']/p['h']-p['q']/p['p']+tau*p['p']/(2*p['g']))
                *mp.im(axis['FT']*mp.conj(axis['FZ'])))
            check(tag+'_frozen_shared_first_y',delta_dy,dy_formula,scale=max(abs(dy_formula),mp.mpf('1e-80')))
            by_old = mp.diff(lambda y:spectrum(chi,y,rH,phase,on,p,contacts=False)['B'],0)
            by_formula = 4j*mp.sqrt(2)*mp.sqrt(axis['s']*p['F2']*axis['U'])/p['lam']**2
            check(tag+'_shared_B_first_y',val['B_y']-by_old,by_formula,
                    scale=max(abs(by_formula),mp.mpf('1e-80')))
            check(tag+'_kinetic_positive',int(axis['Omega']>0 and axis['GTT']>0
                        and axis['determinant']>0 and axis['matter_H']>0),1)
            check(tag+'_charged_positive',int(axis['s']+axis['d']>abs(axis['B'])),1)
            check(tag+'_CW_expansion_first_derivative',val['chi1'],val['linear_d_cw_chi']+val['B2_cw_chi'],tol='1e-35')
            # Independent direct differentiation of the unexpanded determinant.
            direct = mp.diff(lambda c:cw(spectrum(c,0,rH,phase,on,p),p),chi)
            check(tag+'_CW_chain_vs_direct_chi',val['chi1'],direct,tol='1e-65')
            directyy = mp.diff(lambda y:cw(spectrum(chi,y,rH,phase,on,p),p),0,2)
            check(tag+'_CW_chain_vs_direct_yy',val['y2'],directyy,tol='1e-65')
            leading = mp.diff(lambda c:(lambda z:
                z['s']*(2*(z['U']+tau*z['FX2'])/p['lam']**2)
                *(mp.log(z['s']/p['Mi']**2)-1)/(8*mp.pi**2))(
                    spectrum(c,0,rH,phase,on,p)),chi)
            old_cw_chi = mp.diff(lambda c:cw(spectrum(c,0,rH,phase,on,p,contacts=False),p),chi)
            check(tag+'_shared_CW_chi_vs_linear',val['chi1']-old_cw_chi,leading,
                  tol='1e-35',scale=max(abs(leading),mp.mpf('1e-80')))
            denom = 2*axis['U'] if on else mp.mpf(1)
            row = dict(precision=dps, clock_label=label, rH=rH, phase=phase,
                chi=chi, U=axis['U'], tau=tau, A=p['A'], FX2=axis['FX2'],
                s=axis['s'], d=axis['d'], B=axis['B'],
                scalar_splitting_relative=(abs(axis['d'])+abs(axis['B']))/axis['s'],
                charged_CW_value=val['value'], charged_CW_chi=val['chi1'],
                charged_CW_y=val['y1'], charged_CW_yy=val['y2'],
                selected_CW_force_ratio=abs(val['chi1'])/denom if on else None,
                selected_CW_full_gradient_ratio=mp.sqrt(val['chi1']**2+val['y1']**2)/denom if on else None,
                shared_d_linear_CW_chi=leading,
                shared_d_y=delta_dy, full_CW_chi_without_square=old_cw_chi,
                tree_my2=val['my2'], selected_CW_my2=val['loop_my2'],
                selected_loop_to_tree=val['loop_my2']/val['my2'],
                d_yy=val['ddyy'], B_y=val['B_y'], B_yy=val['B_yy'], s_yy=val['s_yy'])
            rows.append({k:(mp.nstr(v,90) if isinstance(v,(mp.mpf,mp.mpc)) else v) for k,v in row.items()})
            print(f'completed {tag}',flush=True)
    return rows


def tail_diagnostic(dps):
    p=reference(dps)
    def ratio(chi):
        force=mp.diff(lambda c:cw(spectrum(c,0,4,0,True,p),p),chi)
        return abs(force)/(2*spectrum(chi,0,4,0,True,p)['U'])
    crossing=mp.findroot(lambda c:mp.log(ratio(c)/mp.mpf('.1')),(150,152),
                         tol=mp.mpf(10)**(-dps+50))
    check(f'tail_dps{dps}_crossing',ratio(crossing),mp.mpf('.1'),tol='1e-90')
    check(f'tail_dps{dps}_positive_tree_at_crossing',int(mp.diff(
        lambda y:spectrum(crossing,y,4,0,True,p)['V'],0,2)>0),1)
    return dict(precision=dps,rH=4,phase=0,
                crossing_chi=mp.nstr(crossing,100),
                gap_after_old_last=mp.nstr(crossing-mp.mpf(str(INPUT['chi_last'])),100),
                ratio_chi_145=mp.nstr(ratio(mp.mpf(145)),80),
                ratio_chi_150=mp.nstr(ratio(mp.mpf(150)),80),
                ratio_chi_155=mp.nstr(ratio(mp.mpf(155)),80),
                scope='Fixed T* and same prescribed U(chi); diagnostic coordinate extension only, no new cosmic trajectory or epoch prediction.')


def main():
    first, second = run(160), run(220)
    outputs = first+second
    for a,b in zip(first,second):
        for key in ('charged_CW_value','charged_CW_chi','charged_CW_yy','tree_my2','selected_CW_my2','d_yy','s_yy'):
            check(f'precision_{a["clock_label"]}_rH{a["rH"]}_{key}',mp.mpf(a[key]),mp.mpf(b[key]),tol='1e-60')
    tails=[tail_diagnostic(160),tail_diagnostic(220)]
    check('tail_precision_crossing',mp.mpf(tails[0]['crossing_chi']),
          mp.mpf(tails[1]['crossing_chi']),tol='1e-90')
    out = dict(scope='Complete selected charged Gaussian determinant at fixed matched T*, including all spectrum y second jets; not full UV/SUGRA loops or cosmological evolution',
               precisions=[160,220],check_count=len(CHECKS),passed=sum(z['passed'] for z in CHECKS),
               all_pass=all(z['passed'] for z in CHECKS),checks=CHECKS,points=outputs,tail_diagnostics=tails,
               implementation='Independent full Kähler matrix and ell_Q Hessian. Frozen canonically normalized current extension. No experiment implementation imported.',
               code_sha256=hashlib.sha256(HERE.read_bytes()).hexdigest(),
               input_path=str(SOURCE.relative_to(REPO)),input_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
               numerical_notes='Unexpanded scalar/fermion determinant and eigenvalue-chain derivatives; soft expansion only a cross-check. All failures retained in output.',
               finite_convention='mu=Mi, scalar/fermion subtraction constant 3/2; unknown finite local counterterms and all other multiplets excluded.')
    HERE.with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ('check_count','passed','all_pass')},indent=2))
    if not out['all_pass']:
        print(json.dumps([z for z in CHECKS if not z['passed']],indent=2))
        raise SystemExit(1)


if __name__=='__main__':
    main()
