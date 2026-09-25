#!/usr/bin/env python3
"""Independent full K-matrix checks; does not import experiment implementations."""
from pathlib import Path
import hashlib
import json
import mpmath as mp

mp.mp.dps = 100
OUT = Path(__file__).with_suffix('.json')
checks = []

def check(name, actual, expected, tolerance='1e-70', scale=None):
    actual, expected = mp.mpf(actual), mp.mpf(expected)
    den = max(abs(expected), mp.mpf('1e-80')) if scale is None else mp.mpf(scale)
    err = abs(actual - expected) / den
    checks.append(dict(name=name, actual=mp.nstr(actual, 35), expected=mp.nstr(expected, 35),
                       relative_error=mp.nstr(err, 12), tolerance=tolerance, passed=bool(err < mp.mpf(tolerance))))

def geom(t, eta, a, A, n):
    delta = t - 1
    return (t + eta*delta**2 + a*delta**4 + A*t**(-n),
            1 + 2*eta*delta + 4*a*delta**3 - n*A*t**(-n-1),
            2*eta + 12*a*delta**2 + n*(n+1)*A*t**(-n-2),
            24*a*delta - n*(n+1)*(n+2)*A*t**(-n-3),
            24*a + n*(n+1)*(n+2)*(n+3)*A*t**(-n-4))

def retune(A, n, a=mp.mpf(1)):
    t = mp.findroot(lambda t: (t-1)*t**(n+3)-n*(n+1)*(n+2)*A/(24*a), (mp.mpf('.999'),mp.mpf('1.001')), tol=mp.mpf('1e-95'))
    eta = -6*a*(t-1)**2-n*(n+1)*A*t**(-n-2)/2
    return t, eta

def hjet(t, hA):
    return 1+hA*t**-3, -3*hA*t**-4, 12*hA*t**-5

def full_potential(t, y, eta, a, A, n, hA, Wc, w, F=mp.mpf('.07'), P=mp.mpf(1)):
    """Exact two-chiral K=-3P log Ω matrix, at v=0 and fixed real χ."""
    g, p, q, _, _ = geom(t,eta,a,A,n)
    h, hp, hpp = hjet(t,hA)
    Om = g-h*F**2*y**2/(3*P)
    Ot = p-hp*F**2*y**2/(3*P)
    Ott = q-hpp*F**2*y**2/(3*P)
    Oy = -2*h*F**2*y/(3*P)
    Oyy = -2*h*F**2/(3*P)
    Oty = -2*hp*F**2*y/(3*P)
    Kt, Ky = -3*P*Ot/Om, -3*P*Oy/Om
    Ktt = -3*P*(Ott/Om-Ot**2/Om**2)
    Kyy = -3*P*(Oyy/Om-Oy**2/Om**2)
    Kty = -3*P*(Oty/Om-Ot*Oy/Om**2)
    r = mp.sqrt(2)*F
    metric = mp.matrix([[Ktt,1j*Kty/r],[-1j*Kty/r,Kyy/r**2]])
    W = Wc+w*mp.exp(-1j*y)
    Wz = -mp.sqrt(2)*w/F*mp.exp(-1j*y)
    deriv = mp.matrix([Kt*W/P, Wz-1j*Ky*W/(r*P)])
    dd = (deriv.conjugate().T * metric**-1 * deriv)[0]
    return mp.re(Om**-3*(dd-3*abs(W)**2/P))

rows=[]
for n in (1,2,3):
    for a in (mp.mpf('.5'),mp.mpf(1),mp.mpf(2)):
        for A in (mp.mpf('-1e-4'),mp.mpf('-1e-8'),mp.mpf('1e-8'),mp.mpf('1e-4')):
            t,eta=retune(A,n,a)
            g,p,q,qt,qtt=geom(t,eta,a,A,n)
            tag=f'n{n}_a{a}_A{A}'
            check(tag+'_q',q,0,scale=1)
            check(tag+'_qt',qt,0,scale=1)
            check(tag+'_metric_positive',int(g>0 and p*p>0),1)
            check(tag+'_radial_positive',int(qtt>0),1)
            # Full K matrix versus exact axial reduction, with a nonzero clock.
            Wc=mp.mpc('.4','.3');w=mp.mpf('-1e-5');F=mp.mpf('.07')
            U=2*w*w/(F*F)
            exact=full_potential(t,0,eta,a,A,n,0,Wc,w,F)
            check(tag+'_axial_clock',exact,U/g**2)
            for y in (mp.mpf('0'),mp.mpf('.13')):
                v=full_potential(t,y,eta,a,A,n,0,Wc,0,F)
                check(tag+f'_uniform_hidden_y{y}',v,0,scale=1)
            rows.append(dict(n=n,a=str(a),A=str(A),t=mp.nstr(t,35),eta=mp.nstr(eta,35),g=mp.nstr(g,35),p=mp.nstr(p,35)))

# Extended clock metric: source-sign geometry -A/t², h=1+κA/t³.
for A in (mp.mpf('1e-3'),mp.mpf('1e-6'),mp.mpf('1e-12')):
    t,eta=retune(-A,2)
    g,p,q,qt,qtt=geom(t,eta,1,-A,2)
    for kappa in (-1,0,1):
        hA=kappa*A
        h,hp,hpp=hjet(t,hA)
        for theta in (mp.mpf(0),mp.pi/2,mp.pi):
            Wc=mp.exp(1j*theta)*g**mp.mpf('1.5')
            F=mp.mpf('.07')
            f=lambda yy: full_potential(t,yy,eta,1,-A,2,hA,Wc,0,F)
            my=mp.diff(f,0,2)/(F**2*h/g)
            formula=2*g**2/p**2*(2*(hp/h)**2-hpp/h)
            tag=f'extended_A{A}_k{kappa}_th{mp.nstr(theta,5)}'
            check(tag+'_mass',my,formula,scale=max(abs(formula),mp.mpf('1e-20')))
            w=mp.mpf('-1e-5')
            fclock=lambda yy:full_potential(t,yy,eta,1,-A,2,hA,Wc,w,F)
            force=mp.diff(fclock,0)
            formula=-4*w*mp.im(Wc)/(g**2)*hp/(p*h)
            check(tag+'_linear_y',force,formula,scale=max(abs(formula),mp.mpf('1e-20')))
            # Physical mass normalization at fixed target m32=1, M_i=100, F_phys=.07, U_phys=1e-12.
            Fphys=mp.mpf('.07'); Uphys=mp.mpf('1e-12')
            Fhol=Fphys*mp.sqrt(g/h);Uhol=g*g*h*Uphys
            wnorm=-Fhol*mp.sqrt(Uhol)/mp.sqrt(2)
            check(tag+'_physical_F',Fhol**2*h/g,Fphys**2)
            check(tag+'_physical_U',Uhol/(g*g*h),Uphys)
            check(tag+'_physical_w',wnorm/g**mp.mpf('1.5'),-Fphys*mp.sqrt(Uphys)/mp.sqrt(2))

# Finite-rho exact vacuum equations in the h=1 control.
for A in (mp.mpf('0'),mp.mpf('-1e-5')):
    t0,eta0=retune(A,2)
    for rho in (mp.mpf('1e-8'),mp.mpf('1e-16')):
        def eq(t,eta):
            g,p,q,qt,_=geom(t,eta,1,A,2);D=p*p-g*q
            return 3*q/(g*g*D)-rho, qt-(2*p*q/g-q*q/p)
        t,eta=mp.findroot(eq,(t0,eta0),tol=mp.mpf('1e-90'))
        g,p,q,qt,_=geom(t,eta,1,A,2);D=p*p-g*q
        V=3*q/(g*g*D)
        check(f'rho_A{A}_{rho}_V',V,rho)
        check(f'rho_A{A}_{rho}_stationary',qt-(2*p*q/g-q*q/p),0,scale=1)
        F=mp.mpf('.07')
        f=lambda yy:full_potential(t,yy,eta,1,A,2,0,1,0,F)
        my=mp.diff(f,0,2)/(F**2/g)
        mg2=1/g**3
        check(f'rho_A{A}_{rho}_y_mass',my,4*rho/3+2*rho**2/(9*mg2))
        Fcomp=-q/(mp.sqrt(g)*D)
        check(f'rho_A{A}_{rho}_Fcomp',Fcomp,-rho/(3*mp.sqrt(mg2)))

# Clock-force radial displacement at rho=0 (uniform metric) tested at small U.
t,eta=retune(mp.mpf('-1e-5'),2)
g,p,q,qt,qtt=geom(t,eta,1,mp.mpf('-1e-5'),2)
U=mp.mpf('1e-20');Wc=mp.mpf(1);F=mp.mpf('.07');w=-F*mp.sqrt(U)/mp.sqrt(2)
f=lambda xx:full_potential(xx,0,eta,1,mp.mpf('-1e-5'),2,0,Wc,w,F)
shift_formula=2*p**3*U/(3*g*abs(Wc+w)**2*qtt)
tclock=mp.findroot(lambda xx:mp.diff(f,xx),(t,t+2*shift_formula),tol=mp.mpf('1e-90'))
check('clock_radial_displacement_leading',tclock-t,shift_formula,tolerance='1e-18')

P=mp.mpf('2.435e27')**2;rho=mp.mpf('2.5181378331717977e-11')
C=mp.zeta(3)/(48*mp.pi**4)
threshold=mp.sqrt(rho/(18*C))
summary=dict(scope='Independent exact two-chiral K-matrix checks and local geometry retuning; not full SUGRA loops or a 5D embedding',precision_decimal_digits=mp.mp.dps,
             check_count=len(checks),passed=sum(x['passed'] for x in checks),all_pass=all(x['passed'] for x in checks),
             source_sign_C=mp.nstr(C,35),rho_dominated_mg_times_mKK_max_eV2=mp.nstr(threshold,35),
             checks=checks,local_vacua=rows,
             code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             derivation_note='All matrices constructed from Ω derivatives directly; no import of primary experiment or older check code.',
             aborted_diagnostic='An initial unrestricted SymPy inversion of functional g(t),h(t) was manually stopped after about 3 minutes without results; finite jet derivation and full numerical K-matrix replaced it. No tolerance changed.')
OUT.write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:summary[k] for k in ('check_count','passed','all_pass','rho_dominated_mg_times_mKK_max_eV2')},indent=2))
if not summary['all_pass']:
    print(json.dumps([x for x in checks if not x['passed']],indent=2))
    raise SystemExit(1)
