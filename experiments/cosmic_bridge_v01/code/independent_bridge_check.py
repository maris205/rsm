#!/usr/bin/env python3
"""Small independent background/units check, not a PM production matrix.

Reads saved products; never imports their background or reference-spectrum
production functions. Uses radiation-scaled physical velocity and Radau.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.integrate import solve_ivp, quad

ROOT = Path(__file__).resolve().parents[1]
H0S = 67.4 / 3.0856775814913673e19
H0GYR = H0S * 365.25 * 86400 * 1e9
TSTAR = 5.391247e-44
OM, OR, AI, APM = .315, 9.2e-5, 1e-7, .02


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run():
    summary_path = ROOT / 'results/background_summary.json'
    summary = json.loads(summary_path.read_text())
    out = {'scope':'Independent Radau background with v=chi_dot/(H0 E_rad), '
           'z=H0 E_rad t; independent D,D_N growth. No new PM production.',
           'models':{}, 'checks':[], 'input_sha256':{}}

    def record(name, value, threshold):
        out['checks'].append(dict(name=name,value=float(value),threshold=threshold,
                                  passed=bool(np.isfinite(value) and abs(value)<=threshold)))

    ref = summary['models']['eps0']
    ref_initial_E = ref['initial_pm']['E']
    ref_initial_f = ref['regular_f_PM']
    ni, npm = np.log(AI), np.log(APM)
    for label, model in summary['models'].items():
        ep, ol = model['epsilon'], model['omega_lambda']
        path = ROOT / f'results/background_{label}.csv'
        saved = np.genfromtxt(path,delimiter=',',names=True)
        ngrid = saved['N']
        if ep == 0:
            def values(n):
                a = np.exp(n)
                e2 = OR/a**4+OM/a**3+ol
                accel = -(2*OR/a**4+1.5*OM/a**3)/e2
                return np.sqrt(e2), accel
            record(label+'_E_analytic',np.max(abs(values(ngrid)[0]/saved['E']-1)),1e-10)
            for a in (.02,.1,.5,1.):
                tau = quad(lambda b:b/np.sqrt(OR+OM*b+ol*b**4),0,a,
                           epsabs=1e-13,epsrel=3e-12)[0]
                from scipy.interpolate import PchipInterpolator
                published = PchipInterpolator(ngrid,saved['t_Gyr'])(np.log(a))
                record(label+f'_age_analytic_a{a:g}',published/(tau/H0GYR)-1,1e-7)
        else:
            cr = np.sqrt(OR/(1-ep))
            eri = cr/AI**2
            matter_ratio = OM*AI/OR
            corr = (1-ep)/(1-2*ep/3)*matter_ratio
            wi_scaled = 1+corr
            qi = 2-corr/2
            e2i = ((1-ep)*(1+matter_ratio)+ol/eri**2+ep*wi_scaled/3)/(1-ep*qi**2/6)
            vi = qi*np.sqrt(e2i)
            zi = .5*(1-corr/3)
            chii = -.5*np.log(4*H0S**2*TSTAR**2*eri**2*wi_scaled)

            def components(n,state):
                u,v,z = state
                er = cr*np.exp(-2*n)
                wscaled = wi_scaled*np.exp(-2*u+4*(n-ni))
                eratio2 = (1-ep)*(1+OM*np.exp(n)/OR)+ol/er**2+ep*(v*v/6+wscaled/3)
                return er, np.sqrt(eratio2), wscaled

            def rhs(n,state):
                er,e,w = components(n,state)
                u,v,z = state
                return (v/e,-v+2*w/e,-2*z+1/e)

            solution = solve_ivp(rhs,(ni,0),(0,vi,zi),method='Radau',
                                 rtol=2e-11,atol=1e-13,max_step=.03,dense_output=True)
            if not solution.success:
                raise RuntimeError(solution.message)

            def values(n):
                state = solution.sol(n)
                er,e,w = components(n,state)
                efull = er*e
                q = state[1]/e
                a = np.exp(n)
                accel = -(2*OR/a**4+1.5*OM/a**3)/efull**2-ep*q*q/2
                return efull,accel

            state = solution.sol(ngrid)
            er,e,w = components(ngrid,state)
            recovered = {'E':er*e,'q':state[1]/e,
                         'chi':chii+state[0],'t_Gyr':state[2]/er/H0GYR}
            for field, val in recovered.items():
                record(label+'_Radau_'+field,np.max(abs(val/saved[field]-1)),1e-7)
            record(label+'_potential_normalization',
                   chii-model['initial_chi'],1e-11)
            out['models'][label] = {'Radau_evaluations':solution.nfev}

        fmatched = ref_initial_f*ref_initial_E/values(npm)[0]
        record(label+'_matched_initial_f',fmatched/model['initial_pm']['f_matched']-1,1e-7)

        def grow(n,state):
            e,accel = values(n)
            d,dn = state
            return dn,1.5*OM*np.exp(-3*n)/e**2*d-(2+accel)*dn

        growth = solve_ivp(grow,(npm,0),(1,fmatched),method='Radau',
                           rtol=2e-11,atol=1e-13,max_step=.025,dense_output=True)
        if not growth.success:
            raise RuntimeError(growth.message)
        gd,gdn = growth.sol(ngrid)
        record(label+'_matched_D_independent',np.max(abs(gd/saved['D_matched']-1)),1e-7)
        record(label+'_matched_f_independent',np.max(abs((gdn/gd)/saved['f_matched']-1)),1e-7)
        out['input_sha256'][str(path.relative_to(ROOT))] = sha(path)

    meta_path=ROOT/'inputs/reference_spectrum_metadata.json'
    meta=json.loads(meta_path.read_text())
    for filename,item in meta['output_files'].items():
        record('CAMB_saved_hash_'+filename,0 if sha(ROOT/'inputs'/filename)==item['sha256'] else 1,0)
    for filename,key in (('code/make_reference_spectrum.py','generator_sha256'),
                         ('protocol.md','protocol_sha256'),
                         ('REFERENCE_SPECTRUM_INPUT.md','input_note_sha256')):
        record('CAMB_provenance_'+key,0 if sha(ROOT/filename)==meta[key] else 1,0)
    csv=np.genfromtxt(ROOT/'inputs/reference_linear_spectrum.csv',delimiter=',',names=True)
    with np.load(ROOT/'inputs/reference_linear_spectrum.npz') as data:
        for field in csv.dtype.names:
            record('CAMB_npz_csv_'+field,np.max(abs(data[field]/csv[field]-1)),1e-14)
    out['input_sha256'][str(summary_path.relative_to(ROOT))]=sha(summary_path)
    out['input_sha256'][str(meta_path.relative_to(ROOT))]=sha(meta_path)
    out['code_sha256']=sha(__file__)
    out['passed']=sum(c['passed'] for c in out['checks'])
    out['total']=len(out['checks'])
    out['failures']=[c for c in out['checks'] if not c['passed']]
    target=ROOT/'reports/independent_bridge_checks.json'
    target.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'passed':out['passed'],'total':out['total'],'failures':out['failures']},indent=2))


if __name__=='__main__':
    run()
