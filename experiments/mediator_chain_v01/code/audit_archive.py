#!/usr/bin/env python3
"""Archive consistency, including deliberately retained scientific failures."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import csv
import hashlib
import json
import re
import subprocess

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
MANIFEST=ROOT/'materials_manifest.json'
AUDIT=ROOT/'results/archive_audit.json'
BASELINE='06188a0eaf90b7db1c65faf125da680fb1e3d5d7'


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p): return json.loads(p.read_text())
def write(p,d): p.write_text(json.dumps(d,indent=2,ensure_ascii=False)+'\n')


def files():
    own=[p for p in ROOT.rglob('*') if p.is_file() and p not in (MANIFEST,AUDIT)
         and '__pycache__' not in p.parts]
    related=list((REPO/'reports').glob('mediator_*_20260926*'))
    related += [REPO/p for p in ('README.md','WORK_PLAN.md','RESEARCH_DIRECTION.md',
                 'experiments/persistent_restoration_v01/results/inputs.json',
                 'experiments/sugra_shift_v01/results/trajectories.csv')]
    return sorted(set(own+related))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--write-manifest',action='store_true')
    args=ap.parse_args()
    if args.write_manifest:
        write(MANIFEST,dict(created_utc=datetime.now(timezone.utc).isoformat(),
            stage='mediator_chain_v01',baseline_git_commit=BASELINE,
            material_type='Conditional local EFT chain, static response and selected spectra; no observations.',
            input='Prior physical parameters and clock coordinates; declared gapped nearest-neighbor mediator action.',
            output='Endpoint matching, complete self terms, scoped finite-F response and Gaussian spectrum.',
            authorship='AI-assisted internal separate implementations; no external peer review.',
            exclusions=['proof of protected locality or bulk-modulus endpoint isolation',
                        'complete UV theory or full supergravity Gaussian determinant',
                        'threshold/RG resummation', 'complete common vacuum rematching',
                        'global or rolling cosmic stability', 'Riemann origin of interactions or inverse-log law'],
            source_read_records=['reports/mediator_chain_sources_20260926_sources.json'],
            retained_failures=[
                'The previous 1 TeV/1e-6 eV example fails finite-F response and charged-mass gap conditions.',
                'The one-loop/tree selected ratio exceeds one at the low-scale fixed-mu control; full loop ordering is not established.',
                'Old-scale heavy C^2 slope approximation: 0.247784 percent exceeds declared 0.2 percent accuracy; exact spectrum is used there.',
                'Algebra initial 975/978 at180 digits retained;220 digits pass unchanged tolerances.',
                'Supergravity initial31/43 implementation failure retained; explicit mpmath endpoint index fixes negative-index error.',
                'Auxiliary first pre-calculation matrix-constructor error recorded in execution_history.'],
            files=[dict(path=str(p.relative_to(REPO)),bytes=p.stat().st_size,sha256=sha(p)) for p in files()],
            not_archived=['manifest and audit excluded from their own inventory',
                          'ignored source fulltexts in build/', 'Python caches']))
    m=read(MANIFEST); checks=[]
    def check(name,passed,detail=''):
        checks.append(dict(name=name,passed=bool(passed),detail=detail))
    def hashes(values,label):
        for path,h in values.items():
            p=REPO/path
            check(label+':'+path,p.is_file() and sha(p)==h)
    for item in m['files']:
        p=REPO/item['path']
        check('hash:'+item['path'],p.is_file() and sha(p)==item['sha256'])
        check('bytes:'+item['path'],p.is_file() and p.stat().st_size==item['bytes'])
        if p.suffix=='.md':
            for link in re.findall(r'!?\[[^\]]*\]\(([^)\n]+)\)',p.read_text()):
                if '://' in link or link.startswith(('#','mailto:')): continue
                target=(p.parent/link.split('#',1)[0].strip('<>')).resolve()
                # This exact output is generated at the end of this same audit.
                check('link:'+item['path']+':'+link,target.exists() or target==AUDIT.resolve())
    check('inventory_complete',set(str(p.relative_to(REPO)) for p in files())
          ==set(f['path'] for f in m['files']))
    check('baseline_record',m['baseline_git_commit']==BASELINE)
    mainresult=read(ROOT/'results/summary.json')
    check('main_count',mainresult['checks_passed']==mainresult['check_count']==124)
    check('main_individual',all(c['passed'] for c in mainresult['checks']))
    hashes(mainresult['source_hashes'],'main_input_hash')
    ind=read(ROOT/'results/independent_comparison.json')
    check('independent_count',ind['passed']==ind['total']==21732)
    check('independent_no_discrepancies',ind['failures']==[])
    check('independent_aggregate',sum(a['passed'] for a in ind['aggregate'].values())==21732)
    check('independent_script_hash',sha(ROOT/'code/independent_compare.py')==ind['code_sha256'])
    hashes(ind['source_hashes'],'independent_input_hash')
    for stem,expected in [('mediator_chain_sources',38),('mediator_chain_algebra',978),
                          ('mediator_sugra_response',359)]:
        p=REPO/'reports'/f'{stem}_20260926_checks.json';r=read(p)
        check(stem+'_count',r['passed']==r['total']==expected)
        if isinstance(r.get('checks'),list):
            check(stem+'_individual',all(c['passed'] for c in r['checks']))
        if 'script_sha256' in r:
            check(stem+'_script_hash',sha(p.with_suffix('.py'))==r['script_sha256'])
    auxp=REPO/'reports/mediator_auxiliary_response_20260926_checks.json';aux=read(auxp)
    check('auxiliary_counts_include_retained_failure',aux['passed']==512 and aux['total']==513)
    check('auxiliary_only_declared_approximation_failure',len(aux['failures'])==1 and
          aux['failures'][0]['name']=='C2_asymptotic_at_resolved_splitting_dV1_heavy_C2_dS')
    check('auxiliary_script_hash',sha(auxp.with_suffix('.py'))==aux['script_sha256'])
    check('auxiliary_input_hash',sha(REPO/'experiments/persistent_restoration_v01/results/inputs.json')==aux['input_sha256'])
    for stem,expected_pass,expected_total in [('mediator_chain_algebra',975,978),
                                             ('mediator_sugra_response',31,43)]:
        p=REPO/'reports'/f'{stem}_20260926_initial_checks.json';r=read(p)
        check(stem+'_initial_preserved',r['passed']==expected_pass and r['total']==expected_total)
    low=read(REPO/'reports/mediator_sugra_response_20260926_150_checks.json')
    check('sugra150_counts',low['passed']==low['total']==161)
    check('sugra150_same_final_code',low['script_sha256']==sha(REPO/'reports/mediator_sugra_response_20260926_checks.py'))
    for name,count in [('chain_matching.csv',777),('response_scan.csv',390),
                       ('representatives.csv',36),('response_limits.csv',6)]:
        p=ROOT/'results'/name
        with p.open() as f: n=sum(1 for _ in csv.DictReader(f))
        check('rows:'+name,n==count,str(n))
        check('LF:'+name,b'\r' not in p.read_bytes())
    for ext in ('.png','.pdf'):
        p=ROOT/'figures'/('chain_matching_and_response'+ext)
        check('figure:'+ext,p.is_file() and p.stat().st_size>1000)
    changed=subprocess.check_output(['git','diff','--name-only',BASELINE,'--','experiments','manuscript_cn'],
                                    cwd=REPO,text=True).splitlines()
    unexpected=[p for p in changed if not p.startswith('experiments/mediator_chain_v01/')]
    check('historical_experiments_and_manuscript_unchanged',not unexpected,str(unexpected))
    result=dict(created_utc=datetime.now(timezone.utc).isoformat(),manifest_sha256=sha(MANIFEST),
                audit_script_sha256=sha(Path(__file__)),passed=sum(c['passed'] for c in checks),
                total=len(checks),checks=checks,
                scope='Material consistency only. Includes correctly preserved failed scientific conditions and approximation tests.')
    write(AUDIT,result)
    print('Archive consistency',result['passed'],'/',result['total'])
    for c in checks:
        if not c['passed']: print(c)
    raise SystemExit(0 if result['passed']==result['total'] else 1)


if __name__=='__main__': main()
