#!/usr/bin/env python3
"""Verify archive consistency without treating it as physical validation."""
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
BASELINE='702102bd86859a106bd5ca31b04e95ac02491d47'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def write(p,d):p.write_text(json.dumps(d,indent=2,ensure_ascii=False)+'\n')

def files():
    own=[p for p in ROOT.rglob('*') if p.is_file() and p not in (MANIFEST,AUDIT)
         and '__pycache__' not in p.parts]
    related=list((REPO/'reports').glob('threshold_*_20260926*'))
    related+=list((REPO/'reports').glob('chain_locality_breaking_20260926*'))
    related += [REPO/p for p in ('README.md','WORK_PLAN.md','RESEARCH_DIRECTION.md',
        'experiments/mediator_chain_v01/results/inputs.json',
        'experiments/persistent_restoration_v01/results/inputs.json',
        'experiments/sugra_shift_v01/results/trajectories.csv',
        'reports/mediator_auxiliary_response_20260926_checks.json')]
    return sorted(set(own+related))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write-manifest',action='store_true')
    args=parser.parse_args()
    if args.write_manifest:
        write(MANIFEST,dict(created_utc=datetime.now(timezone.utc).isoformat(),stage='threshold_rg_v01',
            baseline_git_commit=BASELINE,
            material_type='Selected one-loop leading-source Gaussian matching, low-energy power counting and prescribed deformations.',
            input='Archived chain spectra, fixed finite boundary, old auxiliary-source clock coordinates.',
            output='Threshold/RG scale identities, selected force budgets, light-loop parameters and bypass sensitivity.',
            authorship='AI-assisted project-internal independent derivations and implementations; no external peer review.',
            source_read_records=['reports/threshold_rg_derivation_20260926.md',
                                 'reports/threshold_power_counting_20260926_checks.json',
                                 'reports/chain_locality_breaking_20260926_sources.json'],
            exclusions=['all-order interacting resummation','all UV finite matching or full supergravity determinant',
                        'a computed nonzero bypass or universal quantum floor','complete Higgs stabilization or protected endpoint isolation',
                        'complete new vacuum and cosmic trajectory','observational inference',
                        'Riemann derivation of inverse-log exponent or real interactions'],
            failures_and_limits=[
                'Initial RG derivative implementation93/105 retained; dimensionless differentiation fixes numerical conditioning without tolerance relaxation.',
                'Independent comparison initial18/20 retained; derivative working precision raised without changing1e-80 tolerance.',
                'Historical old-scale approximation failure remains in its immutable archive; current leadingS2 values are explicitly expansions.',
                'A loop/new-tree ratio aboveone remains; actual light-sector repeated-loop parameters are separately computed.',
                'Ordinary spectrum positivity does not protect the tiny endpoint cross coefficient; a minimal Higgs mass origin leaves one flat chiral mode.'],
            files=[dict(path=str(p.relative_to(REPO)),bytes=p.stat().st_size,sha256=sha(p)) for p in files()],
            not_archived=['manifest/audit self-reference','ignored source fulltexts in build/','Python caches']))
    m=read(MANIFEST);checks=[]
    def check(name,ok,detail=''):checks.append(dict(name=name,passed=bool(ok),detail=detail))
    def hashes(values,label):
        for path,h in values.items():
            p=REPO/path;check(label+':'+path,p.is_file() and sha(p)==h)
    for item in m['files']:
        p=REPO/item['path']
        check('hash:'+item['path'],p.is_file() and sha(p)==item['sha256'])
        check('bytes:'+item['path'],p.is_file() and p.stat().st_size==item['bytes'])
        if p.suffix=='.md':
            for link in re.findall(r'!?\[[^\]]*\]\(([^)\n]+)\)',p.read_text()):
                if '://' in link or link.startswith(('#','mailto:')):continue
                dest=(p.parent/link.split('#',1)[0].strip('<>')).resolve()
                check('link:'+item['path']+':'+link,dest.exists() or dest==AUDIT.resolve())
    check('inventory_complete',set(str(p.relative_to(REPO)) for p in files())==set(x['path'] for x in m['files']))
    check('baseline_record',m['baseline_git_commit']==BASELINE)
    for name,count in [('summary',376),('independent_comparison',20)]:
        r=read(ROOT/'results'/f'{name}.json')
        check(name+'_count',r['checks_passed']==r['check_count']==count)
        check(name+'_individual',all(v['passed'] for v in r['checks']))
        hashes(r['source_hashes'],name+'_input_hash')
    ind=read(ROOT/'results/independent_comparison.json')
    check('no_main_source_import_in_independent',ind['main_source_read_or_imported'] is False)
    initial=list((ROOT/'results').glob('independent_failure_*.json'))
    check('independent_failure_preserved',len(initial)==1 and read(initial[0])['checks_passed']==18
          and read(initial[0])['check_count']==20)
    if initial:
        old=read(initial[0])['source_hashes']['experiments/threshold_rg_v01/code/independent_check.py']
        check('independent_initial_source_preserved',sha(ROOT/'code/independent_initial_check.py')==old)
    for stem,count,key in [('threshold_rg_derivation',105,'total'),
                           ('threshold_power_counting',67,'total'),
                           ('chain_locality_breaking',249,'check_count')]:
        p=REPO/'reports'/f'{stem}_20260926_checks.json';r=read(p)
        check(stem+'_counts',r['passed']==r[key]==count)
        check(stem+'_individual',all(v['passed'] for v in r['checks']))
        check(stem+'_script_hash',sha(p.with_suffix('.py'))==r['script_sha256'])
        hashes(r.get('input_hashes',r.get('read_only_inputs',{})),stem+'_input_hash')
        if 'input_sha256' in r:check(stem+'_input',sha(REPO/r['input_path'])==r['input_sha256'])
    r=read(REPO/'reports/threshold_rg_derivation_20260926_initial_checks.json')
    check('rg_initial_failure_preserved',r['passed']==93 and r['total']==105)
    for name,count in [('scale_audit.csv',350),('benchmarks.csv',10),('mass_scan.csv',65),
                       ('clock_profiles.csv',2050),('scale_curve.csv',141)]:
        p=ROOT/'results'/name
        with p.open() as f:n=sum(1 for _ in csv.DictReader(f))
        check('rows:'+name,n==count,str(n));check('LF:'+name,b'\r' not in p.read_bytes())
    for ext in ('.png','.pdf'):
        p=ROOT/'figures'/('threshold_running_and_force'+ext)
        check('figure:'+ext,p.is_file() and p.stat().st_size>1000)
    changed=subprocess.check_output(['git','diff','--name-only',BASELINE,'--','experiments','manuscript_cn'],
                                    cwd=REPO,text=True).splitlines()
    other=[p for p in changed if not p.startswith('experiments/threshold_rg_v01/')]
    check('historical_experiments_and_manuscript_unchanged',not other,str(other))
    result=dict(created_utc=datetime.now(timezone.utc).isoformat(),manifest_sha256=sha(MANIFEST),
                audit_script_sha256=sha(Path(__file__)),passed=sum(v['passed'] for v in checks),total=len(checks),
                checks=checks,scope='Archive consistency and retained records only; not physical verification.')
    write(AUDIT,result)
    print('Archive consistency',result['passed'],'/',result['total'])
    for c in checks:
        if not c['passed']:print(c)
    raise SystemExit(0 if result['passed']==result['total'] else 1)

if __name__=='__main__':main()
