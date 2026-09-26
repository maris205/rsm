#!/usr/bin/env python3
"""Archive consistency only; it is not a test of physical truth."""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import csv
import hashlib
import json
import platform
import re
import subprocess
import numpy,scipy,sympy,mpmath,matplotlib

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
MANIFEST=ROOT/'materials_manifest.json'
AUDIT=ROOT/'results/archive_audit.json'
BASELINE='0fe84705869a46277fa4c538f54e06e111428636'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def write(p,x):p.write_text(json.dumps(x,indent=2,ensure_ascii=False)+'\n')
def inventory():
    paths=[p for p in ROOT.rglob('*') if p.is_file() and p not in (MANIFEST,AUDIT) and '__pycache__' not in p.parts]
    paths+=list((REPO/'reports').glob('higgs_*_20260926*'))
    paths += [REPO/p for p in ('README.md','WORK_PLAN.md','RESEARCH_DIRECTION.md',
        'reports/chain_locality_breaking_20260926.md','experiments/threshold_rg_v01/results/inputs.json',
        'experiments/threshold_rg_v01/protocol.md')]
    return sorted(set(paths))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write-manifest',action='store_true');args=parser.parse_args()
    if args.write_manifest:
        write(MANIFEST,dict(created_utc=datetime.now(timezone.utc).isoformat(),stage='higgs_stabilization_v01',
            baseline_git_commit=BASELINE,material_type='Conditional theoretical derivation and deterministic numerical audit.',
            question='Can either of two specified Higgs stabilizations preserve the actual local endpoint restoration channel?',
            input='Prior chain charges; canonical paired Higgs fields; stated polynomial or extra-gauge extensions and local sources.',
            output='Minimum operator degree, free multiplet masses, current matching, finite-momentum kernel and restricted finite-F stationarity.',
            authorship='AI-assisted internal derivations and independent implementations; no external peer review.',
            source_read_record='reports/higgs_phase_operators_20260926_sources.json',
            versions=dict(python=platform.python_version(),numpy=numpy.__version__,scipy=scipy.__version__,
                          sympy=sympy.__version__,mpmath=mpmath.__version__,matplotlib=matplotlib.__version__),
            limitations=['No universal no-go for all Higgs or UV completions.',
                'The all-heavy free spectrum is at F=D=0; full broken-SUSY Hessian and SUGRA remain open.',
                'Finite-F zero clock response is confined to W=fX, X=Za=0 and no clock F on a smooth rigid branch.',
                'Hidden self-response persists; the prior threshold force budget is not inherited.',
                'No observed constant variation, new cosmic trajectory, or derivation of 1/ln²t from Riemann structure.'],
            retained_failures=['Main NumPy-bool JSON serialization failure, original code and log.',
                'Independent syntax failure plus initial3262/4087 output comparison caused by exact-vs-float log-grid inputs.',
                'Phase-operator initial265/266 SymPy variable-assumption comparison failure.',
                'Extra-gauge initial634/635 expanded-vs-factored structural comparison failure.',
                'Physical gap failures and zero local cross contact are results, not repaired tests.'],
            files=[dict(path=str(p.relative_to(REPO)),bytes=p.stat().st_size,sha256=sha(p)) for p in inventory()],
            exclusions=['manifest/audit self-reference','ignored source caches in build/','Python caches']))
    m=read(MANIFEST);checks=[]
    def check(name,passed,detail=''):checks.append(dict(name=name,passed=bool(passed),detail=detail))
    def hashes(values,prefix):
        for path,h in values.items():
            p=REPO/path;check(prefix+':'+path,p.is_file() and sha(p)==h)
    check('baseline_record',m['baseline_git_commit']==BASELINE)
    for item in m['files']:
        p=REPO/item['path'];check('hash:'+item['path'],p.is_file() and sha(p)==item['sha256'])
        check('size:'+item['path'],p.is_file() and p.stat().st_size==item['bytes'])
        if p.suffix=='.md':
            for link in re.findall(r'!?\[[^\]]*\]\(([^)\n]+)\)',p.read_text()):
                if '://' in link or link.startswith(('#','mailto:')):continue
                dest=(p.parent/link.split('#',1)[0].strip('<>')).resolve()
                check('link:'+item['path']+':'+link,dest.exists() or dest==AUDIT.resolve())
    check('inventory_complete',set(str(p.relative_to(REPO)) for p in inventory())==set(i['path'] for i in m['files']))
    r=read(ROOT/'results/summary.json')
    check('main_counts',r['checks_passed']==r['check_count']==245)
    check('main_checks',all(c['passed'] for c in r['checks']));hashes(r['input_hashes'],'main_input')
    specs=[('higgs_phase_operators',271,'check_count','passed','code_sha256'),
           ('higgs_extra_gauge',594,'total','passed','code_sha256'),
           ('higgs_local_currents',227,'total','passed','script_sha256'),
           ('higgs_finite_auxiliary',678,'total_checks','passed_checks','script_sha256')]
    for stem,total,countkey,passedkey,hashkey in specs:
        p=REPO/'reports'/f'{stem}_20260926_checks.json';v=read(p)
        check(stem+'_counts',v[countkey]==v[passedkey]==total)
        check(stem+'_checks',all(c['passed'] for c in v['checks']))
        check(stem+'_code_hash',sha(p.with_suffix('.py'))==v[hashkey])
        hashes(v.get('input_hashes',v.get('input_sha256',{})),stem+'_inputs')
    ind=read(ROOT/'results/independent_comparison.json')
    check('independent_counts',ind['checks_total']==ind['checks_passed']==4087 and not ind['failures'])
    check('independent_coverage',ind['rows_total']==560 and ind['csv_cells_total']==4009 and ind['comparison_group_count']==70)
    check('independent_no_main_import',ind['main_implementation_not_read'] is True)
    check('independent_code_hash',sha(ROOT/'code/independent_check.py')==ind['script_sha256'])
    hashes(ind['input_hashes'],'independent_inputs')
    for rel in ind['failure_history_files']:
        v=read(ROOT/rel);check('independent_initial_retained:'+rel,v['checks_passed']==3262 and v['checks_total']==4087)
        check('independent_initial_code_retained',sha(ROOT/'code/independent_initial_grid_check.py.txt')==v['script_sha256'])
    check('independent_grid_tolerance_unchanged',ind['grid_comparison_correction']['response_tolerance_relaxed'] is False)
    initial=read(REPO/'reports/higgs_phase_operators_20260926_initial_checks.json')
    check('phase_initial_failure',initial['passed']==265 and initial['check_count']==266)
    check('phase_initial_source',sha(REPO/'reports/higgs_phase_operators_20260926_initial_checks.py')==initial['code_sha256'])
    initial=read(REPO/'reports/higgs_extra_gauge_20260926_initial_checks.json')
    check('gauge_initial_failure',initial['passed']==634 and initial['total']==635)
    check('main_serialization_failure_preserved','not JSON serializable' in (ROOT/'logs/main_initial_failure.log').read_text())
    check('main_original_source_preserved','gap_criterion_lambda_over_g_3=1e11/' in (ROOT/'code/run_stabilization_initial.py').read_text())
    for name,expected in [('phase_lifting.csv',108),('gauge_scan.csv',59),('momentum_profile.csv',393)]:
        p=ROOT/'results'/name
        with p.open() as f:count=sum(1 for _ in csv.DictReader(f))
        check('row_count:'+name,count==expected,str(count));check('LF:'+name,b'\r' not in p.read_bytes())
    for ext,magic in [('png',b'\x89PNG'),('pdf',b'%PDF')]:
        p=ROOT/'figures'/('stabilization_tradeoff.'+ext)
        check('figure:'+ext,p.is_file() and p.stat().st_size>1000 and p.read_bytes().startswith(magic))
    changed=subprocess.check_output(['git','diff','--name-only',BASELINE,'--'],cwd=REPO,text=True).splitlines()
    def permitted(p):return p in ('README.md','WORK_PLAN.md','RESEARCH_DIRECTION.md') or p.startswith('experiments/higgs_stabilization_v01/') or (p.startswith('reports/higgs_') and '_20260926' in p)
    disallowed=[p for p in changed if not permitted(p)]
    check('historical_files_and_manuscript_unchanged',not disallowed,str(disallowed))
    result=dict(created_utc=datetime.now(timezone.utc).isoformat(),manifest_sha256=sha(MANIFEST),
        audit_script_sha256=sha(Path(__file__)),passed=sum(c['passed'] for c in checks),total=len(checks),
        checks=checks,scope='Archive identity, completeness and retained checks; not physical validation.')
    write(AUDIT,result);print('Archive consistency',result['passed'],'/',result['total'])
    for c in checks:
        if not c['passed']:print(c)
    raise SystemExit(0 if result['passed']==result['total'] else 1)

if __name__=='__main__':main()
