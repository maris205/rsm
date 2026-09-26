#!/usr/bin/env python3
"""Freeze and verify this bounded theoretical stage's materials."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import csv
import hashlib
import json
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
MANIFEST = ROOT / 'materials_manifest.json'
AUDIT = ROOT / 'results/archive_audit.json'
BASELINE = '64d06629a7179e6fed8c5ce7a2af0ea8e17d77ce'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read(p):
    return json.loads(p.read_text())


def write(p, data):
    p.write_text(json.dumps(data, indent=2, ensure_ascii=False) + '\n')


def files():
    own = [p for p in ROOT.rglob('*') if p.is_file() and p not in (MANIFEST, AUDIT)
           and '__pycache__' not in p.parts]
    related = list((REPO / 'reports').glob('clock_*_20260926*'))
    related += [REPO / n for n in ['README.md', 'WORK_PLAN.md', 'RESEARCH_DIRECTION.md',
                'experiments/sequestering_test_v01/results/inputs.json',
                'experiments/sugra_shift_v01/results/inputs.json',
                'experiments/sugra_shift_v01/results/trajectories.csv']]
    return sorted(set(own + related))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-manifest', action='store_true')
    args = parser.parse_args()
    if args.write_manifest:
        write(MANIFEST, dict(created_utc=datetime.now(timezone.utc).isoformat(),
              stage='clock_stabilization_v01', baseline_git_commit=BASELINE,
              material_type='Conditional local EFT derivation and deterministic computation; no new observations.',
              input='Archived K/W and clock coordinates, primary source conventions, declared shared-current mediator.',
              output='Finite-interval curvature, selected contact feedback, hierarchy and local valley diagnostics.',
              scope_exclusions=['complete heavy mediator spectrum and loop matching', 'full 5D embedding',
                                'global or rolling cosmic stability', 'stable U=0 vacuum',
                                'observational alpha fit', 'Riemann derivation of interactions or inverse-log exponent'],
              authorship='AI-assisted project-internal separate code paths; no external peer review.',
              source_read_records=['reports/clock_stabilizer_sources_20260926_sources.json'],
              failures_and_changes=['initial 100/160-digit 346/418 run retained; precision raised without relaxing tolerances',
                                    'initial 20/21 comparison retained; analytic valley cancellation corrected',
                                    'signed beta convention corrected; historical absolute budgets unchanged',
                                    'original 1 eV failure and future U=0 instability retained'],
              files=[dict(path=str(p.relative_to(REPO)), bytes=p.stat().st_size, sha256=sha(p)) for p in files()],
              excluded=['this manifest and its audit to avoid self-reference', 'ignored source fulltexts in build/', 'Python caches']))
    manifest = read(MANIFEST)
    checks = []

    def check(name, passed, detail=''):
        checks.append(dict(name=name, passed=bool(passed), detail=detail))

    for item in manifest['files']:
        p = REPO / item['path']
        check('hash:' + item['path'], p.exists() and sha(p) == item['sha256'])
        check('bytes:' + item['path'], p.exists() and p.stat().st_size == item['bytes'])
        if p.suffix == '.md':
            for target in re.findall(r'!?\[[^\]]*\]\(([^)\n]+)\)', p.read_text()):
                if '://' in target or target.startswith(('#', 'mailto:')):
                    continue
                dest = target.split('#', 1)[0].strip('<>')
                check('link:' + item['path'] + ':' + target, (p.parent / dest).exists())
    check('inventory_complete', set(str(p.relative_to(REPO)) for p in files())
          == set(i['path'] for i in manifest['files']))
    check('baseline_record', manifest['baseline_git_commit'] == BASELINE)
    for name, expected, countkey, hashkey in [('summary', 117, 'check_count', 'source_hashes'),
                                            ('independent_summary', 490, 'checks_total', 'hashes'),
                                            ('independent_comparison', 21, 'checks_total', 'hashes')]:
        r = read(ROOT / 'results' / (name + '.json'))
        check(name + '_counts', r['checks_passed'] == r[countkey] == expected)
        check(name + '_individual', all(c['passed'] for c in r['checks']))
        for path, h in r[hashkey].items():
            check(name + '_hash:' + path, sha(REPO / path) == h)
    ind = read(ROOT / 'results/independent_summary.json')
    cmp = read(ROOT / 'results/independent_comparison.json')
    check('independent_precisions', ind['precisions'] == [180, 240])
    check('independent_points_per_precision', ind['high_precision_points_per_precision'] == 36)
    check('no_output_classification_or_missing_value_differences', cmp['differences'] == [])
    for stem, expected, countkey, hashkey in [
        ('clock_quartic_geometry_20260926_checks', 843, 'check_count', 'code_sha256'),
        ('clock_quartic_geometry_20260926_review', 62, 'check_count', 'code_sha256'),
        ('clock_stabilizer_sources_20260926_checks', 25, 'total', 'script_sha256')]:
        p = REPO / 'reports' / (stem + '.json')
        r = read(p)
        check(stem + '_counts', r['passed'] == r[countkey] == expected)
        check(stem + '_individual', all(c['passed'] for c in r['checks']))
        check(stem + '_script_hash', sha(p.with_suffix('.py')) == r[hashkey])
        for path, h in r.get('source_hashes', {}).items():
            check(stem + '_input_hash:' + path, sha(REPO / path) == h)
    for name, count in [('scan.csv', 4802), ('representatives.csv', 54), ('sigma_sensitivity.csv', 50),
                        ('conditional_frontiers.csv', 147), ('independent_high_precision.csv', 72),
                        ('independent_representatives.csv', 108)]:
        with (ROOT / 'results' / name).open() as stream:
            n = sum(1 for _ in csv.DictReader(stream))
        check('rows:' + name, n == count, str(n))
    for p in (ROOT / 'results').glob('*.csv'):
        check('LF:' + p.name, b'\r' not in p.read_bytes())
    for name, passed, total in [('independent_initial_summary', 346, 418), ('independent_initial_comparison', 20, 21)]:
        r = read(ROOT / 'results' / (name + '.json'))
        check('failure_preserved:' + name, r['checks_passed'] == passed and r['checks_total'] == total)
    for name in ['initial_signed_convention_summary.json', 'pre_valley_cancellation_summary.json']:
        check('previous_version_preserved:' + name, (ROOT / 'results' / name).is_file())
    changed = subprocess.check_output(['git', 'diff', '--name-only', BASELINE, '--', 'experiments', 'manuscript_cn'],
                                      cwd=REPO, text=True).splitlines()
    unexpected = [n for n in changed if not n.startswith('experiments/clock_stabilization_v01/')]
    check('historical_experiments_and_manuscript_unchanged', not unexpected, str(unexpected))
    result = dict(created_utc=datetime.now(timezone.utc).isoformat(), manifest_sha256=sha(MANIFEST),
                  audit_script_sha256=sha(Path(__file__)), passed=sum(c['passed'] for c in checks),
                  total=len(checks), checks=checks,
                  scope='Material consistency, local links and recorded verification only; no scientific certification.')
    write(AUDIT, result)
    print('Archive', result['passed'], '/', result['total'])
    for c in checks:
        if not c['passed']:
            print(c)
    raise SystemExit(0 if result['passed'] == result['total'] else 1)


if __name__ == '__main__':
    main()
