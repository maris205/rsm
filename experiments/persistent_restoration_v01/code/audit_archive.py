#!/usr/bin/env python3
"""Check material consistency for this bounded theoretical calculation."""
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
BASELINE = 'e4834cbc8cfde8ccdb00ba362ea68cef6cdcd959'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read(p):
    return json.loads(p.read_text())


def write(p, data):
    p.write_text(json.dumps(data, indent=2, ensure_ascii=False) + '\n')


def files():
    own = [p for p in ROOT.rglob('*') if p.is_file()
           and p not in (MANIFEST, AUDIT) and '__pycache__' not in p.parts]
    related = list((REPO / 'reports').glob('persistent_*_20260926*'))
    related += [REPO / n for n in ['README.md', 'WORK_PLAN.md', 'RESEARCH_DIRECTION.md',
                'experiments/sequestering_test_v01/results/inputs.json',
                'experiments/clock_stabilization_v01/results/inputs.json',
                'experiments/sugra_shift_v01/results/inputs.json',
                'experiments/sugra_shift_v01/results/trajectories.csv']]
    return sorted(set(own + related))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-manifest', action='store_true')
    args = parser.parse_args()
    if args.write_manifest:
        write(MANIFEST, dict(created_utc=datetime.now(timezone.utc).isoformat(),
              stage='persistent_restoration_v01', baseline_git_commit=BASELINE,
              material_type='Conditional local EFT derivation and deterministic computation; no new observations.',
              input='Archived K/W, clock coordinates and matching; added explicit frozen modulus shared current.',
              output='Persistent local mass, correlated soft terms, finite feedback window and selected late-force limits.',
              scope_exclusions=['derivation or protection of tiny tau', 'complete heavy spectrum and loop matching',
                                'full 5D embedding', 'global or rolling cosmic stability',
                                'permanent slow real-clock evolution', 'observational alpha fit',
                                'Riemann derivation of interactions or inverse-log exponent',
                                'physical second-order matching at tuned rH=2'],
              authorship='AI-assisted project-internal separate derivations and implementations; no external peer review.',
              source_read_records=['reports/persistent_current_sources_20260926_sources.json'],
              failures_and_changes=[
                  'Initial 100/160-digit 979/1016 and intermediate 991/1012 runs retained.',
                  'Independent stable CW kernel and positive-curvature root domain corrected; final 160/220 digits without tolerance relaxation.',
                  'Initial main future roots outside valid local-curvature domain retained as diagnostics.',
                  'Frozen canonical radial extension differs from prior hZ(t) extension; mixed soft derivative rederived.',
                  'Negative/no-current controls and large-coupling feedback failures retained.',
                  'Tuned signed tail omits known vacuum soft term and controlled second-order normalization.',
                  'Plotter initially rejected blank invalid-domain budgets; missing values now plotted grey, not replaced by zero.'],
              files=[dict(path=str(p.relative_to(REPO)), bytes=p.stat().st_size,
                          sha256=sha(p)) for p in files()],
              excluded=['this manifest and its audit to avoid self-reference',
                        'ignored downloaded source fulltexts in build/', 'Python caches']))
    manifest = read(MANIFEST)
    checks = []

    def check(name, passed, detail=''):
        checks.append(dict(name=name, passed=bool(passed), detail=detail))

    def hashes(mapping, label):
        for path, expected in mapping.items():
            p = REPO / path
            check(label + ':' + path, p.is_file() and sha(p) == expected)

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
    for name, expected, countkey, hashkey in [('summary', 91, 'check_count', 'source_hashes'),
                                            ('independent_summary', 1012, 'checks_total', 'hashes'),
                                            ('independent_comparison', 27, 'checks_total', 'hashes')]:
        r = read(ROOT / 'results' / (name + '.json'))
        check(name + '_counts', r['checks_passed'] == r[countkey] == expected)
        check(name + '_individual', all(c['passed'] for c in r['checks']))
        hashes(r[hashkey], name + '_hash')
    ind = read(ROOT / 'results/independent_summary.json')
    cmp = read(ROOT / 'results/independent_comparison.json')
    check('independent_precisions', ind['precision_digits'] == [160, 220])
    check('independent_high_precision_points', ind['high_precision_points'] == 144)
    check('no_output_classification_or_missing_value_differences', cmp['differences'] == [])
    for stem, expected, countkey in [
        ('persistent_current_geometry_20260926_checks', 559, 'check_count'),
        ('persistent_charged_loop_20260926_checks', 377, 'check_count'),
        ('persistent_current_sources_20260926_checks', 47, 'total')]:
        p = REPO / 'reports' / (stem + '.json')
        r = read(p)
        check(stem + '_counts', r['passed'] == r[countkey] == expected)
        check(stem + '_individual', all(c['passed'] for c in r['checks']))
        if 'code_sha256' in r:
            check(stem + '_script_hash', sha(p.with_suffix('.py')) == r['code_sha256'])
        if 'input_sha256' in r:
            check(stem + '_input_hash', sha(REPO / r['input_path']) == r['input_sha256'])
    review = read(REPO / 'reports/persistent_current_geometry_20260926_review.json')
    hashes(review['read_sha256'], 'implementation_review_read_hash')
    check('review_keeps_original_checks', review['changes_to_main_or_original_checks'] is False)
    check('review_no_earlier_field_distance_failure',
          review['checks']['additional_field_fraction_roots_before_valid_search_stop'] == 0)
    for name, count in [('scan.csv', 3922), ('representatives.csv', 54),
                        ('future_limits.csv', 8), ('future_curves.csv', 2404),
                        ('independent_high_precision.csv', 144),
                        ('independent_representatives.csv', 108), ('independent_future.csv', 24)]:
        with (ROOT / 'results' / name).open() as stream:
            n = sum(1 for _ in csv.DictReader(stream))
        check('rows:' + name, n == count, str(n))
    for p in (ROOT / 'results').glob('*.csv'):
        check('LF:' + p.name, b'\r' not in p.read_bytes())
    for name, passed, total in [('independent_initial_summary', 979, 1016),
                                ('independent_intermediate_summary', 991, 1012)]:
        r = read(ROOT / 'results' / (name + '.json'))
        check('failure_preserved:' + name,
              r['checks_passed'] == passed and r['checks_total'] == total)
    for name in ['initial_future_domain_summary.json', 'initial_future_limits.csv']:
        check('initial_main_version_preserved:' + name, (ROOT / 'results' / name).is_file())
    for stem in ('persistent_window', 'persistent_future_limit'):
        for ext in ('.png', '.pdf'):
            p = ROOT / 'figures' / (stem + ext)
            check('figure:' + p.name, p.is_file() and p.stat().st_size > 1000)
    changed = subprocess.check_output(['git', 'diff', '--name-only', BASELINE, '--',
                                      'experiments', 'manuscript_cn'], cwd=REPO, text=True).splitlines()
    unexpected = [n for n in changed if not n.startswith('experiments/persistent_restoration_v01/')]
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
