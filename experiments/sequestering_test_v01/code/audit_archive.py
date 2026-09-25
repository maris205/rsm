#!/usr/bin/env python3
"""Verify this stage's archived materials, not scientific validity of the model."""
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
PARENT = '18bf8a90abf5bd5d3b385af99978018ad4b9c054'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def write(path, record):
    path.write_text(json.dumps(record, ensure_ascii=False, indent=2) + '\n')


def files():
    owned = [p for p in ROOT.rglob('*') if p.is_file()
             and p not in (MANIFEST, AUDIT) and '__pycache__' not in p.parts]
    related = list((REPO / 'reports').glob('protection_*_20260925*'))
    related += [REPO / name for name in ('README.md', 'WORK_PLAN.md', 'RESEARCH_DIRECTION.md',
                'experiments/sugra_shift_v01/results/inputs.json',
                'experiments/sugra_shift_v01/results/trajectories.csv')]
    return sorted(set(owned + related))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-manifest', action='store_true')
    args = parser.parse_args()
    if args.write_manifest:
        write(MANIFEST, dict(
            created_utc=datetime.now(timezone.utc).isoformat(),
            stage='sequestering_test_v01', parent_git_commit=PARENT,
            input_type='Primary finite 5D gravity coefficient plus conditional 4D K/W; stored model coordinates, no new observation.',
            output_type='Correlated selected matter contacts and clock curvature; local necessary conditions and independent calculations.',
            scope_exclusions=['full 5D embedding of clock and quartic stabilizer',
                              'complete two-loop matching and finite local potential',
                              'rolling covariant stability or new cosmic trajectory',
                              'Standard Model or observed alpha prediction',
                              'Riemann derivation of interactions or inverse-log exponent'],
            assumptions=['kappa clock-metric extension is explicit input',
                         'nu0=0 and nu0=-c are different finite matching boundaries',
                         'positive local curvature is a necessary condition, not full stability'],
            authorship='AI-assisted project-internal derivation and distinct code paths; no external peer review.',
            failure_provenance=['initial root-residual normalization precision failure retained',
                                'two comparator mass-scaling failures retained; tolerances unchanged',
                                'initial leading-curvature summary retained before finite-clock refinement',
                                'negative-curvature physical branches retained'],
            files=[dict(path=str(p.relative_to(REPO)), bytes=p.stat().st_size, sha256=sha(p)) for p in files()],
            excluded=['this manifest and its audit to avoid self-reference',
                      'ignored primary-source fulltexts in build/', 'Python caches'],
            source_read_records=['reports/protection_5d_20260925_sources.json',
                                 'reports/protection_symmetry_20260925_sources.json']))

    manifest = read(MANIFEST)
    checks = []

    def check(name, passed, detail=''):
        checks.append(dict(name=name, passed=bool(passed), detail=detail))

    for item in manifest['files']:
        p = REPO / item['path']
        check('hash:' + item['path'], p.exists() and sha(p) == item['sha256'])
        check('bytes:' + item['path'], p.exists() and p.stat().st_size == item['bytes'])
    check('inventory_complete', set(str(p.relative_to(REPO)) for p in files())
          == set(i['path'] for i in manifest['files']))
    check('parent_commit_recorded', manifest['parent_git_commit'] == PARENT)
    for item in manifest['files']:
        p = REPO / item['path']
        if p.suffix == '.md':
            for target in re.findall(r'!?\[[^\]]*\]\(([^)\n]+)\)', p.read_text()):
                if '://' in target or target.startswith(('#', 'mailto:')):
                    continue
                dest = target.split('#', 1)[0].strip('<>')
                check('link:' + item['path'] + ':' + target, (p.parent / dest).exists())

    for name, expected in [('summary', 177), ('independent_summary', 153),
                           ('independent_comparison', 22)]:
        record = read(ROOT / 'results' / (name + '.json'))
        check(name + '_checks', record['checks_passed'] == record['check_count'] == expected)
        check(name + '_all_individual_checks', all(c['passed'] for c in record['checks']))
        for path, expected_hash in record['source_hashes'].items():
            check(name + '_input_hash:' + path, sha(REPO / path) == expected_hash)
    independent = read(ROOT / 'results/independent_summary.json')
    compare = read(ROOT / 'results/independent_comparison.json')
    check('independent_precision_points', independent['direct_points_per_precision'] == 72)
    check('independent_precision_digits', independent['precision_digits'] == [140, 220])
    check('classification_agreement', compare['boolean_mismatches'] == [])

    for stem, total, count_key, hash_key in [
        ('protection_symmetry_20260925_checks', 44, 'total', 'script_sha256'),
        ('protection_5d_20260925_checks', 30, 'total', 'script_sha256'),
        ('protection_power_correction_20260925_checks', 404, 'check_count', 'code_sha256'),
    ]:
        path = REPO / 'reports' / (stem + '.json')
        record = read(path)
        check(stem, record['passed'] == record[count_key] == total)
        check(stem + '_code_hash', sha(path.with_suffix('.py')) == record[hash_key])
        check(stem + '_all_individual_checks', all(c['passed'] for c in record['checks']))
    review = read(REPO / 'reports/protection_power_correction_20260925_review.json')
    check('extra_Bjet_review', review['passed'] == review['check_count'] == 24)
    for path, expected_hash in review['source_hashes'].items():
        check('extra_Bjet_review_input_hash:' + path, sha(REPO / path) == expected_hash)

    for filename, expected in [('scan.csv', 19215), ('representatives.csv', 144),
                               ('conditional_limits.csv', 4),
                               ('independent_kahler_140.csv', 72),
                               ('independent_kahler_220.csv', 72)]:
        path = ROOT / 'results' / filename
        with path.open() as stream:
            count = sum(1 for _ in csv.DictReader(stream))
        check('rows:' + filename, count == expected, str(count))
    for p in (ROOT / 'results').glob('*.csv'):
        check('LF:' + p.name, b'\r' not in p.read_bytes())

    for filename, passed, total in [('independent_initial_summary.json', 44, 45),
                                   ('independent_first_comparison.json', 21, 22),
                                   ('independent_second_comparison.json', 21, 22)]:
        record = read(ROOT / 'results' / filename)
        check('failure_preserved:' + filename,
              record['checks_passed'] == passed and record['check_count'] == total)
    check('leading_approximation_record_preserved', (ROOT / 'results/initial_leading_mass_summary.json').exists())

    # Read-only check of the complete historical manuscript and older experiments.
    changed = subprocess.check_output(['git', 'diff', '--name-only', PARENT, '--',
                                      'experiments', 'manuscript_cn'], cwd=REPO, text=True).splitlines()
    unexpected = [name for name in changed if not name.startswith('experiments/sequestering_test_v01/')]
    check('parent_experiments_and_formal_manuscript_unchanged', not unexpected, str(unexpected))

    result = dict(created_utc=datetime.now(timezone.utc).isoformat(),
                  manifest_sha256=sha(MANIFEST), audit_script_sha256=sha(Path(__file__)),
                  passed=sum(c['passed'] for c in checks), total=len(checks), checks=checks,
                  scope='Artifact hashes, local links and current check statuses only; no scientific certification.')
    write(AUDIT, result)
    print('Archive', result['passed'], '/', result['total'])
    for c in checks:
        if not c['passed']:
            print(c)
    raise SystemExit(0 if result['passed'] == result['total'] else 1)


if __name__ == '__main__':
    main()
