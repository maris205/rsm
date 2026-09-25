#!/usr/bin/env python3
"""Freeze and verify material provenance for this bounded research stage."""
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


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, data):
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n')


def files():
    owned = [p for p in ROOT.rglob('*') if p.is_file()
             and p not in (MANIFEST, AUDIT)
             and '__pycache__' not in p.parts]
    related = list((REPO / 'reports').glob('matching_*_20260925*'))
    related += [REPO / p for p in ('README.md', 'WORK_PLAN.md', 'RESEARCH_DIRECTION.md',
                                    'experiments/sugra_shift_v01/results/inputs.json',
                                    'experiments/sugra_shift_v01/results/trajectories.csv',
                                    'reports/noscale_quantum_sources_20260925_sources.json')]
    return sorted(set(owned + related))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-manifest', action='store_true')
    args = parser.parse_args()
    if args.write_manifest:
        write(MANIFEST, dict(created_utc=datetime.now(timezone.utc).isoformat(),
                            stage='eft_matching_v01',
                            parent_git_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip(),
                            input_type='Analytic conditional EFT plus stored model coordinates; no new observations.',
                            output_type='Specified local operator sensitivities, selected RG closure, independent symbolic/numerical audits.',
                            scope_exclusions=['complete O(Q²) SUGRA matching', 'full two-loop potential',
                                              'UV Wilson-coefficient predictions', 'new cosmological or observed alpha fit'],
                            authorship='AI-assisted project-internal derivation and separate code paths; no external peer review.',
                            files=[dict(path=str(p.relative_to(REPO)), bytes=p.stat().st_size, sha256=sha(p)) for p in files()],
                            excluded=['this manifest and its audit (avoid self-reference)', 'ignored source fulltexts in build/',
                                      'Python caches'],
                            source_read_records=['reports/matching_gauge_rg_20260925_sources.json',
                                                 'reports/matching_mixed_sector_20260925_sources.json',
                                                 'reports/noscale_quantum_sources_20260925_sources.json']))

    manifest = json.loads(MANIFEST.read_text())
    checks = []

    def check(name, passed, detail=''):
        checks.append(dict(name=name, passed=bool(passed), detail=detail))

    for item in manifest['files']:
        p = REPO / item['path']
        check('hash:' + item['path'], p.exists() and sha(p) == item['sha256'])
    check('inventory_complete', set(str(p.relative_to(REPO)) for p in files()) == set(i['path'] for i in manifest['files']))
    for item in manifest['files']:
        p = REPO / item['path']
        if p.suffix == '.md':
            for target in re.findall(r'!?\[[^\]]*\]\(([^)\n]+)\)', p.read_text()):
                if '://' in target or target.startswith(('#', 'mailto:')):
                    continue
                dest = target.split('#', 1)[0].strip('<>')
                check('link:' + item['path'] + ':' + target, (p.parent / dest).exists())

    summary = json.loads((ROOT / 'results/summary.json').read_text())
    check('main_checks', summary['checks_passed'] == summary['check_count'] == 300)
    for path, expected in summary['source_hashes'].items():
        check('main_input_hash:' + path, sha(REPO / path) == expected)
    independent = json.loads((ROOT / 'results/independent_summary.json').read_text())
    comparison = json.loads((ROOT / 'results/independent_comparison.json').read_text())
    check('independent_checks', independent['passed'] == independent['total'] == 22)
    check('independent_comparisons', comparison['passed'] == comparison['total'] == 261)
    for path, expected in independent['hashes'].items():
        check('independent_input_hash:' + path, sha(REPO / path) == expected)
    for path, expected in comparison['production_output_hashes'].items():
        check('independent_output_hash:' + path, sha(ROOT / 'results' / path) == expected)
    for stem, expected_count in [('matching_operator_audit', 26), ('matching_gauge_rg', 34),
                                  ('matching_mixed_sector', 122)]:
        base = REPO / 'reports' / (stem + '_20260925_checks')
        record = json.loads(base.with_suffix('.json').read_text())
        passed = record.get('passed', record.get('checks_passed'))
        total = record.get('total', record.get('checks_total'))
        check('topic_checks:' + stem, passed == total == expected_count)
        check('topic_script_hash:' + stem, sha(base.with_suffix('.py')) == record['script_sha256'])
    adversarial = json.loads((REPO / 'reports/matching_gauge_rg_20260925_adversarial.json').read_text())
    check('adversarial_checks', adversarial['passed'] == adversarial['total'] == 10)
    for path, expected in adversarial['source_hashes'].items():
        check('adversarial_source_hash:' + path, sha(REPO / path) == expected)
    for name, expected in [('coefficient_budgets.csv', 12), ('response_kernels.csv', 3075),
                           ('rg_closure.csv', 192), ('independent_kernels.csv', 12300),
                           ('independent_direct_cw.csv', 1344)]:
        with (ROOT / 'results' / name).open() as stream:
            count = sum(1 for _ in csv.DictReader(stream))
        check('rows:' + name, count == expected, str(count))
    # Historical failed diagnostics intentionally remain; they are not current pass records.
    initial = json.loads((ROOT / 'results/independent_initial_diagnostic.json').read_text())
    check('initial_failure_preserved', initial['passed'] == 21 and initial['total'] == 22)
    unchanged = subprocess.run(['git', 'diff', '--quiet', '5f5c6f298dbc63e01057bb71d0f45f9d1a249604', '--',
                                'experiments/noscale_protection_v01', 'experiments/sugra_shift_v01',
                                'manuscript_cn'], cwd=REPO).returncode == 0
    check('parent_experiments_and_formal_manuscript_unchanged', unchanged)
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
