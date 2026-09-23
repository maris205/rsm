"""Split the supplied RSM manuscript into faithful, auditable Markdown chapters.

Requires pandoc on PATH. The source DOCX is never edited. Scholarly review
and corrections belong in separate reports, not in this conversion.
"""
from pathlib import Path
from collections import Counter
import copy
import hashlib
import json
import re
import subprocess
import zipfile
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/'ori_paper/merged_document_v2_with_latex.docx'
CHAPTERS = ROOT/'chapters'
FORMAT = 'markdown+tex_math_dollars-simple_tables-multiline_tables-grid_tables'
SLUGS = ['introduction', 'system_anomaly', 'source_code', 'cosmic_lockstep',
         'mathematical_framework', 'observational_evidence', 'predictions_interpretations',
         'engine_design', 'simulation_case_studies', 'conclusion']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def walk(node):
    if isinstance(node, dict):
        yield node
        for value in node.values():
            yield from walk(value)
    elif isinstance(node, list):
        for value in node:
            yield from walk(value)


def inline_text(node):
    if isinstance(node, dict):
        if node.get('t') == 'Str':
            return node['c']
        if node.get('t') in ('Space', 'SoftBreak', 'LineBreak'):
            return ' '
        if node.get('t') == 'Math':
            return node['c'][1]
        return inline_text(node.get('c', []))
    if isinstance(node, list):
        return ''.join(inline_text(value) for value in node)
    return ''


def counts(ast):
    counter = Counter(n.get('t') for n in walk(ast) if 't' in n)
    return {key: counter[key] for key in ('Header', 'Math', 'InlineMath', 'DisplayMath', 'Table', 'Image', 'Link', 'Note')}


def content_signatures(ast):
    math_order = [n['c'][1] for n in walk(ast) if n.get('t') == 'Math']
    word_order = [word for n in walk(ast) if n.get('t') == 'Str'
                  for word in re.findall(r'\w+', n['c'], flags=re.UNICODE)]
    math, words = Counter(math_order), Counter(word_order)
    images = Counter(n['c'][-1][0] for n in walk(ast) if n.get('t') == 'Image')
    return math, words, images, math_order, word_order


def call(args, **kwargs):
    completed = subprocess.run(args, text=True, encoding='utf-8', stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, check=True, **kwargs)
    if completed.stderr:
        LOG.append(completed.stderr)
    return completed.stdout


LOG = []


def main():
    for directory in (CHAPTERS, ROOT/'reports', ROOT/'build'):
        directory.mkdir(parents=True, exist_ok=True)
    original_hash = sha(SOURCE)
    version = call(['pandoc', '--version']).splitlines()[0]
    source_ast = json.loads(call(['pandoc', str(SOURCE), '-f', 'docx', '-t', 'json',
                                 '--extract-media=assets'], cwd=CHAPTERS))
    (ROOT/'build/source_ast.json').write_text(json.dumps(source_ast, ensure_ascii=False), encoding='utf-8')
    blocks = source_ast['blocks']
    boundaries = [(0, '00_front_matter.md', 'Front matter: title, abstract and original preliminary material')]
    for index, block in enumerate(blocks):
        if block['t'] != 'Header':
            continue
        title = inline_text(block['c'][2])
        match = re.match(r'^Chapter\s+(\d+):', title)
        if match:
            number = int(match.group(1))
            boundaries.append((index, f'{number+1:02d}_chapter_{number:02d}_{SLUGS[number]}.md', title))
        elif block['c'][0] == 2 and title in ('References', 'Acknowledgments'):
            prefix = '11' if title == 'References' else '12'
            boundaries.append((index, f'{prefix}_{title.lower()}.md', title))
    assert len(boundaries) == 13 and [b[0] for b in boundaries] == sorted(b[0] for b in boundaries)
    records, parsed_blocks = [], []
    for order, (start, name, title) in enumerate(boundaries):
        stop = boundaries[order+1][0] if order+1 < len(boundaries) else len(blocks)
        piece = copy.deepcopy(source_ast)
        piece['blocks'] = copy.deepcopy(blocks[start:stop])
        # The original Word chapter headers are Heading 2; promote all source
        # headers by one level, preserving their hierarchy inside each file.
        for node in walk(piece['blocks']):
            if node.get('t') == 'Header':
                node['c'][0] = max(1, node['c'][0]-1)
                node['c'][1] = ['', [], []]
            elif node.get('t') == 'Image':
                node['c'][0] = ['', [], []]
        markdown = call(['pandoc', '-f', 'json', '-t', FORMAT, '--wrap=none', '--atx-headers'],
                        input=json.dumps(piece, ensure_ascii=False), cwd=CHAPTERS)
        comment = (f'<!-- Faithful format conversion; scientific claims are unreviewed source text. '
                   f'Source: ../ori_paper/{SOURCE.name}; Pandoc blocks [{start}, {stop}). -->\n\n')
        output = CHAPTERS/name
        output.write_text(comment+markdown, encoding='utf-8')
        reread = json.loads(call(['pandoc', str(output), '-f', FORMAT, '-t', 'json'], cwd=CHAPTERS))
        before, after = content_signatures(piece), content_signatures(reread)
        record = {'file': str(output.relative_to(ROOT)), 'title': title,
                  'source_block_start_inclusive': start, 'source_block_stop_exclusive': stop,
                  'source_counts': counts(piece), 'markdown_counts': counts(reread),
                  'math_exact_multiset_preserved': before[0] == after[0],
                  'prose_word_multiset_preserved': before[1] == after[1],
                  'image_targets_preserved': before[2] == after[2],
                  'equation_order_preserved': before[3] == after[3],
                  'prose_word_order_preserved': before[4] == after[4],
                  'word_count': sum(before[1].values()), 'sha256': sha(output)}
        if before[1] != after[1]:
            record['prose_word_difference'] = {'missing': dict(before[1]-after[1]), 'added': dict(after[1]-before[1])}
        if before[0] != after[0]:
            record['math_difference'] = {'missing': dict(before[0]-after[0]), 'added': dict(after[0]-before[0])}
        records.append(record)
        parsed_blocks.extend(reread['blocks'])
    xml_ns = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
              'm': 'http://schemas.openxmlformats.org/officeDocument/2006/math',
              'a': 'http://schemas.openxmlformats.org/drawingml/2006/main'}
    with zipfile.ZipFile(SOURCE) as archive:
        xml = ET.fromstring(archive.read('word/document.xml'))
        media = [name for name in archive.namelist() if name.startswith('word/media/')]
        source_inventory = {'omath_nodes': len(xml.findall('.//m:oMath', xml_ns)),
                            'tables': len(xml.findall('.//w:tbl', xml_ns)),
                            'body_image_references': len(xml.findall('.//a:blip', xml_ns)),
                            'embedded_media_entries': len(media),
                            'tracked_insertions': len(xml.findall('.//w:ins', xml_ns)),
                            'tracked_deletions': len(xml.findall('.//w:del', xml_ns))}
        assets = []
        for image in sorted((CHAPTERS/'assets').rglob('*')):
            if not image.is_file():
                continue
            relative = str(image.relative_to(CHAPTERS/'assets'))
            member = 'word/'+relative
            original = archive.read(member)
            assert image.read_bytes() == original
            assets.append({'file': str(image.relative_to(ROOT)), 'docx_member': member,
                           'bytes': len(original), 'sha256': sha(image)})
    combined = dict(source_ast, blocks=parsed_blocks)
    source_counts, output_counts = counts(source_ast), counts(combined)
    assert sha(SOURCE) == original_hash
    partition_ok = records[0]['source_block_start_inclusive'] == 0 and records[-1]['source_block_stop_exclusive'] == len(blocks)
    checks = {'source_unchanged': True, 'all_blocks_partitioned_once': partition_ok,
              'source_omath_to_pandoc_math': source_inventory['omath_nodes'] == source_counts['Math'],
              'equations_preserved': all(r['math_exact_multiset_preserved'] for r in records),
              'prose_words_preserved': all(r['prose_word_multiset_preserved'] for r in records),
              'equation_order_preserved': all(r['equation_order_preserved'] for r in records),
              'prose_word_order_preserved': all(r['prose_word_order_preserved'] for r in records),
              'table_count_preserved': source_counts['Table'] == output_counts['Table'] == source_inventory['tables'],
              'image_count_preserved': source_counts['Image'] == output_counts['Image'] == source_inventory['body_image_references'],
              'image_targets_preserved': all(r['image_targets_preserved'] for r in records),
              'header_count_preserved': source_counts['Header'] == output_counts['Header'],
              'embedded_images_byte_identical': True}
    report = {'source': str(SOURCE.relative_to(ROOT)), 'source_sha256': original_hash,
              'pandoc_version': version, 'markdown_format': FORMAT, 'chapter_file_count': len(records),
              'source_top_level_blocks': len(blocks), 'source_docx_inventory': source_inventory,
              'source_ast_counts': source_counts, 'markdown_roundtrip_counts': output_counts,
              'checks': checks, 'status': 'PASS' if all(checks.values()) else 'REVIEW',
              'chapters': records, 'assets': assets,
              'notes': ['No scientific claims, citations or numerical values were revised or independently verified.',
                        'Original preliminary material and repeated comparison appendices remain in their source positions.',
                        'Word page layout, image sizing and typography are not preserved; heading levels are promoted by one.',
                        'Original Word TOC labels/page numbers are preserved in the front matter; use the new chapter index for working links.',
                        'Only body-referenced media are displayed; package-only media are not inserted as new content.'],
              'conversion_script_sha256': sha(__file__)}
    (ROOT/'reports/conversion_manifest.json').write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    (ROOT/'reports/conversion_log.txt').write_text(''.join(LOG) or 'Pandoc emitted no conversion warnings.\n', encoding='utf-8')
    print(json.dumps({'status': report['status'], 'checks': checks, 'chapters': len(records),
                      'counts': source_counts}, ensure_ascii=False, indent=2))
    if report['status'] != 'PASS':
        raise SystemExit('Conversion requires review; see reports/conversion_manifest.json')


if __name__ == '__main__':
    main()
