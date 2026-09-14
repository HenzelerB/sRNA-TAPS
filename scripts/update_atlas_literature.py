#!/usr/bin/env python3
"""Attach position-matched citations from the curated miRNA literature master TSV.

Run with --master ../bio/report/human_mature_miRNA_m5C_reported_sites_master.tsv
from the repository root. Nonpositional miRNA associations cannot support a site.
"""
import argparse
import csv
import json
import re
from pathlib import Path

DATA = re.compile(r'(<script id="site-data"[^>]*>)(.*?)(</script>)', re.S)
POSITION = re.compile(r'nt (\d+) \(mature\)')
RESOLUTIONS = {'exact_mature_nucleotide', 'CpG_candidate_position_set'}


def annotate(rows, sources):
    index = {}
    for source in sources:
        if source['evidence_resolution'] in RESOLUTIONS:
            key = (source['miRNA_current'], source['mature_position_1based'])
            index.setdefault(key, []).append(source)
    count = 0
    for row in rows:
        if row['biotype'] != 'miRNA' or row['status'] != 'reported':
            continue
        position = POSITION.fullmatch(row['position_label'])
        matches = index.get((row['transcript'], position[1]), []) if position else []
        if not matches:
            raise ValueError('No positional literature for {} {}'.format(row['transcript'], row['position_label']))
        publications = {}
        for source in matches:
            key = source['doi'] or source['pmid'] or source['title']
            if not key:
                raise ValueError('Missing publication identity')
            if key not in publications:
                publications[key] = {
                    'title': source['title'], 'first_author': source['first_author'],
                    'year': source['publication_year'], 'doi': source['doi'],
                    'pmid': source['pmid'], 'pmcid': source['pmcid'],
                    'url': source['source_url'], 'evidence': [],
                }
            evidence = {field: source[field] for field in (
                'record_id', 'miRNA_current', 'mature_position_1based',
                'evidence_resolution', 'evidence_class', 'modification_reported',
                'cell_line_or_tissue', 'source_table_or_figure', 'notes')}
            if evidence not in publications[key]['evidence']:
                publications[key]['evidence'].append(evidence)
        row['literature'] = list(publications.values())
        count += 1
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--master', type=Path, required=True)
    parser.add_argument('--atlas', type=Path, default=Path(__file__).resolve().parents[1] / 'docs/m5c-atlas.html')
    args = parser.parse_args()
    html = args.atlas.read_text()
    match = DATA.search(html)
    if not match:
        raise ValueError('Atlas site-data block missing')
    rows = json.loads(match[2])
    with args.master.open(newline='') as handle:
        count = annotate(rows, list(csv.DictReader(handle, delimiter='\t')))
    data = json.dumps(rows, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    args.atlas.write_text(html[:match.start(2)] + data + html[match.end(2):])
    print('Updated positional literature for {} reported miRNA sites.'.format(count))


if __name__ == '__main__':
    main()
