"""Independently reconstruct SNP chain mappings and gene intervals from source files.

No annotation implementation functions are imported. This verification does not
establish functional effects or sequence uniqueness outside the supplied chains.
"""
from pathlib import Path
from collections import defaultdict
import argparse
import csv
import gzip
import json
import re
from urllib.parse import unquote
from notebook_support import file_digest, read_json

ROOT = Path(__file__).resolve().parents[1]


def read_table(path):
    with path.open(encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def chain_positions(path, positions):
    """Reconstruct one-base positions from raw chain alignment blocks."""
    hits = defaultdict(list)
    queries = defaultdict(list)
    for key, chrom, pos in positions:
        queries[chrom].append((key, int(pos)-1))
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            parts = line.split()
            if not parts:
                continue
            if parts[0] == 'chain':
                if parts[4] != '+':
                    raise ValueError('Unexpected chain source strand')
                source_chrom, t = parts[2], int(parts[5])
                dest_chrom, qsize, strand, q = parts[7], int(parts[8]), parts[9], int(parts[10])
                relevant = [(key, p) for key, p in queries[source_chrom] if t <= p < int(parts[6])]
                continue
            size = int(parts[0])
            for key, p in relevant:
                if t <= p < t+size:
                    offset = q + p-t
                    dest_pos = offset+1 if strand == '+' else qsize-offset
                    hits[key].append((dest_chrom, dest_pos, dest_pos, strand))
            t += size
            q += size
            if len(parts) == 3:
                t += int(parts[1])
                q += int(parts[2])
    return hits


def gene_records(path, provider):
    genes = {}
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            if len(fields) != 9 or fields[2] not in {'gene', 'pseudogene'}:
                continue
            if provider == 'Ensembl':
                gene_id = re.search(r'gene_id "([^"]+)";', fields[8]).group(1)
            else:
                gene_id = unquote(next(x[3:] for x in fields[8].split(';') if x.startswith('ID=')))
            genes[gene_id] = (fields[0], int(fields[3]), int(fields[4]), fields[6])
    return genes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    config = read_json(ROOT / 'config/annotation.json')
    folder = args.output or ROOT / config['output_directory']
    resources = ROOT / config['resource_directory']
    receipt = read_json(folder / 'receipt.private.json')
    for relative, digest in receipt['outputs'].items():
        if file_digest(folder / relative) != digest:
            raise ValueError('Output receipt mismatch')
    variants = read_table(folder / 'variant_crosswalk.private.tsv')
    checks = 0
    for assembly, label in [('CanFam3.1', 'canFam3'), ('CanFam4', 'canFam4')]:
        sites = [(r['record_id'], r['source_chromosome'], r['source_position']) for r in variants if r['source_assembly'] == assembly]
        hits = chain_positions(resources / (label + 'ToROS.chain.gz'), sites)
        native = defaultdict(list)
        for r in read_bed(folder / 'mapping' / (label + '.mapped.bed')):
            native[r[3]].append((r[0], int(r[1])+1, int(r[2]), r[5]))
        for key, _, _ in sites:
            if sorted(hits[key]) != sorted(native[key]):
                raise ValueError(f'Independent forward chain reconstruction failed: {key}')
            checks += 1
    aliases = {r['chromosome']: r['genome_name'] for r in read_json(folder / 'assembly_compatibility.private.json')}
    links = read_table(folder / 'gene_annotations.private.tsv')
    for provider, filename in [('Ensembl', 'ensembl.gtf.gz'), ('RefSeq', 'refseq.gff.gz')]:
        genes = gene_records(resources / filename, provider)
        for link in links:
            if link['provider'] != provider:
                continue
            chrom, start, end, strand = genes[link['gene_id']]
            if provider == 'Ensembl':
                chrom = aliases[chrom]
            if (chrom, start, end, strand) != (link['chromosome'], int(link['start']), int(link['end']), link['strand']):
                raise ValueError('Independent provider gene-interval check failed')
            pos = int(link['ROS_position'])
            dist = start-pos if pos < start else pos-end if pos > end else 0
            if dist != int(link['distance_bp']) or (dist == 0) != (link['relation'] == 'overlap'):
                raise ValueError('Independent gene-distance check failed')
    report = dict(status='passed', independently_reconstructed_forward_SNP_mappings=checks,
                  independently_checked_gene_links=len(links), output_receipts_checked=len(receipt['outputs']),
                  verifier_sha256=file_digest(Path(__file__)),
                  annotation_receipt_sha256=file_digest(folder / 'receipt.private.json'))
    (folder / 'verification.private.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report, indent=2))


def read_bed(path):
    return [line.split() for line in path.read_text().splitlines() if line.strip()]


if __name__ == '__main__':
    main()
