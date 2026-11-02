"""Check predicted codons against independently distributed provider CDS FASTA.

Downloads a small validation cache, verifies published checksums, and matches
transcript/protein IDs. Neither annotation stage nor genome-derived calculations
are rerun or changed.
"""
from pathlib import Path
from collections import defaultdict
import argparse
import gzip
import json
import re
import subprocess
from notebook_support import read_json, file_digest
from download_annotation_assets import published_md5

ROOT = Path(__file__).resolve().parents[1]
NCBI = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9615/106/GCF_014441545.1_ROS_Cfam_1.0/'
ENSEMBL = 'https://ftp.ensembl.org/pub/release-116/fasta/canis_lupus_familiaris/cds/'
CODON_GROUPS = {
    'F': 'TTT TTC', 'L': 'TTA TTG CTT CTC CTA CTG', 'I': 'ATT ATC ATA', 'M': 'ATG',
    'V': 'GTT GTC GTA GTG', 'S': 'TCT TCC TCA TCG AGT AGC', 'P': 'CCT CCC CCA CCG',
    'T': 'ACT ACC ACA ACG', 'A': 'GCT GCC GCA GCG', 'Y': 'TAT TAC', '*': 'TAA TAG TGA',
    'H': 'CAT CAC', 'Q': 'CAA CAG', 'N': 'AAT AAC', 'K': 'AAA AAG', 'D': 'GAT GAC',
    'E': 'GAA GAG', 'C': 'TGT TGC', 'W': 'TGG', 'R': 'CGT CGC CGA CGG AGA AGG', 'G': 'GGT GGC GGA GGG'
}


def fasta(path):
    with gzip.open(path, 'rt') as handle:
        header = None
        parts = []
        for line in handle:
            if line.startswith('>'):
                if header:
                    yield header, ''.join(parts)
                header, parts = line[1:].strip(), []
            else:
                parts.append(line.strip())
        if header:
            yield header, ''.join(parts)


def main():
    import csv
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--annotation', type=Path)
    args = parser.parse_args()
    config = read_json(ROOT / 'config/annotation.json')
    folder = args.annotation or ROOT / config['output_directory']
    resources = ROOT / config['resource_directory']
    cache = resources / 'coding-validation'
    cache.mkdir(exist_ok=True)
    urls = {'refseq.cds.fa.gz': NCBI + 'GCF_014441545.1_ROS_Cfam_1.0_cds_from_genomic.fna.gz',
            'ensembl.cds.fa.gz': ENSEMBL + 'Canis_lupus_familiaris.ROS_Cfam_1.0.cds.all.fa.gz',
            'ensembl.CHECKSUMS': ENSEMBL + 'CHECKSUMS'}
    lock = cache / 'assets.json'
    if lock.exists():
        saved = read_json(lock)
        if saved['urls'] != urls or any(file_digest(cache / name) != digest for name, digest in saved['files'].items()):
            raise ValueError('CDS validation cache differs')
    else:
        for name, url in urls.items():
            path = cache / name
            if not path.exists():
                partial = path.with_suffix(path.suffix + '.partial')
                subprocess.run(['curl', '--fail', '--location', '--silent', '--show-error', '--retry', '3',
                                '--connect-timeout', '30', '--max-time', '300', url, '-o', str(partial)], check=True)
                partial.replace(path)
        published_md5(cache / 'refseq.cds.fa.gz', resources / 'ncbi.md5checksums.txt',
                      'GCF_014441545.1_ROS_Cfam_1.0_cds_from_genomic.fna.gz')
        line = next(x for x in (cache / 'ensembl.CHECKSUMS').read_text().splitlines()
                    if x.split()[-1] == 'Canis_lupus_familiaris.ROS_Cfam_1.0.cds.all.fa.gz')
        actual = subprocess.check_output(['sum', '-r', str(cache / 'ensembl.cds.fa.gz')], text=True)
        if [int(x) for x in line.split()[:2]] != [int(x) for x in actual.split()[:2]]:
            raise ValueError('Ensembl CDS checksum mismatch')
        lock.write_text(json.dumps(dict(urls=urls, files={name: file_digest(cache / name) for name in urls}), indent=2) + '\n')
    path = folder / 'coding/coding_consequences.private.tsv'
    with path.open(newline='') as handle:
        predictions = list(csv.DictReader(handle, delimiter='\t'))
    protein_ids = defaultdict(set)
    with gzip.open(resources / 'refseq.gff.gz', 'rt') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            if len(fields) != 9 or fields[2] != 'CDS':
                continue
            a = dict(x.split('=', 1) for x in fields[8].split(';') if '=' in x)
            for tid in a.get('Parent', '').split(','):
                if a.get('protein_id'):
                    protein_ids[tid].add(a['protein_id'])
    sequences = {'Ensembl': defaultdict(set), 'RefSeq': defaultdict(set)}
    wanted = {'Ensembl': {r['transcript_id'] for r in predictions if r['provider'] == 'Ensembl'},
              'RefSeq': set().union(*(protein_ids[r['transcript_id']] for r in predictions if r['provider'] == 'RefSeq'))}
    for header, seq in fasta(cache / 'ensembl.cds.fa.gz'):
        key = header.split()[0].split('.')[0]
        if key in wanted['Ensembl']:
            sequences['Ensembl'][key].add(seq.upper())
    for header, seq in fasta(cache / 'refseq.cds.fa.gz'):
        m = re.search(r'\[protein_id=([^\]]+)\]', header)
        if m and m[1] in wanted['RefSeq']:
            sequences['RefSeq'][m[1]].add(seq.upper())
    checks = []
    for r in predictions:
        if r['prediction_status'] != 'calculated':
            continue
        keys = {r['transcript_id']} if r['provider'] == 'Ensembl' else protein_ids[r['transcript_id']]
        options = set().union(*(sequences[r['provider']][key] for key in keys))
        if len(options) != 1:
            raise ValueError('Provider CDS sequence missing or ambiguous')
        seq = next(iter(options))
        offset = int(r['coding_position'])-1
        start = (offset//3)*3
        if seq[start:start+3] != r['reference_codon']:
            raise ValueError(f'Independent provider CDS codon mismatch: {r["transcript_id"]}')
        for prefix in ['reference', 'alternate']:
            if r[prefix + '_codon'] not in CODON_GROUPS.get(r[prefix + '_amino_acid'], '').split():
                raise ValueError('Independent amino-acid lookup disagrees with prediction')
        checks.append(dict(provider=r['provider'], transcript_id=r['transcript_id'],
                           CDS_ids=sorted(keys), reference_codon=seq[start:start+3], status='source_CDS_codon_matches'))
    report = dict(status='passed', provider_CDS_checks=checks, predictions_sha256=file_digest(path),
                  source_CDS_cache_sha256=file_digest(lock), verifier_sha256=file_digest(Path(__file__)))
    (folder / 'coding/verification.private.json').write_text(json.dumps(report, indent=2) + '\n')
    print(f'Provider CDS sequences independently confirm {len(checks)} transcript codons.')


if __name__ == '__main__':
    main()
