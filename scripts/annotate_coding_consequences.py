"""Calculate reference-relative SNP codon changes for CDS-overlap transcripts.

Both providers remain separate. This small, independently checkable calculation
uses the ROS genome and source CDS phases; it is not a general VEP replacement.
It reports transcript predictions, not disease effects or association direction.
"""
from pathlib import Path
from collections import defaultdict
import argparse
import json
from annotate_candidates import table, write_table, feature_rows, save_json
from notebook_support import file_digest, read_json
from reference_alignment import TwoBit, COMP

ROOT = Path(__file__).resolve().parents[1]
BASES = 'TCAG'
AMINO = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODE = {a+b+c: amino for (a, b, c), amino in zip(((a, b, c) for a in BASES for b in BASES for c in BASES), AMINO)}


def consequence(reference, alternate):
    a, b = CODE.get(reference, '?'), CODE.get(alternate, '?')
    if '?' in (a, b):
        return a, b, 'unresolved_codon'
    label = 'synonymous_variant' if a == b else 'stop_gained' if b == '*' else 'stop_lost' if a == '*' else 'missense_variant'
    return a, b, label


def coding_sequence(genome, chromosome, strand, blocks, position):
    blocks = sorted(set(blocks), reverse=strand == '-')
    first_phase = blocks[0][2]
    sequence = ''
    offset = None
    for i, (start, end, phase) in enumerate(blocks):
        expected = first_phase if i == 0 else (3-(len(sequence)-first_phase) % 3) % 3
        if phase != expected:
            raise ValueError('Inconsistent CDS phases')
        part = genome.get(chromosome, start-1, end)
        if len(part) != end-start+1:
            raise ValueError('CDS lies outside reference sequence')
        if strand == '-':
            part = part.translate(COMP)[::-1]
        if start <= position <= end:
            if offset is not None:
                raise ValueError('Overlapping CDS blocks')
            offset = len(sequence) + (position-start if strand == '+' else end-position) - first_phase
        sequence += part
    return sequence[first_phase:], offset


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--annotation', type=Path)
    args = parser.parse_args()
    config = read_json(ROOT / 'config/annotation.json')
    folder = args.annotation or ROOT / config['output_directory']
    resources = ROOT / config['resource_directory']
    output = folder / 'coding'
    receipt_path = folder / 'receipt.private.json'
    receipt = read_json(receipt_path)
    source_paths = [Path(__file__), ROOT / 'scripts/annotate_candidates.py', ROOT / 'scripts/reference_alignment.py',
                    ROOT / 'scripts/notebook_support.py', resources / 'ensembl.gtf.gz', resources / 'refseq.gff.gz',
                    resources / 'ros.2bit', folder / 'variant_crosswalk.private.tsv', folder / 'transcript_annotations.private.tsv',
                    folder / 'assembly_compatibility.private.json', receipt_path]
    sources = {str(p.relative_to(ROOT)): file_digest(p) for p in source_paths}
    for name in ['variant_crosswalk.private.tsv', 'transcript_annotations.private.tsv', 'assembly_compatibility.private.json']:
        if file_digest(folder / name) != receipt['outputs'][name]:
            raise ValueError('Annotation stage changed')
    if (output / 'receipt.private.json').exists():
        saved = read_json(output / 'receipt.private.json')
        if saved['sources'] != sources or file_digest(output / 'coding_consequences.private.tsv') != saved['output_sha256']:
            raise ValueError('Coding inputs/output changed; use a new annotation folder')
        print('Verified completed CDS predictions.')
        return
    if output.exists():
        raise FileExistsError('Incomplete coding folder; use a new annotation folder')
    genome = TwoBit(resources / 'ros.2bit')
    aliases = {}
    for r in read_json(folder / 'assembly_compatibility.private.json'):
        for key in [r['chromosome'], 'chr'+r['chromosome'], r['genbank'], r['refseq'], r['genome_name']]:
            aliases[key] = r['genome_name']
    variants = {r['record_id']: r for r in table(folder / 'variant_crosswalk.private.tsv')}
    selected = [r for r in table(folder / 'transcript_annotations.private.tsv') if r['region_annotation'] == 'CDS_overlap']
    results = []
    for provider, filename in [('Ensembl', 'ensembl.gtf.gz'), ('RefSeq', 'refseq.gff.gz')]:
        needed = {r['transcript_id'] for r in selected if r['provider'] == provider}
        blocks = defaultdict(list)
        exceptional = set()
        # Read source phase directly; the interval-stage iterator intentionally
        # omits phase because it performs no coding calculations.
        import gzip
        from annotate_candidates import attributes
        with gzip.open(resources / filename, 'rt') as handle:
            for line in handle:
                if line.startswith('#'):
                    continue
                fields = line.rstrip().split('\t')
                if len(fields) != 9 or fields[0] not in aliases:
                    continue
                a = attributes(fields[8], provider)
                tids = [a.get('transcript_id', '')] if provider == 'Ensembl' else a.get('Parent', '').split(',')
                if provider == 'RefSeq' and a.get('ID') in needed and ('exception' in a or 'transl_except' in a):
                    exceptional.add(a['ID'])
                for tid in set(tids) & needed:
                    if 'exception' in a or 'transl_except' in a:
                        exceptional.add(tid)
                    if fields[2] == 'CDS':
                        blocks[tid].append((int(fields[3]), int(fields[4]), int(fields[7])))
        for row in selected:
            if row['provider'] != provider:
                continue
            variant = variants[row['record_id']]
            r = dict(record_id=row['record_id'], marker=variant['marker'], provider=provider, release=row['release'],
                     gene_id=row['gene_id'], symbol=row['symbol'], transcript_id=row['transcript_id'], strand=row['strand'],
                     ROS_reference=variant['ROS_reference'], ROS_other_allele=variant['ROS_other_allele'],
                     prediction_status='not_calculated', coding_position='', amino_acid_position='',
                     reference_codon='', alternate_codon='', reference_amino_acid='', alternate_amino_acid='',
                     consequence='', note='Reference-relative prediction; no historical/risk allele inferred.')
            tid = row['transcript_id']
            if tid in exceptional or not blocks[tid]:
                r['prediction_status'] = 'exceptional_or_missing_CDS'
            else:
                try:
                    sequence, offset = coding_sequence(genome, variant['ROS_chromosome'], row['strand'], blocks[tid], int(variant['ROS_position']))
                    if offset is None or offset < 0:
                        raise ValueError('Variant is in an incomplete start codon')
                    ref, alt = variant['ROS_reference'], variant['ROS_other_allele']
                    if row['strand'] == '-':
                        ref, alt = ref.translate(COMP), alt.translate(COMP)
                    if ref not in 'ACGT' or alt not in 'ACGT' or len(alt) != 1 or sequence[offset] != ref:
                        raise ValueError('Transcript reference or allele is unresolved')
                    start = (offset // 3) * 3
                    codon = sequence[start:start+3]
                    if len(codon) != 3:
                        raise ValueError('Variant is in an incomplete end codon')
                    alternate = codon[:offset % 3] + alt + codon[offset % 3+1:]
                    aa, aa_alt, label = consequence(codon, alternate)
                    if CODE.get(sequence[:3]) != 'M' and start == 0:
                        raise ValueError('Uncertain initiation codon')
                    if start == 0 and codon == 'ATG' and alternate != 'ATG':
                        label = 'start_lost'
                    r.update(prediction_status='calculated', coding_position=offset+1, amino_acid_position=offset//3+1,
                             reference_codon=codon, alternate_codon=alternate, reference_amino_acid=aa,
                             alternate_amino_acid=aa_alt, consequence=label)
                except (ValueError, IndexError) as error:
                    r.update(prediction_status='unresolved', note=str(error))
            results.append(r)
    output.mkdir()
    fields = ['record_id', 'marker', 'provider', 'release', 'gene_id', 'symbol', 'transcript_id', 'strand',
              'ROS_reference', 'ROS_other_allele', 'prediction_status', 'coding_position', 'amino_acid_position',
              'reference_codon', 'alternate_codon', 'reference_amino_acid', 'alternate_amino_acid', 'consequence', 'note']
    write_table(output / 'coding_consequences.private.tsv', results, fields)
    save_json(output / 'receipt.private.json', dict(status='complete', sources=sources, transcript_predictions=len(results),
                                                   output_sha256=file_digest(output / 'coding_consequences.private.tsv')))
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
