"""Attach existing MyGene functional context by exact stable IDs, never position/symbol.

Requires pyarrow only for this optional stage. The source Parquet is read-only;
no expanded or refreshed gene-annotation resource is created here.
"""
from pathlib import Path
from collections import defaultdict
import argparse
import csv
import json
import pyarrow
import pyarrow.parquet as pq
from notebook_support import file_digest, read_json

ROOT = Path(__file__).resolve().parents[1]


def table(path):
    with path.open(encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def write_table(path, records):
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mygene', type=Path, required=True)
    parser.add_argument('--annotation', type=Path)
    args = parser.parse_args()
    config = read_json(ROOT / 'config/annotation.json')
    folder = args.annotation or ROOT / config['output_directory']
    source = args.mygene.resolve()
    output = folder / 'mygene'
    sources = dict(parquet_sha256=file_digest(source), gene_links_sha256=file_digest(folder / 'gene_annotations.private.tsv'),
                   script_sha256=file_digest(Path(__file__)), annotation_script_sha256=file_digest(ROOT / 'scripts/annotate_candidates.py'),
                   notebook_support_sha256=file_digest(ROOT / 'scripts/notebook_support.py'), pyarrow_version=pyarrow.__version__)
    if (output / 'receipt.private.json').exists():
        saved = read_json(output / 'receipt.private.json')
        if saved['sources'] != sources:
            raise ValueError('Supplement inputs/code changed; choose a new annotation folder')
        for name, digest in saved['outputs'].items():
            if file_digest(output / name) != digest:
                raise ValueError('Supplement output changed')
        print('Verified completed MyGene supplement.')
        return
    if output.exists():
        raise FileExistsError('Incomplete/existing MyGene supplement folder')
    schema = pq.read_schema(source)
    requested = ['symbol', 'name', 'mygene_id', 'entrezgene', 'ensembl.gene', 'type_of_gene', 'summary',
                 'alias', 'other_names', 'go.BP', 'go.CC', 'go.MF', 'interpro', 'pfam', 'homologene',
                 'refseq.rna', 'uniprot.Swiss-Prot', 'uniprot.TrEMBL', 'pathway.reactome']
    rows = pq.read_table(source, columns=[name for name in requested if name in schema.names]).to_pylist()
    entrez = defaultdict(list)
    ensembl = defaultdict(list)
    for r in rows:
        if r.get('entrezgene'):
            entrez[str(r['entrezgene'])].append(r)
        ids = r.get('ensembl.gene') or []
        if isinstance(ids, str):
            ids = [ids]
        for key in ids:
            ensembl[key].append(r)
    links = table(folder / 'gene_annotations.private.tsv')
    genes = {(r['provider'], r['gene_id']): r for r in links}
    matching = []
    records = []
    for (provider, gid), g in sorted(genes.items()):
        candidates = entrez[g['entrezgene']] if g['entrezgene'] else ensembl[gid]
        candidates = list({str(r['mygene_id']): r for r in candidates}.values())
        status = 'unique_exact_ID' if len(candidates) == 1 else 'ambiguous_ID' if candidates else 'no_exact_ID_match'
        matching.append(dict(provider=provider, gene_id=gid, annotation_symbol=g['symbol'], annotation_biotype=g['biotype'],
                             entrezgene=g['entrezgene'], match_status=status,
                             mygene_ids=';'.join(str(r['mygene_id']) for r in candidates),
                             matching_basis='entrezgene' if g['entrezgene'] else 'ensembl_gene_id'))
        if len(candidates) == 1:
            records.append(dict(provider=provider, annotation_gene_id=gid, functional_context=candidates[0]))
    output.mkdir()
    write_table(output / 'matching.private.tsv', matching)
    save_json(output / 'functional_context.private.json', records)
    summary = dict(source_rows=len(rows), distinct_provider_genes=len(genes), exact_matches=len(records),
                   unmatched=sum(r['match_status'] == 'no_exact_ID_match' for r in matching),
                   ambiguous=sum(r['match_status'] == 'ambiguous_ID' for r in matching),
                   source_gene_types=sorted({r.get('type_of_gene', '') for r in rows}),
                   limitations=['Existing protein-coding resource; absence of a match is not biological absence.',
                                'No symbol or coordinate matching; Ensembl ID changes can leave records unmatched.',
                                'Functional annotations are source-provided context, not evidence for a causal variant.'])
    save_json(output / 'summary.private.json', summary)
    save_json(output / 'receipt.private.json', dict(status='complete', source_path=str(source), sources=sources,
                                                   outputs={p.name: file_digest(p) for p in output.iterdir() if p.is_file()}))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
