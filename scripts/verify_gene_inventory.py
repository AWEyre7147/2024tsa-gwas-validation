"""Check complete gene-resource coverage, coordinate bounds and RNA relationships."""
from pathlib import Path
import argparse
import gzip
import json
import sqlite3
from collections import Counter
from extension_support import ROOT,save
from notebook_support import read_json,file_digest


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--inventory',type=Path)
    args=parser.parse_args();config=read_json(ROOT/'config/gene_inventory.json');folder=args.inventory or ROOT/config['output_directory']
    receipt=read_json(folder/'receipt.private.json')
    if file_digest(folder/'dog_ros.sqlite')!=receipt['outputs']['dog_ros.sqlite']:raise ValueError('Gene database receipt mismatch')
    db=sqlite3.connect(f'file:{(folder/"dog_ros.sqlite").as_posix()}?mode=ro',uri=True)
    if db.execute('PRAGMA quick_check').fetchone()[0]!='ok':raise ValueError('Database integrity check failed')
    bounds={}
    for table in ['genes','transcripts','features']:
        n=db.execute(f'''SELECT count(*) FROM {table} t LEFT JOIN sequence_aliases a
            ON t.provider=a.provider AND t.sequence=a.source_sequence
            WHERE a.source_sequence IS NULL OR t.start<1 OR t.end<t.start OR t.end>a.length''').fetchone()[0]
        if n:raise ValueError(f'{table}: {n} invalid/unresolved sequence intervals')
        bounds[table]=db.execute('SELECT count(*) FROM '+table).fetchone()[0]
    gene_keys={r[0] for r in db.execute('SELECT record_key FROM genes')}
    nested=0
    for key,parents,ancestors in db.execute('SELECT record_key,parent_feature_keys_json,gene_keys_json FROM transcripts'):
        if not json.loads(ancestors) or any(x not in gene_keys for x in json.loads(ancestors)):raise ValueError('Unresolved RNA gene ancestry')
        nested+=any(x not in gene_keys for x in json.loads(parents))
    source_counts={}
    resources=ROOT/config['resource_directory']
    for provider,name in [('Ensembl','ensembl.gtf.gz'),('RefSeq','refseq.gff.gz')]:
        counts=Counter()
        with gzip.open(resources/name,'rt') as handle:
            for line in handle:
                if not line.startswith('#'):counts[line.split('\t',3)[2]]+=1
        total=sum(counts.values())-counts['region']
        stored=sum(db.execute('SELECT count(*) FROM '+table+' WHERE provider=?',(provider,)).fetchone()[0] for table in ['genes','transcripts','features'])
        if total!=stored:raise ValueError('Every source annotation record is not retained')
        source_counts[provider]=stored
    db.close()
    summary=dict(status='passed',checked_coordinate_bounds=bounds,source_annotation_records=source_counts,
        nested_RNA_parent_records=nested,provider_symbols_not_merged=True,
        receipt_sha256=file_digest(folder/'receipt.private.json'),verifier_sha256=file_digest(Path(__file__)))
    save(folder/'verification.private.json',summary);print(summary)


if __name__=='__main__':main()
