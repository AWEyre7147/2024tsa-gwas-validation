"""Build a reusable, complete ROS coding/noncoding gene and transcript inventory.

All provider feature relationships are retained in SQLite. No study dog data or
candidate selection enters this resource. Source feature IDs and stable gene IDs
remain separate; non-identical GenBank/RefSeq sequences are never joined by name.
"""
from pathlib import Path
from collections import Counter
import argparse
import csv
import gzip
import json
import re
import sqlite3
from urllib.parse import unquote
from extension_support import ROOT, save, sources, start, finish
from notebook_support import file_digest, read_json


def attrs(value, provider):
    if provider == 'Ensembl':
        return dict(re.findall(r'(\S+) "([^"]*)";', value))
    return {key: unquote(v) for item in value.split(';') if '=' in item for key,v in [item.split('=',1)]}


def gene_class(kind):
    if 'pseudogene' in kind or kind == 'pseudo':
        return 'pseudogene'
    if kind == 'protein_coding':
        return 'protein_coding'
    if kind.startswith(('IG_', 'TR_')) or kind in {'V_segment','C_region','D_segment','J_segment'}:
        return 'immune_gene_segment'
    if kind in {'lncRNA','lincRNA','miRNA','snoRNA','snRNA','rRNA','tRNA','misc_RNA','scaRNA',
                'Y_RNA','ribozyme','vault_RNA','guide_RNA','SRP_RNA','RNase_P_RNA','RNase_MRP_RNA','ncRNA'}:
        return 'noncoding'
    return 'other_or_unclassified'


def export_table(db, path, query):
    cursor = db.execute(query)
    with gzip.open(path, 'wt', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow([x[0] for x in cursor.description])
        writer.writerows(cursor)


def resolve_transcript_genes(db):
    """Preserve precursor/product parents and resolve their gene ancestors."""
    columns={r[1] for r in db.execute('PRAGMA table_info(transcripts)')}
    if 'parent_gene_keys_json' in columns:
        db.execute('ALTER TABLE transcripts RENAME COLUMN parent_gene_keys_json TO parent_feature_keys_json')
    if 'gene_keys_json' not in columns:
        db.execute('ALTER TABLE transcripts ADD COLUMN gene_keys_json TEXT')
    genes={r[0] for r in db.execute('SELECT record_key FROM genes')}
    parents={r[0]:json.loads(r[1]) for r in db.execute('SELECT record_key,parent_feature_keys_json FROM transcripts')}
    resolved={}
    def walk(key,stack):
        if key in genes:return {key}
        if key in resolved:return resolved[key]
        if key in stack:raise ValueError('Cyclic RNA parent relationship')
        if key not in parents:raise ValueError('Unknown RNA parent feature')
        result=set()
        for parent in parents[key]:result.update(walk(parent,stack|{key}))
        resolved[key]=result;return result
    for key in parents:walk(key,set())
    db.executemany('UPDATE transcripts SET gene_keys_json=? WHERE record_key=?',
                   [(json.dumps(sorted(values)),key) for key,values in resolved.items()])
    return sum(any(parent in parents for parent in values) for values in parents.values())


def resume_inventory(output,plan,resources):
    """Finalize a fully parsed but incomplete build after alias validation repair."""
    original=read_json(output/'plan.private.json')
    if original['config']!=plan['config']:
        raise ValueError('Incomplete inventory configuration differs')
    for relative,digest in original['sources'].items():
        if relative!='scripts/build_gene_inventory.py' and plan['sources'].get(relative)!=digest:
            raise ValueError('Incomplete inventory source changed')
    db=sqlite3.connect(output/'dog_ros.sqlite')
    if db.execute('PRAGMA integrity_check').fetchone()[0]!='ok':raise ValueError('Incomplete SQLite integrity check failed')
    headers={};stats={}
    for provider,name in [('Ensembl','ensembl.gtf.gz'),('RefSeq','refseq.gff.gz')]:
        counts=Counter();metadata=[]
        with gzip.open(resources/name,'rt') as handle:
            for line in handle:
                if line.startswith('#'):
                    if line.startswith('#!'):metadata.append(line.strip())
                else:counts[line.split('\t',3)[2]]+=1
        headers[provider]=metadata;stats[provider]=dict(counts)
        genes=db.execute('SELECT count(*) FROM genes WHERE provider=?',(provider,)).fetchone()[0]
        tx=db.execute('SELECT count(*) FROM transcripts WHERE provider=?',(provider,)).fetchone()[0]
        parts=db.execute('SELECT count(*) FROM features WHERE provider=?',(provider,)).fetchone()[0]
        if genes!=counts['gene']+counts['pseudogene'] or genes+tx+parts+counts['region']!=sum(counts.values()):
            raise ValueError('Incomplete inventory does not contain every source feature')
    # Assembly report names the main Y molecule by scaffold name, while Ensembl
    # uses assigned molecule Y. Only assembled-molecule aliases get this name;
    # the separate Y unlocalized scaffolds keep their own identifiers.
    main_y=db.execute('SELECT * FROM sequence_aliases WHERE provider="Ensembl" AND assigned_molecule="Y" AND role="assembled-molecule"').fetchall()
    if len(main_y)!=1:raise ValueError('Main Y sequence alias is ambiguous')
    r=main_y[0];db.execute('INSERT OR IGNORE INTO sequence_aliases VALUES (?,?,?,?,?,?,?,?,?)',(r[0],'Y',*r[2:]))
    gaps=db.execute('SELECT DISTINCT g.provider,g.sequence FROM genes g LEFT JOIN sequence_aliases a ON g.provider=a.provider AND g.sequence=a.source_sequence WHERE a.source_sequence IS NULL').fetchall()
    if gaps:raise ValueError('Unresolved source sequence aliases')
    nested=resolve_transcript_genes(db)
    db.commit()
    export_table(db,output/'genes.tsv.gz','SELECT * FROM genes ORDER BY provider,sequence,start,record_key')
    export_table(db,output/'transcripts.tsv.gz','SELECT * FROM transcripts ORDER BY provider,sequence,start,record_key')
    summary=dict(status='complete_gene_inventory',gene_counts=dict(db.execute('SELECT provider,count(*) FROM genes GROUP BY provider')),
        transcript_counts=dict(db.execute('SELECT provider,count(*) FROM transcripts GROUP BY provider')),
        gene_classes=[dict(provider=p,gene_class=c,count=n) for p,c,n in db.execute('SELECT provider,gene_class,count(*) FROM genes GROUP BY provider,gene_class')],
        feature_rows=db.execute('SELECT count(*) FROM features').fetchone()[0],source_feature_counts=stats,
        transcript_orphan_parent_counts={'Ensembl':0,'RefSeq':0},coordinate_convention='1-based inclusive',
        all_source_sequences_included=True,provider_records_merged=False,source_resource_modified=False,
        nested_RNA_parent_records=nested)
    db.close();save(output/'initial_plan.private.json',original);save(output/'plan.private.json',plan)
    save(output/'repair.private.json',dict(reason='Resolve assigned-molecule Y alias; independently recount every source feature',
        original_script_sha256=original['sources']['scripts/build_gene_inventory.py']))
    finish(output,plan,summary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--resume-incomplete', action='store_true', help='Validate and finalize a fully parsed incomplete build')
    args=parser.parse_args(); config=read_json(ROOT/'config/gene_inventory.json')
    resources=ROOT/config['resource_directory']; output=args.output or ROOT/config['output_directory']
    assets=read_json(resources/'assets.json')
    paths=[Path(__file__), ROOT/'scripts/extension_support.py',ROOT/'scripts/notebook_support.py',ROOT/'config/gene_inventory.json',
           resources/'ensembl.gtf.gz',resources/'refseq.gff.gz',resources/'assembly_report.txt',resources/'assets.json']
    for name in ['ensembl.gtf.gz','refseq.gff.gz','assembly_report.txt']:
        if file_digest(resources/name)!=assets['files'][name]['sha256']:
            raise ValueError('Source annotation differs from verified asset receipt')
    plan=dict(config=config,sources=sources(paths),sqlite_version=sqlite3.sqlite_version)
    if args.resume_incomplete and output.exists() and not (output/'receipt.private.json').exists():
        resume_inventory(output,plan,resources);return
    if not start(output,plan):return
    db=sqlite3.connect(output/'dog_ros.sqlite')
    db.executescript('''
    PRAGMA foreign_keys=ON;
    CREATE TABLE metadata(key TEXT PRIMARY KEY, value TEXT NOT NULL);
    CREATE TABLE sequence_aliases(provider TEXT, source_sequence TEXT, assigned_molecule TEXT, genbank_accession TEXT,
        refseq_accession TEXT, relationship TEXT, length INTEGER, role TEXT, shared_sequence_key TEXT,
        PRIMARY KEY(provider,source_sequence));
    CREATE TABLE genes(record_key TEXT PRIMARY KEY,provider TEXT,assembly TEXT,release INTEGER,feature_id TEXT,stable_gene_id TEXT,
        gene_version TEXT,symbol TEXT,name TEXT,biotype TEXT,gene_class TEXT,sequence TEXT,start INTEGER,end INTEGER,strand TEXT,
        entrez_ids_json TEXT,attributes_json TEXT);
    CREATE TABLE transcripts(record_key TEXT PRIMARY KEY,provider TEXT,assembly TEXT,release INTEGER,feature_id TEXT,
        transcript_id TEXT,transcript_version TEXT,parent_gene_keys_json TEXT,biotype TEXT,sequence TEXT,start INTEGER,end INTEGER,
        strand TEXT,attributes_json TEXT);
    CREATE TABLE features(row_id INTEGER PRIMARY KEY,provider TEXT,sequence TEXT,type TEXT,start INTEGER,end INTEGER,
        strand TEXT,phase TEXT,parent_keys_json TEXT,attributes_json TEXT);
    ''')
    aliases=[]
    for line in (resources/'assembly_report.txt').read_text().splitlines():
        if not line or line.startswith('#'):continue
        r=line.split('\t'); names={'Ensembl':r[2] if r[1]=='assembled-molecule' and r[2]!='na' else r[0], 'RefSeq':r[6]}
        for provider,name in names.items():
            if name=='na':continue
            shared=r[6] if r[5]=='=' and r[6]!='na' else ''
            aliases.append((provider,name,r[2],r[4],r[6],r[5],int(r[8]),r[1],shared))
    db.executemany('INSERT INTO sequence_aliases VALUES (?,?,?,?,?,?,?,?,?)',aliases)
    stats={}; source_headers={}; source_gene_count={}; orphan_counts={}
    for provider,filename,assembly,release in [('Ensembl','ensembl.gtf.gz',config['ensembl_accession'],config['ensembl_release']),
        ('RefSeq','refseq.gff.gz',config['refseq_accession'],config['refseq_release'])]:
        headers=[]; counts=Counter(); batch=[]; genes=[]; transcripts=[]; seen_gene_ids=set(); seen_tx_ids=set()
        with gzip.open(resources/filename,'rt',encoding='utf-8') as handle:
            for line in handle:
                if line.startswith('#'):
                    if line.startswith('#!'):headers.append(line.strip())
                    continue
                f=line.rstrip('\n').split('\t')
                if len(f)!=9:raise ValueError('Invalid annotation source line')
                seq,kind,start_bp,end_bp,strand=f[0],f[2],int(f[3]),int(f[4]),f[6]; a=attrs(f[8],provider)
                if start_bp<1 or end_bp<start_bp:raise ValueError('Invalid inclusive annotation interval')
                counts[kind]+=1
                if kind in {'gene','pseudogene'}:
                    feature_id=a['gene_id'] if provider=='Ensembl' else a['ID']
                    entrez=[x.split(':',1)[1] for x in a.get('Dbxref','').split(',') if x.startswith('GeneID:')]
                    stable=feature_id if provider=='Ensembl' else entrez[0] if len(entrez)==1 else ''
                    biotype=a.get('gene_biotype','pseudogene' if kind=='pseudogene' else '')
                    key=provider+':'+feature_id
                    if key in seen_gene_ids:raise ValueError('Duplicate provider gene feature ID')
                    seen_gene_ids.add(key)
                    genes.append((key,provider,assembly,release,feature_id,stable,a.get('gene_version',''),
                        a.get('gene_name',a.get('gene',a.get('Name',''))),a.get('description',''),biotype,gene_class(biotype),seq,
                        start_bp,end_bp,strand,json.dumps(entrez),json.dumps(a,sort_keys=True)))
                elif (provider=='Ensembl' and kind=='transcript') or (provider=='RefSeq' and a.get('Parent') and kind not in {'exon','CDS','region','cDNA_match','match','match_part'} and a.get('ID')):
                    feature_id=a['transcript_id'] if provider=='Ensembl' else a['ID']
                    tid=a.get('transcript_id',a.get('Name',feature_id))
                    parents=[provider+':'+a['gene_id']] if provider=='Ensembl' else [provider+':'+x for x in a['Parent'].split(',')]
                    key=provider+':'+feature_id
                    if key in seen_tx_ids:raise ValueError('Duplicate transcript feature ID')
                    seen_tx_ids.add(key)
                    transcripts.append((key,provider,assembly,release,feature_id,tid,a.get('transcript_version',''),json.dumps(parents),
                        a.get('transcript_biotype',kind),seq,start_bp,end_bp,strand,json.dumps(a,sort_keys=True)))
                elif kind!='region':
                    parents=[provider+':'+a['transcript_id']] if provider=='Ensembl' and a.get('transcript_id') else [provider+':'+x for x in a.get('Parent','').split(',') if x]
                    batch.append((provider,seq,kind,start_bp,end_bp,strand,f[7],json.dumps(parents),json.dumps(a,sort_keys=True)))
                if len(batch)>=5000:
                    db.executemany('INSERT INTO features(provider,sequence,type,start,end,strand,phase,parent_keys_json,attributes_json) VALUES (?,?,?,?,?,?,?,?,?)',batch);batch=[]
        db.executemany('INSERT INTO genes VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)',genes)
        db.executemany('INSERT INTO transcripts VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)',transcripts)
        db.executemany('INSERT INTO features(provider,sequence,type,start,end,strand,phase,parent_keys_json,attributes_json) VALUES (?,?,?,?,?,?,?,?,?)',batch)
        known=seen_gene_ids|seen_tx_ids
        orphan_counts[provider]=sum(parent not in known for r in transcripts for parent in json.loads(r[7]))
        if orphan_counts[provider]:raise ValueError('Transcript parent gene was not found')
        source_headers[provider]=headers; stats[provider]=dict(counts)
        source_gene_count[provider]=len(genes)
        db.commit()
        print(provider,'genes',len(genes),'transcripts',len(transcripts),flush=True)
    db.executescript('CREATE INDEX gene_position ON genes(provider,sequence,start,end); CREATE INDEX transcript_position ON transcripts(provider,sequence,start,end); CREATE INDEX feature_position ON features(provider,sequence,start,end);')
    db.execute('INSERT INTO metadata VALUES (?,?)',('configuration',json.dumps(config)))
    db.execute('INSERT INTO metadata VALUES (?,?)',('source_headers',json.dumps(source_headers)))
    if db.execute('PRAGMA integrity_check').fetchone()[0]!='ok':raise ValueError('SQLite integrity check failed')
    gene_counts={provider:count for provider,count in db.execute('SELECT provider,count(*) FROM genes GROUP BY provider')}
    if gene_counts!=source_gene_count:raise ValueError('Inventory row counts disagree with source parsing')
    alias_gaps=db.execute('SELECT DISTINCT g.provider,g.sequence FROM genes g LEFT JOIN sequence_aliases a ON g.provider=a.provider AND g.sequence=a.source_sequence WHERE a.source_sequence IS NULL').fetchall()
    # Ensembl GTF uses GenBank accession names for some scaffolds. Resolve exact
    # accession aliases, never approximate coordinates or symbols.
    for provider,seq in alias_gaps:
        r=db.execute('SELECT * FROM sequence_aliases WHERE provider=? AND genbank_accession=?',(provider,seq)).fetchone()
        if r:
            db.execute('INSERT INTO sequence_aliases VALUES (?,?,?,?,?,?,?,?,?)',(provider,seq,*r[2:]))
    db.commit()
    alias_gaps=db.execute('SELECT DISTINCT g.provider,g.sequence FROM genes g LEFT JOIN sequence_aliases a ON g.provider=a.provider AND g.sequence=a.source_sequence WHERE a.source_sequence IS NULL').fetchall()
    if alias_gaps:raise ValueError(f'Unresolved source sequence aliases: {alias_gaps[:5]}')
    nested=resolve_transcript_genes(db);db.commit()
    export_table(db,output/'genes.tsv.gz','SELECT * FROM genes ORDER BY provider,sequence,start,record_key')
    export_table(db,output/'transcripts.tsv.gz','SELECT * FROM transcripts ORDER BY provider,sequence,start,record_key')
    summary=dict(status='complete_gene_inventory',gene_counts=gene_counts,
        transcript_counts=dict(db.execute('SELECT provider,count(*) FROM transcripts GROUP BY provider')),
        gene_classes=[dict(provider=p,gene_class=c,count=n) for p,c,n in db.execute('SELECT provider,gene_class,count(*) FROM genes GROUP BY provider,gene_class')],
        feature_rows=db.execute('SELECT count(*) FROM features').fetchone()[0],source_feature_counts=stats,
        transcript_orphan_parent_counts=orphan_counts,coordinate_convention='1-based inclusive',all_source_sequences_included=True,
        provider_records_merged=False,source_resource_modified=False,nested_RNA_parent_records=nested)
    db.close();finish(output,plan,summary)


if __name__=='__main__':main()
