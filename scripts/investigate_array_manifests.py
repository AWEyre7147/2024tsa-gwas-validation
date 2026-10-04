"""Compare public Illumina canine manifests with local marker labels.

Compatibility counts are evidence of content overlap, not proof of export strand.
No genotype calls are altered and no pilot exclusions are reversed.
"""
from pathlib import Path
from collections import Counter
import concurrent.futures
import csv
import hashlib
import json
import shutil
import subprocess

ROOT=Path(__file__).resolve().parents[1]
DEST=ROOT/'data/reference/array-manifests'
URLS={
 'caninehd_b.csv':'https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/caninehd/caninehd_b.csv',
 'CanineHTS-24_20095584_B1.csv':'https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/caninehd/CanineHTS-24_20095584_B1.csv',
}
COMP=str.maketrans('ACGT','TGCA')

def fetch(item):
 name,url=item;p=DEST/name
 if not p.exists():
  temp=p.with_name(p.name+'.partial');curl=shutil.which('curl.exe') or shutil.which('curl')
  subprocess.run([curl,'--fail','--location','--silent','--show-error','--retry','2','--max-time','180',url,'--output',str(temp)],check=True);temp.replace(p)
 return name,p

def parse(p):
 header=[];data={};duplicates=0
 with p.open(encoding='utf-8-sig',newline='') as f:
  reader=csv.reader(f)
  for r in reader:
   if r and r[0]=='[Assay]':break
   header.append(r)
  columns=next(reader)
  for r in reader:
   if r and r[0].startswith('['):break
   if not r:continue
   d=dict(zip(columns,r));name=d.get('Name')
   if name in data:duplicates+=1
   if name:data[name]=d
 return header,columns,data,duplicates

def main():
 DEST.mkdir(parents=True,exist_ok=True)
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:paths=dict(pool.map(fetch,URLS.items()))
 source=next((ROOT/'data/raw').rglob('*.bim'));bim={r[1]:r for r in [x.split() for x in source.read_text().splitlines()]}
 excluded={}
 for b in ['GSP','LR']:
  p=ROOT/f'data/archive/first-imputation-pilot-20261003/analysis/{b}/harmonization.private.tsv'
  excluded[b]=[r['marker'] for r in csv.DictReader(p.open(),delimiter='\t') if r['status']=='palindromic_excluded']
 summaries={}
 for name,p in paths.items():
  header,columns,data,duplicates=parse(p);matched=set(bim)&set(data);counts=Counter();builds=Counter();strandvalues=Counter();pal=[]
  for s in matched:
   d=data[s];r=bim[s];snp=d.get('SNP','');pair=set(snp.strip('[]').split('/'));observed=set(r[4:6]);builds[d.get('GenomeBuild','')]+=1
   for col in ['IlmnStrand','RefStrand','SourceStrand']:strandvalues[col+':'+d.get(col,'')]+=1
   if len(pair)!=2 or not pair<=set('ACGT'):continue
   if pair in ({'A','T'},{'C','G'}):continue
   counts['nonpalindromic_compared']+=1
   plus={a.translate(COMP) for a in pair} if d.get('RefStrand')=='-' else pair
   top={a.translate(COMP) for a in pair} if d.get('IlmnStrand')=='BOT' else pair
   counts['matches_design_allele_pair']+=observed==pair
   if d.get('RefStrand') in ('+','-'):counts['matches_plus_allele_pair']+=observed==plus
   if d.get('IlmnStrand') in ('TOP','BOT'):counts['matches_top_allele_pair']+=observed==top
  for b,ids in excluded.items():
   for s in ids:
    d=data.get(s,{})
    pal.append(dict(breed=b,marker=s,in_manifest=bool(d),IlmnID=d.get('IlmnID',''),SNP=d.get('SNP',''),IlmnStrand=d.get('IlmnStrand',''),RefStrand=d.get('RefStrand',''),SourceStrand=d.get('SourceStrand',''),GenomeBuild=d.get('GenomeBuild',''),Chr=d.get('Chr',''),MapInfo=d.get('MapInfo',''),SourceSeq=d.get('SourceSeq',''),TopGenomicSeq=d.get('TopGenomicSeq',''),AlleleA_ProbeSeq=d.get('AlleleA_ProbeSeq','')))
  with (DEST/(name+'.pilot_palindromic_records.private.json')).open('w') as f:json.dump(pal,f,indent=2)
  summaries[name]=dict(url=URLS[name],sha256=hashlib.sha256(p.read_bytes()).hexdigest(),header=header,columns=columns,manifest_unique_markers=len(data),duplicate_names=duplicates,raw_markers=len(bim),raw_marker_ID_matches=len(matched),manifest_markers_absent_from_raw=len(set(data)-set(bim)),raw_markers_absent_from_manifest=len(set(bim)-set(data)),genome_build_values=dict(builds),strand_values=dict(strandvalues),nonpalindromic_pair_compatibility=dict(counts),pilot_palindromic_coverage={b:dict(excluded=len(ids),manifest_matches=sum(s in data for s in ids),with_reference_strand=sum(data.get(s,{}).get('RefStrand') in ('+','-') for s in ids),with_probe_sequence=sum(bool(data.get(s,{}).get('AlleleA_ProbeSeq')) for s in ids)) for b,ids in excluded.items()})
 report=dict(status='investigation_only',raw_bim_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),manifests=summaries,limitation='Marker overlap and allele-pair compatibility do not establish the historical export convention or exact array revision. No genotypes changed.')
 (DEST/'manifest_comparison.private.json').write_text(json.dumps(report,indent=2)+'\n')
 for name,s in summaries.items():print(name,json.dumps({k:s[k] for k in ['manifest_unique_markers','raw_marker_ID_matches','raw_markers_absent_from_manifest','genome_build_values','nonpalindromic_pair_compatibility','pilot_palindromic_coverage']},indent=2))

if __name__=='__main__':main()
