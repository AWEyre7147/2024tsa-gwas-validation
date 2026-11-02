"""Align QC SNPs to CanFam4; preserve source files and verify every retained call.

CanineHD is investigator-confirmed. TOP reporting is supported by manifest
non-palindromic allele pairs. Unresolved sites are reported rather than guessed.
Run in Windows or Linux Python; UCSC liftOver uses the existing Linux environment.
"""
from pathlib import Path
from collections import Counter
import bisect
import csv
import gzip
import hashlib
import json
import mmap
import os
import re
import struct
import subprocess
import numpy as np
from investigate_array_manifests import parse

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/02-aligned'
REF = ROOT / 'data/02-reference/genomes'
COMP = str.maketrans('ACGT', 'TGCA')


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        while b := f.read(1024*1024): h.update(b)
    return h.hexdigest()


def rows(p):
    return [line.split() for line in p.read_text().splitlines() if line.strip()]


def dump(p, obj):
    p.write_text(json.dumps(obj, indent=2)+'\n', encoding='utf-8')


class TwoBit:
    """Read UCSC two-bit sequence ranges, including N blocks, without extraction."""
    def __init__(self, path):
        self.file = path.open('rb')
        self.data = mmap.mmap(self.file.fileno(), 0, access=mmap.ACCESS_READ)
        self.endian = '>' if self.data[:4] == bytes.fromhex('1a412743') else '<'
        assert self.u32(0) == 0x1a412743 and self.u32(4) == 0
        self.index = {}; offset = 16
        for _ in range(self.u32(8)):
            size = self.data[offset]; offset += 1
            name = self.data[offset:offset+size].decode(); offset += size
            record = self.u32(offset); offset += 4
            length = self.u32(record); n = self.u32(record+4)
            starts = [self.u32(record+8+4*i) for i in range(n)]
            sizes = [self.u32(record+8+4*n+4*i) for i in range(n)]
            p = record+8+8*n; masks = self.u32(p)
            packed = p+4+8*masks+4
            self.index[name] = length, packed, starts, sizes

    def u32(self, offset):
        return struct.unpack_from(self.endian+'I', self.data, offset)[0]

    def get(self, chrom, start, end):
        if chrom not in self.index: return ''
        length, packed, starts, sizes = self.index[chrom]
        if start < 0 or end > length: return ''
        seq = ['TCAG'[(self.data[packed+i//4] >> (6-2*(i%4))) & 3] for i in range(start,end)]
        j = max(0,bisect.bisect_right(starts,start)-1)
        for k in range(j,len(starts)):
            if starts[k] >= end: break
            for i in range(max(start,starts[k]),min(end,starts[k]+sizes[k])): seq[i-start]='N'
        return ''.join(seq)


def context(genome, chrom, pos, sequence):
    m = re.fullmatch(r'([A-Z]+)\[([ACGT])/([ACGT])\]([A-Z]+)', sequence.upper())
    if not m: return None, {}
    left,a,b,right = m.groups(); scores = {}
    for strand,l,r in [('+',left,right),('-',right.translate(COMP)[::-1],left.translate(COMP)[::-1])]:
        l,r=l[-25:],r[:25]
        reference=genome.get(chrom,pos-1-len(l),pos-1)+genome.get(chrom,pos,pos+len(r))
        pairs=[(x,y) for x,y in zip(l+r,reference) if x in 'ACGT' and y in 'ACGT']
        scores[strand]=(len(pairs),sum(x!=y for x,y in pairs))
    passing=[s for s,(n,e) in scores.items() if n>=45 and e<=2]
    return (passing[0] if len(passing)==1 else None),scores


def wslpath(p):
    p=p.resolve().as_posix(); return '/mnt/'+p[0].lower()+p[2:]


def main():
    if (OUT/'manifest.private.json').exists(): raise FileExistsError('Completed alignment already exists')
    OUT.mkdir(exist_ok=True); reports=OUT/'reports'; reports.mkdir(exist_ok=True)
    executed_source=Path(__file__).read_bytes()
    (reports/'reference_alignment.executed.py').write_bytes(executed_source)
    qcpath=ROOT/'data/02-qc/manifest.private.json'; qc=json.loads(qcpath.read_text())
    assert qc['status']=='complete'
    original={p:h for ds in qc['datasets'] for p,h in ds['files'].items()}
    assert all(sha(ROOT/p)==h for p,h in original.items())
    datasets=[]; union={}
    for ds in qc['datasets']:
        prefix=next(ROOT/p for p in ds['files'] if p.endswith('.fam')).with_suffix('')
        variants=rows(prefix.with_suffix('.bim')); datasets.append((ds,prefix,variants))
        for r in variants:
            assert r[1] not in union or union[r[1]]==r
            union[r[1]]=r
    bed=reports/'qc_sites.canfam3.bed'
    bed.write_text(''.join(f'chr{r[0]}\t{int(r[3])-1}\t{r[3]}\t{s}\t0\t+\n' for s,r in union.items()))
    assets=ROOT/'data/02-reference/imputation-pilot'
    mapped=reports/'qc_sites.canfam4.bed'; unmapped=reports/'unmapped.canfam4.bed'
    if os.name == 'nt':
        command=['wsl.exe','-d','Ubuntu','--','/home/aweyr/anaconda3/envs/tsa2024-imputation/bin/liftOver',wslpath(bed),wslpath(assets/'canFam3ToCanFam4.over.chain.gz'),wslpath(mapped),wslpath(unmapped)]
    else:
        command=[str(Path.home()/'anaconda3/envs/tsa2024-imputation/bin/liftOver'),str(bed),str(assets/'canFam3ToCanFam4.over.chain.gz'),str(mapped),str(unmapped)]
    result=subprocess.run(command,capture_output=True,text=True)
    dump(reports/'liftover_command.private.json',dict(command=command,returncode=result.returncode,stdout=result.stdout,stderr=result.stderr))
    result.check_returncode()
    hits={r[3]:r for r in rows(mapped)}; destinations=Counter((r[0],r[1]) for r in hits.values())
    cf3,cf4=TwoBit(REF/'canFam3.2bit'),TwoBit(REF/'canFam4.2bit')
    manifestpath=ROOT/'data/02-reference/array-manifests/caninehd_b.csv'
    _,_,manifest,_=parse(manifestpath)
    decisions={}; convention=Counter()
    for s,r in union.items():
        pair=set(r[4:6]); pal=pair in ({'A','T'},{'C','G'})
        h=hits.get(s); d=manifest.get(s)
        decision=dict(marker=s,source_chr=r[0],source_pos=int(r[3]),palindromic=pal,in_caninehd=bool(d),status='unmapped')
        decisions[s]=decision
        if d and not pal:
            design=set(d['SNP'].strip('[]').split('/'))
            top={a.translate(COMP) for a in design} if d['IlmnStrand']=='BOT' else design
            convention['compared']+=1; convention['top_matches']+=pair==top
            if pair!=top: decision['status']='manifest_allele_disagreement';continue
        if not h: continue
        chrom,pos=h[0],int(h[1])+1
        decision.update(target_chr=chrom,target_pos=pos,chain_strand=h[5])
        if chrom not in [f'chr{i}' for i in range(1,39)]: decision['status']='nonautosomal_destination';continue
        if destinations[(chrom,h[1])]!=1: decision['status']='duplicate_destination';continue
        ref3=cf3.get('chr'+r[0],int(r[3])-1,int(r[3])); ref4=cf4.get(chrom,pos-1,pos)
        decision.update(source_reference=ref3,target_reference=ref4)
        if ref3 not in 'ACGT' or len(ref3)!=1 or ref4 not in 'ACGT' or len(ref4)!=1:
            decision['status']='unknown_reference_base';continue
        source_strand=None
        if d:
            source_strand,score=context(cf3,'chr'+r[0],int(r[3]),d.get('TopGenomicSeq',''))
            decision['source_context_scores']=score
        if pal and not source_strand:
            decision['status']='palindromic_context_unresolved' if d else 'palindromic_manifest_absent';continue
        if not source_strand:
            if ref3 in pair: source_strand='+'
            elif ref3 in {a.translate(COMP) for a in pair}: source_strand='-'
            else: decision['status']='source_reference_mismatch';continue
            decision['orientation_method']='nonpalindromic_reference_alleles'
        else: decision['orientation_method']='caninehd_top_sequence'
        plus3=[a.translate(COMP) for a in r[4:6]] if source_strand=='-' else r[4:6]
        if ref3 not in plus3: decision['status']='source_reference_mismatch';continue
        strand4=source_strand if h[5]=='+' else ('-' if source_strand=='+' else '+')
        a1,a2=[a.translate(COMP) for a in r[4:6]] if strand4=='-' else r[4:6]
        if ref4 not in (a1,a2): decision['status']='target_reference_mismatch';continue
        if pal:
            observed,score=context(cf4,chrom,pos,d['TopGenomicSeq'])
            decision['target_context_scores']=score
            if observed!=strand4: decision['status']='palindromic_target_context_unresolved';continue
        decision.update(status='aligned',source_to_plus=source_strand,source_to_target=strand4,REF=ref4,ALT=a2 if a1==ref4 else a1,swap_codes=(a1==ref4))
    assert convention['compared']>100000 and convention['compared']==convention['top_matches'], 'TOP convention evidence disagrees'
    dump(reports/'marker_alignment.private.json',list(decisions.values()))
    outputdatasets=[]
    lookup=np.array([sum(({0:3,1:1,2:2,3:0}[(b>>(2*j))&3])<<(2*j) for j in range(4)) for b in range(256)],dtype=np.uint8)
    for ds,prefix,variants in datasets:
        folder=OUT/ds['scenario']/ds['breed'];folder.mkdir(parents=True,exist_ok=True)
        dest=folder/prefix.name; fam=prefix.with_suffix('.fam').read_bytes(); n=len(rows(prefix.with_suffix('.fam')));stride=(n+3)//4
        raw=prefix.with_suffix('.bed').read_bytes();assert raw[:3]==b'\x6c\x1b\x01' and len(raw)==3+len(variants)*stride
        retained=sorted([(i,r,decisions[r[1]]) for i,r in enumerate(variants) if decisions[r[1]]['status']=='aligned'],key=lambda t:(int(t[2]['target_chr'][3:]),t[2]['target_pos']))
        dest.with_suffix('.fam').write_bytes(fam)
        with dest.with_suffix('.bed').open('wb') as f:
            f.write(raw[:3])
            for i,r,d in retained:
                block=np.frombuffer(raw,dtype=np.uint8,count=stride,offset=3+i*stride)
                block=lookup[block] if d['swap_codes'] else block.copy()
                if n%4: block[-1]&=(1<<(2*(n%4)))-1
                f.write(block.tobytes())
        dest.with_suffix('.bim').write_text(''.join(f"{d['target_chr'][3:]}\t{r[1]}\t0\t{d['target_pos']}\t{d['ALT']}\t{d['REF']}\n" for i,r,d in retained))
        # Export VCF directly from source calls, independently of output BED recoding.
        vcf=dest.with_suffix('.vcf.gz');samples=[r[1] for r in rows(prefix.with_suffix('.fam'))]
        with gzip.open(vcf,'wt') as f:
            f.write('##fileformat=VCFv4.2\n##reference=UU_Cfam_GSD_1.0\n')
            for chrom in sorted({d['target_chr'] for _,_,d in retained},key=lambda c:int(c[3:])): f.write(f'##contig=<ID={chrom},length={cf4.index[chrom][0]}>\n')
            f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Unphased genotype">\n')
            f.write('\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',*samples])+'\n')
            for i,r,d in retained:
                aa=[a.translate(COMP) for a in r[4:6]] if d['source_to_target']=='-' else r[4:6]
                labels=['0' if a==d['REF'] else '1' for a in aa]
                coding={0:labels[0]+'/'+labels[0],1:'./.',2:'/'.join(sorted(labels)),3:labels[1]+'/'+labels[1]}
                calls=[coding[(raw[3+i*stride+j//4]>>(2*(j%4)))&3] for j in range(n)]
                f.write('\t'.join([d['target_chr'],str(d['target_pos']),r[1],d['REF'],d['ALT'],'.','PASS','.','GT',*calls])+'\n')
        output=dest.with_suffix('.bed').read_bytes();js=np.arange(n)
        with gzip.open(vcf,'rt') as f:
            records=(line.rstrip().split('\t') for line in f if not line.startswith('#'))
            for k,((i,r,d),v) in enumerate(zip(retained,records,strict=True)):
                old=(np.frombuffer(raw,dtype=np.uint8,count=stride,offset=3+i*stride)[js//4]>>(2*(js%4)))&3
                actual=(np.frombuffer(output,dtype=np.uint8,count=stride,offset=3+k*stride)[js//4]>>(2*(js%4)))&3
                expected=np.array([3,1,2,0],dtype=np.uint8)[old] if d['swap_codes'] else old
                assert np.array_equal(actual,expected)
                genotypes=np.array(['1/1','./.','0/1','0/0'])[actual].tolist()
                assert genotypes==v[9:] and v[3:5]==[d['REF'],d['ALT']]
        counts=Counter(decisions[r[1]]['status'] for r in variants)
        record=dict(scenario=ds['scenario'],breed=ds['breed'],dogs=n,input_variants=len(variants),aligned_variants=len(retained),palindromic_aligned=sum(d['palindromic'] for _,_,d in retained),statuses=dict(counts),all_calls_verified=True,files={p.relative_to(ROOT).as_posix():sha(p) for p in [dest.with_suffix(e) for e in ('.bed','.bim','.fam','.vcf.gz')]})
        outputdatasets.append(record);print(json.dumps({k:v for k,v in record.items() if k!='files'}),flush=True)
    assert all(sha(ROOT/p)==h for p,h in original.items())
    manifest=dict(status='complete',source_assembly='CanFam3.1',target_assembly='CanFam4 / UU_Cfam_GSD_1.0',platform='CanineHD, investigator confirmed',allele_convention='TOP, supported by exact non-palindromic manifest allele comparisons; palindromic interpretation conditional on this convention',convention_check=dict(convention),source_qc_unchanged=True,source_manifest_sha256=sha(qcpath),script_sha256=hashlib.sha256(executed_source).hexdigest(),resources={str(p.relative_to(ROOT)):sha(p) for p in [REF/'canFam3.2bit',REF/'canFam4.2bit',manifestpath,assets/'canFam3ToCanFam4.over.chain.gz']},datasets=outputdatasets,limitations=['Coordinate-anchored sequence context, not a genome-wide probe uniqueness test.','Reference-panel variant compatibility is a subsequent imputation step.','ROS coordinates remain a separate annotation crosswalk; this is not a ROS-aligned genotype dataset.'])
    dump(OUT/'manifest.private.json',manifest)


if __name__=='__main__': main()
