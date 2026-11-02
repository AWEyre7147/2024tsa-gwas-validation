"""Resumable autosomal Beagle workflow from verified CanFam4 aligned genotypes.

Linux/WSL: python scripts/imputation_workflow.py --prepare-only
Linux/WSL: python scripts/imputation_workflow.py --execute
Preparation splits all chromosomes and freezes the plan without downloading the
whole reference panel or running Beagle. --execute performs the complete plan.
"""
from pathlib import Path
from collections import Counter
import argparse
import csv
import gzip
import hashlib
import json
import os
import random
import re
import subprocess
import sys
import tarfile
import time
import numpy as np
from reference_alignment import TwoBit
from cache_dog10k_panel import ensure_panel

ROOT=Path(__file__).resolve().parents[1]
ASSETS=ROOT/'data/02-reference/imputation-pilot'
PANEL='https://kiddlabshare.med.umich.edu/dog10K/phased-imputation-panel/AutoAndXPAR.Dog10K.phased.bcf'


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        while chunk:=f.read(1048576):h.update(chunk)
    return h.hexdigest()


def dump(p,o):
    temp=p.with_name(p.name+'.partial');temp.write_text(json.dumps(o,indent=2)+'\n');temp.replace(p)


def vcf(path):
    with gzip.open(path,'rt') as f:
        for line in f:
            if not line.startswith('#'):yield line.rstrip().split('\t')


def count_gt(gt):
    return np.nan if '.' in gt else sum(int(a) for a in gt.replace('|','/').split('/'))


def correlation(x,y):
    ok=np.isfinite(x)&np.isfinite(y);x=x[ok];y=y[ok]
    return None if len(x)<3 or np.std(x)==0 or np.std(y)==0 else float(np.corrcoef(x,y)[0,1]**2)


def same_location(a,b):
    """Accept Windows-mounted path aliases only when they identify the same file."""
    if a==b:return True
    if not all(isinstance(p,str) and re.match(r'^/mnt/[a-zA-Z]/',p) for p in (a,b)):return False
    try:return os.path.samefile(a,b)
    except OSError:return False


def plans_equal(saved,current):
    if set(saved)!=set(current):return False
    for key in saved:
        if key in ('tools','reference_source'):
            if not same_location(saved[key],current[key]):return False
        elif key=='resource_hashes':
            old,new=saved[key],current[key]
            if len(old)!=len(new):return False
            remaining=dict(new)
            for path,digest in old.items():
                matches=[p for p,h in remaining.items() if h==digest and same_location(path,p)]
                if len(matches)!=1:return False
                del remaining[matches[0]]
        elif saved[key]!=current[key]:return False
    return True


class Workflow:
    def __init__(self,args):
        self.args=args;self.out=args.output.resolve();self.out.mkdir(parents=True,exist_ok=True)
        self.logs=self.out/'logs';self.logs.mkdir(exist_ok=True)
        self.tools=args.tools;self.config=json.loads(args.config.read_text())
        self.env=dict(os.environ);self.env['PATH']=str(self.tools)+':'+self.env.get('PATH','')
        self.inputs=json.loads((ROOT/'data/02-aligned/manifest.private.json').read_text())
        assert self.inputs['status']=='complete'
        self.sources={p:h for d in self.inputs['datasets'] for p,h in d['files'].items()}
        self.verify_sources()
        self.prefixes={d['breed']:next(ROOT/p for p in d['files'] if p.endswith('.vcf.gz')) for d in self.inputs['datasets'] if d['scenario']=='not-behelim-vs-behelim'}
        self.scheme_a={d['breed']:[line.split()[1] for line in (next(ROOT/p for p in d['files'] if p.endswith('.fam'))).read_text().splitlines()] for d in self.inputs['datasets'] if d['scenario']=='success-vs-behelim'}
        self.regions=args.regions or ['chr'+str(c) for c in self.config['chromosomes']]
        assert len(set(self.regions))==len(self.regions)
        assert all(re.fullmatch(r'chr(?:[1-9]|[12][0-9]|3[0-8])(?::[1-9][0-9]*-[1-9][0-9]*)?',r) for r in self.regions)
        self.keys={r:r.replace(':','_').replace('-','_') for r in self.regions}
        assets=json.loads((ASSETS/'assets.json').read_text())
        pinned={str(ASSETS/name):assets[name]['sha256'] for name in ['beagle.27Feb25.75f.jar','campbell_sex_average_canFam4.tar.gz','samples.txt','AutoAndXPAR.Dog10K.phased.bcf.csi','panel.md5sums.txt']}
        pinned[str(ROOT/'data/02-reference/genomes/canFam4.2bit')]=self.inputs['resources']['data\\02-reference\\genomes\\canFam4.2bit'] if 'data\\02-reference\\genomes\\canFam4.2bit' in self.inputs['resources'] else self.inputs['resources']['data/02-reference/genomes/canFam4.2bit']
        assert all(sha(Path(p))==h for p,h in pinned.items()),'Reference or software asset changed'
        versions={name:subprocess.check_output([str(self.tools/name),'--version'],env=self.env,text=True).splitlines()[0] for name in ['bcftools','java']}
        versions.update(python=sys.version.split()[0],numpy=np.__version__)
        self.plan=dict(config=self.config,regions=self.regions,source_manifest_sha256=sha(ROOT/'data/02-aligned/manifest.private.json'),reference_source=args.reference_source or PANEL,tools=str(self.tools),versions=versions,resource_hashes=pinned,script_sha256=sha(Path(__file__)))
        self.plan['dependency_script_hashes']={name:sha(ROOT/'scripts'/name) for name in ['reference_alignment.py','investigate_array_manifests.py','cache_dog10k_panel.py']}
        if args.reference_source and Path(args.reference_source).is_file():self.plan['local_reference_sha256']=sha(Path(args.reference_source))
        planpath=self.out/'plan.private.json'
        if planpath.exists():assert plans_equal(json.loads(planpath.read_text()),self.plan),'Plan changed: use a new output directory'
        else:dump(planpath,self.plan)
        (self.out/'implementation.py').write_bytes(Path(__file__).read_bytes())
        self.genome=TwoBit(ROOT/'data/02-reference/genomes/canFam4.2bit')
        self.statepath=self.out/'manifest.private.json'
        self.state=json.loads(self.statepath.read_text()) if self.statepath.exists() else dict(status='preparing',assembly='CanFam4 / UU_Cfam_GSD_1.0',jobs={},association_ready=False,reference_independence_review='pending; genotype similarity screen is not proof of independence')

    def verify_sources(self):
        assert all(sha(ROOT/p)==h for p,h in self.sources.items()),'Aligned input changed'

    def save(self):dump(self.statepath,self.state)

    def run(self,label,args,outputs=()):
        args=[str(x) for x in args];record=self.logs/(label+'.json')
        if record.exists():
            old=json.loads(record.read_text())
            if old['returncode']==0 and old['argv']==args and all(p.exists() and sha(p)==old['outputs'].get(str(p)) for p in outputs):return
        start=time.time()
        with (self.logs/(label+'.txt')).open('w') as log:r=subprocess.run(args,stdout=log,stderr=subprocess.STDOUT,env=self.env)
        dump(record,dict(argv=args,returncode=r.returncode,elapsed_seconds=time.time()-start,outputs={str(p):sha(p) for p in outputs if p.exists()}))
        r.check_returncode()

    def capture(self,label,args,path):
        argv=[str(x) for x in args];record=self.logs/(label+'.json')
        if record.exists():
            old=json.loads(record.read_text())
            if old['returncode']==0 and old['argv']==argv and path.exists() and sha(path)==old['outputs'].get(str(path)):return
        with path.open('w') as out,(self.logs/(label+'.txt')).open('w') as err:r=subprocess.run(argv,stdout=out,stderr=err,env=self.env)
        dump(record,dict(argv=argv,returncode=r.returncode,outputs={str(path):sha(path)}));r.check_returncode()

    def compress(self,label,source,dest):
        self.run(label,[self.tools/'bcftools','view','-Oz','-o',dest,source],[dest])
        self.run(label+'_index',[self.tools/'bcftools','index','-f',dest],[Path(str(dest)+'.csi')])

    def prepare(self):
        maps=self.out/'maps';maps.mkdir(exist_ok=True)
        with tarfile.open(ASSETS/'campbell_sex_average_canFam4.tar.gz') as archive:
            for chrom in {r.split(':')[0] for r in self.regions}:
                member=next(m for m in archive.getmembers() if m.isfile() and m.name.endswith('_'+chrom+'_map.txt'))
                rows=[r.split() for r in archive.extractfile(member).read().decode().splitlines()[1:] if r.strip()]
                assert len(rows)>1 and all(int(b[0])>int(a[0]) and float(b[2])>=float(a[2]) for a,b in zip(rows,rows[1:]))
                (maps/(chrom+'.map')).write_text(''.join(f'{chrom}\t.\t{r[2]}\t{r[0]}\n' for r in rows))
        for breed,source in self.prefixes.items():
            headers=[];handles={};counts=Counter()
            for region in self.regions:
                folder=self.out/'targets'/breed/self.keys[region];folder.mkdir(parents=True,exist_ok=True)
                handles[region]=(folder/'aligned.vcf').open('w')
            with gzip.open(source,'rt') as f:
                for line in f:
                    if line.startswith('#'):headers.append(line);continue
                    r=line.split('\t',3);chrom,pos=r[0],int(r[1])
                    for region,out in handles.items():
                        c=region.split(':')[0]
                        if c!=chrom:continue
                        if ':' in region:
                            start,end=map(int,region.split(':')[1].split('-'))
                            if not start<=pos<=end:continue
                        if counts[region]==0:out.writelines(headers)
                        out.write(line);counts[region]+=1
            for region,out in handles.items():
                out.close();assert counts[region]>20,'Insufficient target markers'
                folder=self.out/'targets'/breed/self.keys[region]
                self.compress(breed+'_'+self.keys[region]+'_target',folder/'aligned.vcf',folder/'aligned.vcf.gz')
                (folder/'aligned.vcf').unlink()
                self.state['jobs'].setdefault(region,dict(status='prepared',breeds={}))['breeds'].setdefault(breed,{})['aligned_markers']=counts[region]
        self.state['status']='prepared';self.save()
        print('Prepared',len(self.regions),'regions/chromosomes for both breeds.',flush=True)

    def reference(self,region):
        key=self.keys[region];folder=self.out/'reference'/key;folder.mkdir(parents=True,exist_ok=True)
        source=self.effective_reference
        raw=folder/'eligible.bcf'
        self.run(key+'_reference',[self.tools/'bcftools','view','-r',region,'-m2','-M2','-v','snps','-p','-g','^miss','-Ob','-o',raw,source],[raw])
        self.run(key+'_rawindex',[self.tools/'bcftools','index','-f',raw],[Path(str(raw)+'.csi')])
        sites=folder/'sites.tsv';self.capture(key+'_sites',[self.tools/'bcftools','query','-f','%CHROM\t%POS\t%REF\t%ALT\n',raw],sites)
        ploidy=folder/'ploidy.tsv';self.capture(key+'_ploidy',[self.tools/'bcftools','+check-ploidy',raw],ploidy)
        ploidies=[r.split() for r in ploidy.read_text().splitlines() if not r.startswith('#')]
        assert ploidies and all(r[-1]=='2' for r in ploidies),'Reference contains non-diploid genotypes'
        samples=folder/'samples.txt';self.capture(key+'_samples',[self.tools/'bcftools','query','-l',raw],samples)
        ids=samples.read_text().splitlines();assert len(ids)==1929 and len(set(ids))==1929
        assert {r[0] for r in ploidies}==set(ids)
        rows=[r.split() for r in sites.read_text().splitlines()];counts=Counter((r[0],r[1]) for r in rows)
        good={};excluded=Counter()
        for chrom,p,ref,alt in rows:
            pos=int(p)
            if counts[(chrom,p)]!=1:excluded['duplicate_position']+=1
            elif ref not in 'ACGT' or alt not in 'ACGT' or len(ref)!=1 or len(alt)!=1:excluded['noncanonical_alleles']+=1
            elif self.genome.get(chrom,pos-1,pos)!=ref:excluded['reference_mismatch']+=1
            else:good[pos]=(ref,alt)
        assert len(good)>100
        positions=folder/'verified_positions.bed';positions.write_text(''.join(f'{region.split(":")[0]}\t{p-1}\t{p}\n' for p in sorted(good)))
        selected=folder/'verified.bcf'
        self.run(key+'_verified',[self.tools/'bcftools','view','-T',positions,'-Ob','-o',selected,raw],[selected])
        ref=folder/'reference.vcf.gz'
        self.run(key+'_GT_only',[self.tools/'bcftools','annotate','-x','^FORMAT/GT','-Oz','-o',ref,selected],[ref])
        self.run(key+'_refindex',[self.tools/'bcftools','index','-f',ref],[Path(str(ref)+'.csi')])
        self.state['jobs'][region]['reference']=dict(samples=len(ids),variants=len(good),exclusions=dict(excluded),ref_verified=True,phased_nonmissing_diploid=True)
        self.save();return ref,good,ids

    def targets(self,breed,region,good,refids):
        key=self.keys[region];folder=self.out/'targets'/breed/key
        source=folder/'aligned.vcf.gz';records=[];headers=[];reasons=Counter()
        with gzip.open(source,'rt') as f:
            for line in f:
                if line.startswith('#'):
                    headers.append(line)
                    if line.startswith('#CHROM'):samples=line.rstrip().split('\t')[9:]
                    continue
                r=line.rstrip().split('\t');p=int(r[1])
                if p not in good:reasons['absent_or_ineligible_panel_site']+=1
                elif tuple(r[3:5])!=good[p]:reasons['panel_allele_disagreement']+=1
                else:records.append(r);reasons['compatible']+=1
        assert len(records)>20 and not set(samples)&set(refids)
        full=folder/'compatible.vcf';full.write_text(''.join(headers)+''.join('\t'.join(r)+'\n' for r in records))
        self.compress(breed+'_'+key+'_compatible',full,folder/'compatible.vcf.gz');full.unlink()
        (folder/'compatible_positions.bed').write_text(''.join(f'{r[0]}\t{int(r[1])-1}\t{r[1]}\n' for r in records))
        validation=[v for v in self.config['validation_regions'] if v['chromosome']==region.split(':')[0]]
        eligible=[r for r in records if any(v['start']<=int(r[1])<=v['end'] for v in validation)]
        truth={}
        if eligible:
            assert len(eligible)>=100,'Too few validation markers for the predefined core'
            rng=random.Random(self.config['seed']+int(region.split(':')[0][3:]))
            selected=rng.sample(eligible,max(10,round(len(eligible)*self.config['mask_fraction'])))
            truth={r[1]:dict(marker=r[2],ref=r[3],alt=r[4],gt=r[9:]) for r in selected}
            masked=folder/'masked.vcf';masked.write_text(''.join(headers)+''.join('\t'.join(r)+'\n' for r in records if r[1] not in truth))
            self.compress(breed+'_'+key+'_masked',masked,folder/'masked.vcf.gz');masked.unlink()
            assert not set(truth)&{r[1] for r in vcf(folder/'masked.vcf.gz')}
            dump(folder/'truth.private.json',truth)
        self.state['jobs'][region]['breeds'][breed].update(panel_compatibility=dict(reasons),dogs=len(samples),masked_markers=len(truth));self.save()
        return records,samples,truth

    def beagle(self,breed,region,mode,reference):
        key=self.keys[region];folder=self.out/'targets'/breed/key
        dest=(self.out/'not-behelim-vs-behelim'/breed/key if mode=='full' else self.out/'validation'/breed/key);dest.mkdir(parents=True,exist_ok=True)
        source=folder/('compatible.vcf.gz' if mode=='full' else 'masked.vcf.gz')
        prefix=dest/('imputed' if mode=='full' else 'masked');output=Path(str(prefix)+'.vcf.gz')
        args=[self.tools/'java',f'-Xmx{self.config["java_memory_gb"]}g','-jar',ASSETS/'beagle.27Feb25.75f.jar',f'gt={source}',f'ref={reference}',f'map={self.out/"maps"/(region.split(":")[0]+".map")}',f'out={prefix}','impute=true',f'gp={str(self.config["genotype_probabilities"]).lower()}',f'nthreads={self.config["threads"]}',f'seed={self.config["seed"]}',f'ne={self.config["initial_ne"]}',f'em={str(self.config["estimate_parameters"]).lower()}']
        self.run(breed+'_'+key+'_'+mode+'_beagle',args,[output,Path(str(prefix)+'.log')])
        self.run(breed+'_'+key+'_'+mode+'_index',[self.tools/'bcftools','index','-f',output],[Path(str(output)+'.csi')])
        return output

    def verify(self,breed,region,records,samples,truth,full,masked):
        key=self.keys[region];folder=self.out/'not-behelim-vs-behelim'/breed/key
        subset=folder/'phased.vcf.gz'
        positions=self.out/'targets'/breed/key/'compatible_positions.bed'
        self.run(breed+'_'+key+'_phased',[self.tools/'bcftools','view','-T',positions,'-Oz','-o',subset,full],[subset])
        before={r[1]:r for r in records};preserved=0;seen=set()
        for path in [full,subset,*([masked] if masked else [])]:
            with gzip.open(path,'rt') as f:
                observed=next(line.rstrip().split('\t')[9:] for line in f if line.startswith('#CHROM'))
            assert observed==samples,'Output sample order changed'
        for r in vcf(subset):
            assert r[1] in before and r[3:5]==before[r[1]][3:5]
            seen.add(r[1]);gtindex=r[8].split(':').index('GT')
            for a,b in zip(before[r[1]][9:],r[9:],strict=True):
                gt=b.split(':')[gtindex];assert gt in ('0|0','0|1','1|0','1|1')
                if '.' not in a:assert count_gt(a)==count_gt(gt);preserved+=1
        assert seen==set(before)
        self.run(breed+'_'+key+'_phaseindex',[self.tools/'bcftools','index','-f',subset],[Path(str(subset)+'.csi')])
        summary=dict(observed_calls_preserved=preserved,all_input_sites_recovered=True)
        if truth:
            metrics=[];pooled=[];seen=set()
            for r in vcf(masked):
                if r[1] not in truth:continue
                t=truth[r[1]];assert r[3:5]==[t['ref'],t['alt']];seen.add(r[1]);fmt=r[8].split(':')
                out=[dict(zip(fmt,g.split(':'))) for g in r[9:]]
                obs=np.array([count_gt(g) for g in t['gt']]);pred=np.array([float(g['DS']) for g in out]);calls=np.array([count_gt(g['GT']) for g in out])
                assert len(out)==len(samples) and np.all(np.isfinite(pred)) and np.all((pred>=0)&(pred<=2)) and np.all(np.isfinite(calls))
                ok=np.isfinite(obs);maf=float(np.nanmean(obs)/2);maf=min(maf,1-maf)
                metrics.append(dict(marker=t['marker'],position=int(r[1]),maf=maf,calls=int(ok.sum()),dosage_r2=correlation(obs,pred),concordance=float(np.mean(obs[ok]==calls[ok])),reported_DR2=next((float(x[4:]) for x in r[7].split(';') if x.startswith('DR2=')),None)))
                pooled.extend(zip(obs[ok],pred[ok],calls[ok]))
            assert seen==set(truth),'Held-out variants missing from imputed output'
            a=np.array(pooled);nonref=(a[:,0]>0)|(a[:,2]>0)
            summary.update(masked_markers=len(truth),called_truth_genotypes=len(a),hard_call_concordance=float(np.mean(a[:,0]==a[:,2])),nonreference_union_concordance=float(np.mean(a[nonref,0]==a[nonref,2])) if nonref.any() else None,pooled_dosage_r2=correlation(a[:,0],a[:,1]),dosage_rmse=float(np.sqrt(np.mean((a[:,0]-a[:,1])**2))))
            summary['frequency_bins']={label:dict(markers=len(part),median_dosage_r2=float(np.median([x['dosage_r2'] for x in part if x['dosage_r2'] is not None])) if any(x['dosage_r2'] is not None for x in part) else None) for label,part in [(label,[v for v in metrics if low<=v['maf']<high]) for label,low,high in [('<0.05',0,.05),('0.05-0.10',.05,.10),('>=0.10',.10,1)]]}
            dump(self.out/'validation'/breed/key/'masked_variant_metrics.private.json',metrics)
        dump(folder/'validation.private.json',summary)
        self.state['jobs'][region]['breeds'][breed]['verification']=summary
        paths=[full,subset,*([masked] if masked else [])]
        self.state['jobs'][region]['breeds'][breed]['outputs']={str(p.relative_to(ROOT)):sha(p) for p in [*paths,*[Path(str(p)+'.csi') for p in paths]]}
        self.save()

    def scheme_a_outputs(self,breed,region,full):
        key=self.keys[region];folder=self.out/'success-vs-behelim'/breed/key;folder.mkdir(parents=True,exist_ok=True)
        ids=folder/'samples.txt';ids.write_text('\n'.join(self.scheme_a[breed])+'\n')
        dest=folder/'imputed.vcf.gz'
        self.run(breed+'_'+key+'_scheme_a',[self.tools/'bcftools','view','-S',ids,'-Oz','-o',dest,full],[dest])
        self.run(breed+'_'+key+'_scheme_a_index',[self.tools/'bcftools','index','-f',dest],[Path(str(dest)+'.csi')])
        with gzip.open(dest,'rt') as f:observed=next(line.rstrip().split('\t')[9:] for line in f if line.startswith('#CHROM'))
        assert observed==self.scheme_a[breed]
        self.state['jobs'][region]['breeds'][breed]['outputs'][str(dest.relative_to(ROOT))]=sha(dest)
        self.state['jobs'][region]['breeds'][breed]['outputs'][str(Path(str(dest)+'.csi').relative_to(ROOT))]=sha(Path(str(dest)+'.csi'))

    def similarity(self,breed,region,records,samples,reference,refids):
        key=self.keys[region];folder=self.out/'targets'/breed/key
        n=min(len(records),self.config['reference_similarity_anchors_per_chromosome'])
        selected=[records[i] for i in np.linspace(0,len(records)-1,n,dtype=int)]
        positions=folder/'similarity_positions.bed';positions.write_text(''.join(f'{r[0]}\t{int(r[1])-1}\t{r[1]}\n' for r in selected))
        text=folder/'reference_anchor_GT.tsv'
        self.capture(breed+'_'+key+'_similarity_GT',[self.tools/'bcftools','query','-R',positions,'-f','%POS[\t%GT]\n',reference],text)
        target={r[1]:r for r in selected};matches=np.zeros((len(samples),len(refids)),dtype=np.int32);called=np.zeros(len(samples),dtype=np.int32);seen=set()
        with text.open() as f:
            for line in f:
                row=line.rstrip().split('\t');pos=row[0];assert pos in target and pos not in seen;seen.add(pos)
                ref=np.array([count_gt(g) for g in row[1:]],dtype=np.float32);assert len(ref)==len(refids) and np.isfinite(ref).all()
                study=np.array([count_gt(g) for g in target[pos][9:]],dtype=np.float32);valid=np.isfinite(study)
                matches+=(study[:,None]==ref[None,:])&valid[:,None];called+=valid
        assert seen==set(target)
        np.savez_compressed(folder/'reference_similarity.private.npz',matches=matches,called=called,samples=samples,reference_samples=refids)

    def summarize_similarity(self):
        report={}
        for breed in self.prefixes:
            total=None;called=None;ids=None;refs=None
            for region in self.regions:
                p=self.out/'targets'/breed/self.keys[region]/'reference_similarity.private.npz'
                with np.load(p) as data:
                    if total is None:total=data['matches'].copy();called=data['called'].copy();ids=data['samples'].tolist();refs=data['reference_samples'].tolist()
                    else:
                        assert ids==data['samples'].tolist() and refs==data['reference_samples'].tolist()
                        total+=data['matches'];called+=data['called']
            assert np.all(called>0)
            best=total.argmax(axis=1);values=total[np.arange(len(ids)),best]/called
            enough=called>=self.config['reference_similarity_min_genotypes']
            flags=enough&(values>=self.config['reference_similarity_flag_concordance'])
            private=[dict(study_id=ids[i],reference_id=refs[best[i]],called_anchors=int(called[i]),maximum_concordance=float(values[i]),flag=bool(flags[i])) for i in range(len(ids))]
            dump(self.out/(breed+'_reference_similarity.private.json'),private)
            report[breed]=dict(dogs=len(ids),dogs_with_enough_anchors=int(enough.sum()),flagged_dogs=int(flags.sum()),minimum_anchor_calls=int(called.min()),maximum_concordance=float(values.max()))
        dump(self.out/'reference_similarity_summary.private.json',dict(breeds=report,interpretation='Similarity flags prompt identity review. This screen does not exclude close relatives or establish independent reference membership.'))
        self.state['reference_similarity']=report

    def execute(self):
        self.effective_reference=self.args.reference_source or str(ensure_panel())
        if not self.args.reference_source:self.state['reference_panel']=json.loads((ROOT/'data/02-reference/dog10k/panel_assets.json').read_text())
        self.state['status']='running';self.save()
        for region in self.regions:
            if self.state['jobs'][region]['status']=='complete':
                stored={p:h for d in self.state['jobs'][region]['breeds'].values() for p,h in d['outputs'].items()}
                assert all(sha(ROOT/p)==h for p,h in stored.items());continue
            print('Processing',region,flush=True)
            reference,good,ids=self.reference(region)
            for breed in self.prefixes:
                records,samples,truth=self.targets(breed,region,good,ids)
                full=self.beagle(breed,region,'full',reference)
                masked=self.beagle(breed,region,'masked',reference) if truth else None
                self.verify(breed,region,records,samples,truth,full,masked)
                self.scheme_a_outputs(breed,region,full)
                self.similarity(breed,region,records,samples,reference,ids)
                print(breed,region,'phasing, imputation, and call checks finished.',flush=True)
            self.state['jobs'][region]['status']='complete';self.save()
        self.summarize_similarity()
        self.verify_sources();self.state.update(status='complete',source_hashes_unchanged=True,association_ready=False);self.save()
        print('Selected plan complete. Scientific validation and reference-independence review precede association use.',flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',type=Path,default=ROOT/'config/imputation.json')
    ap.add_argument('--output',type=Path,default=ROOT/'data/04-imputation/autosomal-20261003')
    ap.add_argument('--tools',type=Path,default=Path.home()/'anaconda3/envs/tsa2024-imputation/bin')
    ap.add_argument('--regions',nargs='+');ap.add_argument('--reference-source')
    group=ap.add_mutually_exclusive_group(required=True);group.add_argument('--prepare-only',action='store_true');group.add_argument('--execute',action='store_true')
    args=ap.parse_args();assert os.name!='nt','Run this workflow in Linux/WSL'
    for name in ('python','bcftools','java'):assert (args.tools/name).exists()
    wf=Workflow(args)
    try:
        wf.prepare()
        if args.execute:wf.execute()
    except Exception as e:
        wf.state.update(status='failed',error=str(e));wf.save();raise


if __name__=='__main__':main()
