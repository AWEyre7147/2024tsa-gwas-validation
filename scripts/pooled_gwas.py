"""Exploratory pooled measured GWAS with breed/sex adjustment and one residual-genetic LOCO matrix."""
from pathlib import Path
import argparse,csv,json,os,subprocess
import numpy as np
from breed_gwas import decode,dump
from cross_breed_meta import aligned_family
from replication_models import rows
from notebook_support import read_json,file_digest,study_path,verify_datasets

ROOT=Path(__file__).resolve().parents[1]

def run(output,tools,chromosomes):
    cfg=read_json(ROOT/'config/cross_breed.json');alignment=read_json(ROOT/'data/02-aligned/manifest.private.json');verify_datasets(ROOT,alignment)
    sources={p:h for d in alignment['datasets'] for p,h in d['files'].items() if Path(p).suffix in ['.bed','.bim','.fam']}
    pruned=set();broad=[]
    for breed in ['LR','GSP']:
        prefix=ROOT/f'data/03-structure/20261002-pca-umap-tsne/{breed}/02_pruned_genotypes'
        pruned.update(r[1] for r in rows(prefix.with_suffix('.bim')))
        sources[str(prefix.with_suffix('.bim').relative_to(ROOT))]=file_digest(prefix.with_suffix('.bim'))
        fam,bim,g=decode(ROOT/f'data/02-aligned/not-behelim-vs-behelim/{breed}/2024TSA_{breed}')
        broad.append((fam,bim,g,{r[1]:i for i,r in enumerate(bim)}))
    common=aligned_family('not-behelim-vs-behelim');background=sorted(common&pruned)
    zparts=[];background_keep=np.ones(len(background),bool)
    for fam,bim,g,index in broad:
        measured=g[[index[k] for k in background]];af=np.nanmean(measured,axis=1)/2
        background_keep&=(af>0)&(af<1)
        with np.errstate(invalid='ignore',divide='ignore'):
            zparts.append(np.nan_to_num((measured-2*af[:,None])/np.sqrt(2*af*(1-af))[:,None],nan=0,posinf=0,neginf=0))
    bg=np.array(background)[background_keep];z=np.column_stack(zparts)[background_keep]
    bgchrom=np.array([int(broad[0][1][broad[0][3][k]][0]) for k in bg])
    broad_ids=[r[1] for item in broad for r in item[0]];indexdog={k:i for i,k in enumerate(broad_ids)}
    assert len(indexdog)==len(broad_ids)
    versions=dict(GEMMA=subprocess.check_output([str(tools/'gemma'),'-h'],text=True).splitlines()[:4],GMMAT=subprocess.check_output([str(tools/'Rscript'),'-e','cat(as.character(packageVersion("GMMAT")))'],text=True).strip(),gemma_sha256=file_digest(tools/'gemma'))
    plan=dict(config=cfg,chromosomes=chromosomes,sources=sources,versions=versions,scripts={n:file_digest(ROOT/'scripts'/n) for n in ['pooled_gwas.py','pooled_gwas.R','breed_gwas.py','cross_breed_meta.py','replication_models.py','notebook_support.py']},background_markers=len(bg))
    output.mkdir(parents=True,exist_ok=True)
    if (output/'plan.private.json').exists():assert read_json(output/'plan.private.json')==plan,'Specification changed: choose a new output folder'
    else:dump(output/'plan.private.json',plan)
    jobs=[];groups=[]
    for scenario in cfg['schemes']:
        shared=aligned_family(scenario);blocks=[];maps=[]
        for breed in ['LR','GSP']:
            fam,bim,g=decode(ROOT/f'data/02-aligned/{scenario}/{breed}/2024TSA_{breed}')
            blocks.append((fam,bim,g));maps.append({r[1]:i for i,r in enumerate(bim)})
        keys=sorted(shared,key=lambda k:(int(blocks[0][1][maps[0][k]][0]),int(blocks[0][1][maps[0][k]][3]),k))
        arrays=[g[[mp[k] for k in keys]] for (_,_,g),mp in zip(blocks,maps)]
        afs=[np.nanmean(g,axis=1)/2 for g in arrays];eligible=np.logical_and.reduce([np.isfinite(af)&(np.minimum(af,1-af)>=.025) for af in afs])
        filled=np.column_stack([np.where(np.isnan(g),2*af[:,None],g) for g,af in zip(arrays,afs)])
        observed=np.column_stack(arrays);fam=[r for block in blocks for r in block[0]];gsp=np.r_[np.zeros(len(blocks[0][0]),int),np.ones(len(blocks[1][0]),int)]
        assert all(r[4] in ['1','2'] and r[5] in ['1','2'] for r in fam)
        case=np.array([r[5]=='2' for r in fam]);indices=[indexdog[r[1]] for r in fam]
        chrom=np.array([int(blocks[0][1][maps[0][k]][0]) for k in keys]);groups.append(dict(breed='pooled',scenario=scenario,dogs=len(fam),cases=int(case.sum()),controls=int((~case).sum()),candidate_markers=len(keys),eligible_markers=int(eligible.sum()),background_markers=len(bg)))
        for ch in chromosomes:
            folder=output/scenario/f'chr{ch:02d}';folder.mkdir(parents=True,exist_ok=True);receipt=folder/'receipt.private.json'
            job=dict(breed='pooled',scenario=scenario,chromosome=ch,folder=str(folder.relative_to(ROOT)))
            if receipt.exists():
                assert all(file_digest(folder/p)==h for p,h in read_json(receipt)['files'].items());jobs.append(job);continue
            matrix=z[bgchrom!=ch][:,indices];K=matrix.T@matrix/len(matrix)
            assert np.linalg.eigvalsh(K).min()>-1e-8;assert np.max(np.abs(K-K.T))<1e-12
            np.savetxt(folder/'kinship.txt',K,fmt='%.16g',delimiter='\t')
            with (folder/'dogs.private.tsv').open('w',encoding='utf-8',newline='') as h:
                w=csv.writer(h,delimiter='\t');w.writerow(['id','y','male','gsp']);w.writerows([r[1],int(r[5])-1,int(r[4]=='1'),int(b)] for r,b in zip(fam,gsp))
            np.savetxt(folder/'phenotype.txt',case.astype(int),fmt='%d');X=np.column_stack([np.ones(len(fam)),[int(r[4]=='1') for r in fam],gsp]);assert np.linalg.matrix_rank(X)==3
            np.savetxt(folder/'covariates.txt',X,fmt='%d',delimiter='\t')
            selected=np.where(eligible&(chrom==ch))[0]
            with (folder/'markers.private.tsv').open('w',encoding='utf-8',newline='') as mh,(folder/'genotypes.txt').open('w',encoding='utf-8') as gh,(folder/'gemma_genotypes.txt').open('w',encoding='utf-8') as lh:
                w=csv.writer(mh,delimiter='\t');w.writerow(['marker','chromosome','position','ALT','REF','AF','MAF','called_dogs','case_called','control_called','case_alt_copies','control_alt_copies','case_0','case_1','case_2','control_0','control_1','control_2'])
                for i in selected:
                    r=blocks[0][1][maps[0][keys[i]]];a=observed[i];values=[format(v,'.16g') for v in filled[i]];frequency=filled[i].mean()/2
                    gh.write(keys[i]+'\t'+'\t'.join(values)+'\n');lh.write(','.join([keys[i],r[4],r[5],*values])+'\n')
                    w.writerow([keys[i],r[0],r[3],r[4],r[5],frequency,min(frequency,1-frequency),int(np.isfinite(a).sum()),int(np.isfinite(a[case]).sum()),int(np.isfinite(a[~case]).sum()),float(np.nansum(a[case])),float(np.nansum(a[~case])),*[int(np.sum(a[m]==c)) for m in [case,~case] for c in [0,1,2]]])
            dump(folder/'input.private.json',dict(breed='pooled',scenario=scenario,chromosome=ch,dogs=len(fam),cases=int(case.sum()),markers=len(selected),background_markers=len(matrix),gemma=str(tools/'gemma'),minimum_eigenvalue=float(np.linalg.eigvalsh(K).min()),cross_breed_mean=float(K[np.ix_(gsp==0,gsp==1)].mean())))
            print(f'Exploratory pooled {scenario} chr{ch}/38: {len(selected)} markers, {len(fam)} dogs',flush=True)
            with (folder/'model.console.txt').open('w',encoding='utf-8') as h:
                subprocess.run([str(tools/'Rscript'),str(ROOT/'scripts/pooled_gwas.R'),str(folder)],stdout=h,stderr=subprocess.STDOUT,check=True,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
            dump(receipt,dict(files={str(p.relative_to(folder)):file_digest(p) for p in folder.rglob('*') if p.is_file()}));jobs.append(job)
    assert all(file_digest(study_path(ROOT,p))==h for p,h in sources.items())
    dump(output/'manifest.private.json',dict(status='exploratory_pooled_scans_complete',config=cfg,groups=groups,jobs=jobs,source_hashes_unchanged=True,publication_inclusion_undecided=True))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'data/07-gwas/pooled-measured-20261004');p.add_argument('--tools',type=Path,default=Path.home()/'anaconda3/envs/tsa2024-association/bin');p.add_argument('--chromosomes',nargs='+',type=int,default=list(range(1,39)))
    a=p.parse_args();assert set(a.chromosomes)<=set(range(1,39));run(a.output,a.tools,a.chromosomes)
