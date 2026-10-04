"""Paired measured-genotype GWAS, separated by breed, phenotype scheme and chromosome."""
from pathlib import Path
import argparse
import csv
import json
import os
import subprocess
import numpy as np
from notebook_support import read_json, file_digest, study_path, verify_datasets
from replication_models import standardized_background, rows

ROOT=Path(__file__).resolve().parents[1]


def dump(path,value):
    path.write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8')


def decode(prefix):
    fam=rows(prefix.with_suffix('.fam'));bim=rows(prefix.with_suffix('.bim'))
    raw=prefix.with_suffix('.bed').read_bytes();assert raw[:3]==b'\x6c\x1b\x01'
    packed=np.frombuffer(raw,np.uint8,offset=3).reshape(len(bim),(len(fam)+3)//4)
    ix=np.arange(len(fam));codes=(packed[:,ix//4]>>(2*(ix%4)))&3
    return fam,bim,np.array([2,np.nan,1,0])[codes]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'data/gwas/measured-20261004')
    parser.add_argument('--tools',type=Path,default=Path.home()/'anaconda3/envs/tsa2024-association/bin')
    parser.add_argument('--chromosomes',type=int,nargs='+',default=list(range(1,39)))
    args=parser.parse_args();assert set(args.chromosomes)<=set(range(1,39))
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    cfg=read_json(ROOT/'config/gwas.json')
    alignment=read_json(ROOT/'data/aligned/manifest.private.json');verify_datasets(ROOT,alignment)
    sources={}
    for dataset in alignment['datasets']:
        for path,digest in dataset['files'].items():
            if Path(path).suffix in ['.bed','.bim','.fam']:sources[path]=digest
    for breed in cfg['breeds']:
        bg=ROOT/'data/structure/20261002-pca-umap-tsne'/breed/'02_pruned_genotypes'
        sources.update({str(bg.with_suffix(s).relative_to(ROOT)):file_digest(bg.with_suffix(s)) for s in ['.bed','.bim','.fam']})
    versions=dict(GEMMA=subprocess.check_output([str(args.tools/'gemma'),'-h'],text=True).splitlines()[:4],
        R=subprocess.check_output([str(args.tools/'Rscript'),'--version'],text=True).strip(),
        GMMAT=subprocess.check_output([str(args.tools/'Rscript'),'-e','cat(as.character(packageVersion("GMMAT")))'],text=True).strip(),
        gemma_sha256=file_digest(args.tools/'gemma'))
    plan=dict(config=cfg,chromosomes=args.chromosomes,sources=sources,versions=versions,
        scripts={name:file_digest(ROOT/'scripts'/name) for name in ['breed_gwas.py','breed_gwas.R','replication_models.py','notebook_support.py']})
    if (out/'plan.private.json').exists():assert read_json(out/'plan.private.json')==plan,'Run specification changed; use a new output folder'
    else:dump(out/'plan.private.json',plan)
    jobs=[];groups=[]
    for breed in cfg['breeds']:
        broad=ROOT/'data/aligned'/cfg['sensitivity_cohort']/breed/('2024TSA_'+breed)
        background=ROOT/'data/structure/20261002-pca-umap-tsne'/breed/'02_pruned_genotypes'
        bg_ids,bg_chromosomes,z=standardized_background(background,broad.with_suffix('.bim'))
        idpos={value:i for i,value in enumerate(bg_ids)}
        for scenario in [cfg['primary_cohort'],cfg['sensitivity_cohort']]:
            prefix=ROOT/'data/aligned'/scenario/breed/('2024TSA_'+breed)
            fam,bim,g=decode(prefix)
            assert all(r[4] in ['1','2'] and r[5] in ['1','2'] for r in fam)
            if scenario==cfg['primary_cohort']:
                broadfam=rows(broad.with_suffix('.fam'));assert {r[1] for r in fam if r[5]=='2'}=={r[1] for r in broadfam if r[5]=='2'}
            indices=np.array([idpos[r[1]] for r in fam]);case=np.array([r[5]=='2' for r in fam]);control=~case
            af=np.nanmean(g,axis=1)/2;maf=np.minimum(af,1-af)
            eligible=np.isfinite(maf)&(maf>=cfg['minimum_cohort_maf'])
            groups.append(dict(breed=breed,scenario=scenario,dogs=len(fam),cases=int(case.sum()),controls=int(control.sum()),
                candidate_markers=len(bim),eligible_markers=int(eligible.sum()),aligned_background_markers=len(z)))
            for chrom in args.chromosomes:
                folder=out/breed/scenario/f'chr{chrom:02d}';folder.mkdir(parents=True,exist_ok=True)
                done=folder/'receipt.private.json'
                if done.exists():
                    receipt=read_json(done)
                    assert all(file_digest(folder/p)==h for p,h in receipt['files'].items()),'Completed output changed'
                    jobs.append(dict(breed=breed,scenario=scenario,chromosome=chrom,folder=str(folder.relative_to(ROOT))));continue
                selected=[i for i,r in enumerate(bim) if int(r[0])==chrom and eligible[i]]
                matrix=z[bg_chromosomes!=chrom][:,indices];K=matrix.T@matrix/len(matrix)
                assert np.linalg.eigvalsh(K).min()>-1e-8
                np.savetxt(folder/'kinship.txt',K,fmt='%.16g',delimiter='\t')
                with (folder/'dogs.private.tsv').open('w',encoding='utf-8',newline='') as h:
                    writer=csv.writer(h,delimiter='\t');writer.writerow(['id','y','male']);writer.writerows([r[1],int(r[5])-1,int(r[4]=='1')] for r in fam)
                np.savetxt(folder/'phenotype.txt',case.astype(int),fmt='%d')
                np.savetxt(folder/'covariates.txt',np.column_stack([np.ones(len(fam)),[int(r[4]=='1') for r in fam]]),fmt='%d',delimiter='\t')
                with (folder/'markers.private.tsv').open('w',encoding='utf-8',newline='') as mh, (folder/'genotypes.txt').open('w',encoding='utf-8') as gh, (folder/'gemma_genotypes.txt').open('w',encoding='utf-8') as lh:
                    writer=csv.writer(mh,delimiter='\t');writer.writerow(['marker','chromosome','position','ALT','REF','AF','MAF','called_dogs','case_called','control_called','case_alt_copies','control_alt_copies','case_0','case_1','case_2','control_0','control_1','control_2'])
                    for i in selected:
                        r=bim[i];observed=g[i];filled=np.where(np.isnan(observed),2*af[i],observed);values=[format(v,'.16g') for v in filled]
                        gh.write(r[1]+'\t'+'\t'.join(values)+'\n');lh.write(','.join([r[1],r[4],r[5],*values])+'\n')
                        writer.writerow([r[1],r[0],r[3],r[4],r[5],af[i],maf[i],int(np.isfinite(observed).sum()),int(np.isfinite(observed[case]).sum()),int(np.isfinite(observed[control]).sum()),int(np.nansum(observed[case])),int(np.nansum(observed[control])),*[int(np.sum(observed[mask]==count)) for mask in [case,control] for count in [0,1,2]]])
                dump(folder/'input.private.json',dict(breed=breed,scenario=scenario,chromosome=chrom,dogs=len(fam),cases=int(case.sum()),markers=len(selected),background_markers=len(matrix),gemma=str(args.tools/'gemma')))
                print(f'{breed} {scenario} chromosome {chrom}/38: {len(selected)} measured markers; both models',flush=True)
                env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
                with (folder/'model.console.txt').open('w',encoding='utf-8') as log:
                    subprocess.run([str(args.tools/'Rscript'),str(ROOT/'scripts/breed_gwas.R'),str(folder)],stdout=log,stderr=subprocess.STDOUT,check=True,env=env)
                files=[p for p in folder.rglob('*') if p.is_file()]
                dump(done,dict(files={str(p.relative_to(folder)):file_digest(p) for p in files}))
                jobs.append(dict(breed=breed,scenario=scenario,chromosome=chrom,folder=str(folder.relative_to(ROOT))))
    assert all(file_digest(study_path(ROOT,p))==h for p,h in sources.items())
    dump(out/'manifest.private.json',dict(status='paired_measured_scans_complete',chromosomes=args.chromosomes,
        groups=groups,jobs=jobs,config=cfg,source_hashes_unchanged=True,discovery_claims_validated=False))
    print('Paired measured scans complete.',flush=True)


if __name__=='__main__':main()
