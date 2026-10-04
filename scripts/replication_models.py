"""Run paired measured-target mixed models and a limited null-calibration screen."""
from pathlib import Path
import argparse
import csv
import json
import os
import subprocess
import sys
import numpy as np
from notebook_support import read_json, file_digest, study_path, verify_datasets

ROOT=Path(__file__).resolve().parents[1]


def rows(path):return [line.split() for line in path.read_text(encoding='utf-8').splitlines() if line.strip()]


def dump(path,value):path.write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8')


def standardized_background(prefix, aligned_bim):
    """Standardize observed background alleles; missing background calls become centered zero."""
    fam=rows(prefix.with_suffix('.fam'));bim=rows(prefix.with_suffix('.bim'))
    n=len(fam);m=len(bim);raw=prefix.with_suffix('.bed').read_bytes()
    assert raw[:3]==bytes([108,27,1]);stride=(n+3)//4
    packed=np.frombuffer(raw,np.uint8,offset=3).reshape(m,stride)
    indices=np.arange(n)
    codes=(packed[:,indices//4]>>(2*(indices%4)))&3
    genotypes=np.array([2,np.nan,1,0])[codes]
    frequency=np.nanmean(genotypes,axis=1)/2
    locations={row[1]:int(row[0]) for row in rows(aligned_bim)}
    # Use retained aligned locations for chromosome exclusion, even when a locus moved.
    polymorphic=(frequency>0)&(frequency<1)&np.array([row[1] in locations for row in bim])
    genotypes=genotypes[polymorphic];frequency=frequency[polymorphic]
    z=(genotypes-2*frequency[:,None])/np.sqrt(2*frequency*(1-frequency))[:,None]
    z=np.nan_to_num(z,nan=0)
    chromosomes=np.array([locations[row[1]] for row,ok in zip(bim,polymorphic) if ok])
    return [row[1] for row in fam],chromosomes,z


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,default=ROOT/'data/replication/models-aligned-20261004')
    ap.add_argument('--tools',type=Path,default=Path.home()/'anaconda3/envs/tsa2024-association/bin')
    ap.add_argument('--calibration-replicates',type=int)
    args=ap.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    cfg=read_json(ROOT/'config/replication.json')
    if args.calibration_replicates is not None:cfg['calibration_replicates']=args.calibration_replicates
    assert cfg['calibration_replicates']>=0
    prep=ROOT/'data/replication/target-preparation-20261004'
    alignment=read_json(ROOT/'data/aligned/manifest.private.json');verify_datasets(ROOT,alignment)
    background=ROOT/'data/structure/20261002-pca-umap-tsne/LR/02_pruned_genotypes'
    aligned_bim=ROOT/'data/aligned/not-behelim-vs-behelim/LR/2024TSA_LR.bim'
    source={str(p.relative_to(ROOT)):file_digest(p) for p in [*background.parent.glob('02_pruned_genotypes.b*'),background.with_suffix('.fam'),aligned_bim,prep/'allele_counts.private.tsv',*prep.glob('*_dosages.private.tsv')]}
    versions=dict(gemma=subprocess.check_output([str(args.tools/'gemma'),'-h'],text=True).splitlines()[:4],
        R=subprocess.check_output([str(args.tools/'Rscript'),'--version'],text=True).strip(),
        GMMAT=subprocess.check_output([str(args.tools/'Rscript'),'-e','cat(as.character(packageVersion("GMMAT")))'],text=True).strip(),
        gemma_sha256=file_digest(args.tools/'gemma'))
    plan=dict(config=cfg,sources=source,script_sha256=file_digest(Path(__file__)),
              R_script_sha256=file_digest(ROOT/'scripts/replication_models.R'),versions=versions)
    planfile=out/'plan.private.json'
    if planfile.exists():assert read_json(planfile)==plan,'Plan changed; use a new model output directory'
    else:dump(planfile,plan)
    ids,chromosomes,z=standardized_background(background,aligned_bim);idpos={value:i for i,value in enumerate(ids)}
    with (prep/'allele_counts.private.tsv').open(encoding='utf-8',newline='') as handle:
        inventory=list(csv.DictReader(handle,delimiter='\t'))
    jobs=[]
    for scenario in ['success-vs-behelim','not-behelim-vs-behelim']:
        with (prep/(scenario+'_dosages.private.tsv')).open(encoding='utf-8',newline='') as handle:
            dogs=list(csv.DictReader(handle,delimiter='\t'))
        for number,target in enumerate([r for r in inventory if r['scenario']==scenario]):
            marker=target['marker'];chrom=int(target['canFam4_chromosome'][3:])
            folder=out/scenario/marker;folder.mkdir(parents=True,exist_ok=True)
            complete=[row for row in dogs if row[marker]!='']
            indices=np.array([idpos[row['study_id']] for row in complete])
            outside=chromosomes!=chrom;matrix=z[outside][:,indices]
            kinship=matrix.T@matrix/len(matrix)
            assert np.allclose(kinship,kinship.T) and np.linalg.eigvalsh(kinship).min()>-1e-8
            np.savetxt(folder/'kinship.txt',kinship,fmt='%.16g',delimiter='\t')
            with (folder/'dogs.private.tsv').open('w',encoding='utf-8',newline='') as handle:
                writer=csv.writer(handle,delimiter='\t');writer.writerow(['id','y','male','genotype'])
                writer.writerows([row['study_id'],int(row['phenotype'])-1,int(row['sex']=='1'),row[marker]] for row in complete)
            (folder/'genotypes.txt').write_text(marker+'\t'+'\t'.join(row[marker] for row in complete)+'\n',encoding='utf-8')
            (folder/'gemma_genotypes.txt').write_text(','.join([marker,target['alternate_allele'],target['reference_allele'],*[row[marker] for row in complete]])+'\n',encoding='utf-8')
            (folder/'covariates.txt').write_text(''.join(f"1\t{int(row['sex']=='1')}\n" for row in complete),encoding='utf-8')
            np.savetxt(folder/'phenotype.txt',[int(row['phenotype'])-1 for row in complete],fmt='%d')
            metadata=dict(target,complete_dogs=len(complete),case_count=sum(row['phenotype']=='2' for row in complete),
                          background_markers=len(matrix),kinship_min_eigenvalue=float(np.linalg.eigvalsh(kinship).min()),
                          seed=cfg['seed']+number+(1000 if scenario=='not-behelim-vs-behelim' else 0),
                          replicates=cfg['calibration_replicates'],gemma=str(args.tools/'gemma'))
            dump(folder/'input.private.json',metadata)
            done=folder/'result.private.json'
            receipt=folder/'result-receipt.private.json'
            if done.exists():
                assert receipt.exists() and read_json(receipt)['sha256']==file_digest(done),'Result changed or lacks a receipt; use a new output directory'
            if not done.exists():
                print(f'{scenario}: {marker}, {len(complete)} dogs; running both models and null screen',flush=True)
                command=[str(args.tools/'Rscript'),str(ROOT/'scripts/replication_models.R'),str(folder)]
                env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
                with (folder/'model.console.txt').open('w',encoding='utf-8') as log:
                    subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,env=env,check=True)
                dump(receipt,dict(sha256=file_digest(done)))
            jobs.append(dict(scenario=scenario,marker=marker,result=str(done.relative_to(ROOT))))
    assert all(file_digest(study_path(ROOT,p))==h for p,h in source.items())
    dump(out/'manifest.private.json',dict(status='paired_models_complete',association_ready_for_claims=False,
        historical_direction_pending=True,config=cfg,jobs=jobs,source_hashes_unchanged=True))
    print('Both models and their initial null screen completed at 12 targets in each cohort.',flush=True)


if __name__=='__main__':main()
