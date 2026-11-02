"""Independently reconstruct pooled LOCO matrices and spot-check numerical dosage inputs."""
from pathlib import Path
import csv,json
import numpy as np
from breed_gwas import decode
from replication_models import rows
from notebook_support import read_json,study_path,file_digest

ROOT=Path(__file__).resolve().parents[1]
run=ROOT/'data/07-gwas/pooled-measured-20261004';manifest=read_json(run/'manifest.private.json');plan=read_json(run/'plan.private.json')
assert all(file_digest(study_path(ROOT,p))==h for p,h in plan['sources'].items())
pruned=set();blocks=[]
for breed in ['LR','GSP']:
    pruned.update(r[1] for r in rows(ROOT/f'data/03-structure/20261002-pca-umap-tsne/{breed}/02_pruned_genotypes.bim'))
    fam,bim,g=decode(ROOT/f'data/02-aligned/not-behelim-vs-behelim/{breed}/2024TSA_{breed}');blocks.append((fam,bim,g,{r[1]:i for i,r in enumerate(bim)}))
keys=sorted(pruned&set(blocks[0][3])&set(blocks[1][3]));standard=[];keep=np.ones(len(keys),bool)
for fam,bim,g,mp in blocks:
    genotype=g[[mp[k] for k in keys]];mu=np.nanmean(genotype,axis=1);frequency=mu/2;scale=np.sqrt(2*frequency*(1-frequency));keep&=(scale>0)
    with np.errstate(divide='ignore',invalid='ignore'):
        centered=np.where(np.isfinite(genotype),genotype-mu[:,None],0);standard.append(np.divide(centered,scale[:,None],out=np.zeros_like(centered),where=scale[:,None]>0))
Z=np.column_stack(standard)[keep];keys=np.array(keys)[keep];chrom=np.array([int(blocks[0][1][blocks[0][3][k]][0]) for k in keys])
assert len(keys)==plan['background_markers'];ids=[r[1] for block in blocks for r in block[0]];mp={k:i for i,k in enumerate(ids)}
total=Z.T@Z;sums={c:Z[chrom==c].T@Z[chrom==c] for c in range(1,39)};checks=[]
cohorts={}
for scenario in plan['config']['schemes']:
    cohorts[scenario]=[decode(ROOT/f'data/02-aligned/{scenario}/{breed}/2024TSA_{breed}') for breed in ['LR','GSP']]
for job in manifest['jobs']:
    folder=study_path(ROOT,job['folder']);dogs=list(csv.DictReader((folder/'dogs.private.tsv').open(encoding='utf-8'),delimiter='\t'));ix=[mp[d['id']] for d in dogs]
    c=job['chromosome'];expected=((total-sums[c])/(len(Z)-int((chrom==c).sum())))[np.ix_(ix,ix)];actual=np.loadtxt(folder/'kinship.txt')
    assert np.max(np.abs(expected-actual))<1e-10
    data=cohorts[job['scenario']];expected_ids=[r[1] for block in data for r in block[0]];assert expected_ids==[d['id'] for d in dogs]
    markers=list(csv.DictReader((folder/'markers.private.tsv').open(encoding='utf-8'),delimiter='\t'));selected={0,len(markers)//2,len(markers)-1}
    with (folder/'genotypes.txt').open(encoding='utf-8') as h:
        for i,line in enumerate(h):
            if i not in selected:continue
            fields=line.split();marker=markers[i]['marker'];assert fields[0]==marker;parts=[]
            for fam,bim,g in data:
                index={r[1]:j for j,r in enumerate(bim)};r=bim[index[marker]];assert r[4]==markers[i]['ALT'] and r[5]==markers[i]['REF']
                observed=g[index[marker]];assert min(np.nanmean(observed)/2,1-np.nanmean(observed)/2)>=.025
                parts.append(np.where(np.isnan(observed),np.nanmean(observed),observed))
            assert np.max(np.abs(np.r_[*parts]-np.array(fields[1:],float)))<1e-12
    checks.append(dict(folder=job['folder'],relationship_verified=True,selected_dosages_verified=len(selected)))
(run/'input-check.private.json').write_text(json.dumps(dict(status='pooled_relationships_and_selected_inputs_verified',checks=checks,background_markers=len(Z)),indent=2)+'\n',encoding='utf-8')
print(len(checks),'pooled matrices and selected dosage inputs independently verified.')
