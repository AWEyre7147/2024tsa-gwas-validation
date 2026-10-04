"""Independent checks of source calls, PCA subspace, embeddings and notebook images."""
from pathlib import Path
import base64
import csv
import hashlib
import json
import os
import sys
os.environ.setdefault('NUMBA_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import numpy as np
import nbformat
from umap import UMAP
from sklearn.manifold import TSNE
from threadpoolctl import threadpool_limits
from array_qc import verify_calls

def rows(p):return [r.split() for r in p.read_text().splitlines()]
def packed(p,n,m):return np.frombuffer(p.with_suffix('.bed').read_bytes(),np.uint8,offset=3).reshape(m,(n+3)//4)
def main():
    root=Path(__file__).resolve().parents[1]
    out=root/'data/structure/20261002-pca-umap-tsne'
    manifest=json.loads((out/'manifest.private.json').read_text())
    checks={}
    for scope in ['combined','GSP','LR']:
        p=out/scope/'02_pruned_genotypes'; fam=rows(p.with_suffix('.fam'));bim=rows(p.with_suffix('.bim'));n=len(fam);m=len(bim)
        j=np.arange(n);calls=(packed(p,n,m)[:,j//4]>>(2*(j%4)).astype(np.uint8))&3
        # The pooled merge must reproduce every call and allele from each breed.
        for breed in ([scope] if scope!='combined' else ['GSP','LR']):
            source=root/f'data/qc/not-behelim-vs-behelim/{breed}/2024TSA_{breed}'
            if scope!='combined':verify_calls(source,p)
            else:
                sf=rows(source.with_suffix('.fam'));sb=rows(source.with_suffix('.bim'));vi={r[1]:i for i,r in enumerate(sb)};si={r[1]:i for i,r in enumerate(sf)}
                samplepos=[i for i,r in enumerate(fam) if r[1] in si];inds=np.array([si[fam[i][1]] for i in samplepos]);vind=np.array([vi[r[1]] for r in bim])
                assert all(r==sb[vi[r[1]]] for r in bim)
                original=packed(source,len(sf),len(sb))[vind]
                expected=(original[:,inds//4]>>(2*(inds%4)).astype(np.uint8))&3
                assert np.array_equal(expected,calls[:,samplepos])
        geno=np.array([2,np.nan,1,0])[calls];mean=np.nanmean(geno,axis=1);freq=mean/2
        maf=np.minimum(freq,1-freq)
        assert maf.min()>=.025
        z=np.where(np.isnan(geno),0,(geno-mean[:,None])/np.sqrt(2*freq*(1-freq))[:,None])
        grm=z.T@z/m
        vals,vec=np.linalg.eigh(grm);vec=vec[:,::-1]
        pcrows=rows(out/scope/'03_pca.eigenvec')
        assert [r[1] for r in pcrows[1:]]==[r[1] for r in fam]
        pcs=np.array([[float(v) for v in r[2:]] for r in pcrows[1:]])
        basis=np.linalg.qr(pcs[:,:20])[0]
        similarity=np.linalg.svd(vec[:,:20].T@basis,compute_uv=False).min()
        # Small differences are expected from PLINK's missingness corrections.
        assert similarity>.98, f'{scope}: independent PCA disagrees'
        checks[scope]=dict(source_calls_verified=True,minimum_maf=float(maf.min()),monomorphic_count=int((maf==0).sum()),
                           independent_first20_subspace_min_similarity=float(similarity),
                           leading_PC_correlations=[float(abs(np.corrcoef(vec[:,k],pcs[:,k])[0,1])) for k in range(2)])
    scope='LR';path=out/scope/'coordinates.private.tsv'
    with path.open(encoding='utf-8',newline='') as f:
        data=list(csv.DictReader(f,delimiter='\t'))
    x=np.array([[float(r[f'PC{k+1}']) for k in range(20)] for r in data])
    old_u=np.array([[float(r['UMAP1']),float(r['UMAP2'])] for r in data]);old_t=np.array([[float(r['TSNE1']),float(r['TSNE2'])] for r in data])
    with threadpool_limits(limits=1):
        u=UMAP(n_components=2,n_neighbors=15,min_dist=.1,metric='euclidean',random_state=20261002,n_jobs=1).fit_transform(x)
        t=TSNE(n_components=2,perplexity=15,early_exaggeration=12,learning_rate='auto',max_iter=1500,init='pca',metric='euclidean',random_state=20261002,method='barnes_hut',n_jobs=1).fit_transform(x)
    assert np.array_equal(u,old_u) and np.array_equal(t,old_t)
    checks['repeated_LR_embeddings']={'umap_identical':True,'tsne_identical':True}
    public=nbformat.read(root/'notebooks/03_structure_and_model_design.ipynb',as_version=4)
    review=nbformat.read(root/'results/private/03_structure/03_structure_and_model_design.ipynb',as_version=4)
    nbformat.validate(public);nbformat.validate(review)
    assert all(c.get('outputs',[])==[] and c.get('execution_count') is None for c in public.cells if c.cell_type=='code')
    attachments={name:base64.b64decode(data['image/png']) for c in review.cells for name,data in c.get('attachments',{}).items()}
    assert len(attachments)==9
    for name,data in attachments.items():assert data==(root/'results/private/03_structure/figures'/name).read_bytes()
    for path,h in manifest['figures'].items():assert hashlib.sha256((root/path).read_bytes()).hexdigest()==h
    checks['notebook']={'embedded_figures':9,'public_outputs_cleared':True,'private_attachments_verified':True}
    (out/'verification.private.json').write_text(json.dumps(checks,indent=2)+'\n')
    print(json.dumps(checks,indent=2))

if __name__=='__main__':main()
