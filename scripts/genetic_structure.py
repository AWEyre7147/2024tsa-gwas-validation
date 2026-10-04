"""PLINK PCA/IBD and descriptive UMAP/t-SNE for the unique QC study dogs.

Run in the isolated environment with requirements-structure.txt installed.
Images and row-level results remain private. No association or automatic
sample exclusion is performed. A completed run cannot be overwritten.
"""
from pathlib import Path
import argparse
import base64
import csv
import hashlib
from importlib import metadata
import json
import os
import subprocess
import warnings

os.environ.setdefault('NUMBA_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('MPLBACKEND', 'Agg')
_plot_cache = Path(__file__).resolve().parents[1] / 'data/processed/private/matplotlib-cache'
_plot_cache.mkdir(parents=True, exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(_plot_cache))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.manifold import TSNE, trustworthiness
from threadpoolctl import threadpool_limits
from umap import UMAP
import nbformat

SEED = 20261002
SCOPES = ['combined', 'GSP', 'LR']
COLORS = {'GSP':'#246B8C', 'LR':'#B45B34', 'Validated / Certified':'#287F82',
          'Other elimination':'#A7ADB5', 'Behavioral elimination':'#C45636'}

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def input_path(root, relative): return root/relative.replace('\\','/')
def rows(p): return [s.split() for s in p.read_text().splitlines() if s.strip()]
def report(p):
    rs=rows(p); return [dict(zip(rs[0],r)) for r in rs[1:]]
def dump(p,v): p.write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8')
def save_table(p, headers, records):
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.writer(f,delimiter='\t'); w.writerow(headers); w.writerows(records)

def load_pca_scores(prefix, components=20):
    """Read PLINK components and weight each by its eigenvalue's square root."""
    records=report(prefix.with_suffix('.eigenvec'))
    eigenvalues=np.loadtxt(prefix.with_suffix('.eigenval'))
    assert np.isfinite(eigenvalues).all() and np.all(eigenvalues>=-1e-7)
    eigenvalues=np.maximum(eigenvalues,0)
    if components>len(eigenvalues):
        raise ValueError('Requested more components than the PCA provides')
    vectors=np.array([[float(row[f'PC{k+1}']) for k in range(components)] for row in records])
    scores=vectors*np.sqrt(eigenvalues[:components])
    assert np.isfinite(scores).all()
    return [row['IID'] for row in records], scores, eigenvalues

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    ap.add_argument('--environment',type=Path,required=True)
    ap.add_argument('--output',type=Path,help='Separate run directory inside the repository')
    ap.add_argument('--figures',type=Path,help='Separate figure directory inside the repository')
    args=ap.parse_args(); root=args.root.resolve()
    out=(args.output or root/'data/structure/20261002-pca-umap-tsne').resolve()
    figures=(args.figures or root/'results/private/03_structure/figures').resolve()
    assert out.is_relative_to(root) and figures.is_relative_to(root),'Keep generated outputs inside the study repository'
    private_notebook=root/'results/private/03_structure/03_structure_and_model_design.ipynb'
    if (out/'manifest.private.json').exists(): raise FileExistsError('Completed structure run exists')
    out.mkdir(parents=True,exist_ok=True); figures.mkdir(parents=True,exist_ok=True)
    env=json.loads(args.environment.read_text()); exe=Path(env['executable'])
    assert sha(exe)==env['executable_sha256']
    qm=root/'data/qc/manifest.private.json'; qc=json.loads(qm.read_text())
    assert qc['status']=='complete'
    tracked={}; prefixes={}
    for d in qc['datasets']:
        for p,h in d['files'].items(): assert sha(input_path(root,p))==h; tracked[p]=h
        if d['scenario']=='not-behelim-vs-behelim':
            prefixes[d['breed']]=next(input_path(root,p) for p in d['files'] if p.endswith('.fam')).with_suffix('')
    primary={d['breed']: {r[1] for r in rows(next(input_path(root,p) for p in d['files'] if p.endswith('.fam')))}
             for d in qc['datasets'] if d['scenario']=='success-vs-behelim'}
    commands_path=out/'commands.private.json'
    commands=json.loads(commands_path.read_text()) if commands_path.exists() else []
    def run(label,source,prefix,flags):
        cmd=[str(exe),'--bfile',str(source),'--dog','--keep-allele-order','--threads','1','--memory','2048',*map(str,flags),'--out',str(prefix)]
        old=next((c for c in commands if c['label']==label and c['returncode']==0),None)
        if old:
            assert old['argv']==cmd and prefix.with_suffix('.log').exists()
            return
        proc=subprocess.run(cmd,capture_output=True,text=True)
        prefix.with_suffix('.console.txt').write_text(proc.stdout+proc.stderr,encoding='utf-8')
        commands.append(dict(label=label,argv=cmd,returncode=proc.returncode)); dump(commands_path,commands)
        if proc.returncode: raise RuntimeError(f'PLINK failed at {label}; see private reports')
    variants={b:{r[1]:r for r in rows(p.with_suffix('.bim'))} for b,p in prefixes.items()}
    shared=set(variants['GSP']) & set(variants['LR'])
    assert all(variants['GSP'][s]==variants['LR'][s] for s in shared)
    pruned={}; plotdata={}; summary={}; metrics=[]; sensitivity={}; warnings_record=[]
    for breed in ['GSP','LR']:
        folder=out/breed; folder.mkdir(exist_ok=True)
        prune=folder/'01_ld_pruning'
        run(f'{breed} LD pruning',prefixes[breed],prune,['--indep-pairwise','50','5','0.2'])
        ids=set(prune.with_suffix('.prune.in').read_text().split())
        assert ids and ids<=set(variants[breed]); pruned[breed]=ids
        dest=folder/'02_pruned_genotypes'
        run(f'{breed} pruned genotypes',prefixes[breed],dest,['--extract',prune.with_suffix('.prune.in'),'--make-bed'])
    pooled_ids=shared & pruned['GSP'] & pruned['LR']
    assert len(pooled_ids)>100
    folder=out/'combined'; folder.mkdir(exist_ok=True)
    pooled_list=folder/'common_within_breed_pruned.txt'
    pooled_list.write_text('\n'.join(sorted(pooled_ids))+'\n')
    mergeprefixes=[]
    for breed in ['GSP','LR']:
        dest=folder/f'01_{breed}_shared'
        run(f'{breed} pooled shared genotypes',prefixes[breed],dest,['--extract',pooled_list,'--make-bed'])
        mergeprefixes.append(dest)
    pooled=folder/'02_pruned_genotypes'
    run('combined merge',mergeprefixes[0],pooled,['--bmerge',mergeprefixes[1],'--make-bed'])
    assert len(rows(pooled.with_suffix('.fam')))==164
    label_info={}
    for breed,prefix in prefixes.items():
        for r in rows(prefix.with_suffix('.fam')):
            category='Behavioral elimination' if r[5]=='2' else 'Validated / Certified' if r[1] in primary[breed] else 'Other elimination'
            assert r[1] not in label_info
            label_info[r[1]]=dict(breed=breed,outcome=category,sex=int(r[4]),primary=r[1] in primary[breed])
    for scope in SCOPES:
        folder=out/scope; source=folder/'02_pruned_genotypes'
        fam=rows(source.with_suffix('.fam')); n=len(fam)
        print(f'{scope}: PCA on {n} dogs',flush=True)
        pca=folder/'03_pca'
        run(f'{scope} PCA',source,pca,['--maf','0.025','--pca',str(n-1),'header','tabs','--make-rel','square'])
        ids,scores,eig=load_pca_scores(pca)
        assert set(ids)=={r[1] for r in fam} and len(ids)==n
        assert len(eig)==n-1 and np.isfinite(eig).all() and np.all(eig>=-1e-7)
        eig=np.maximum(eig,0); explained=eig/eig.sum()
        assert scores.shape==(n,20) and np.isfinite(scores).all()
        um=UMAP(n_components=2,n_neighbors=15,min_dist=.1,metric='euclidean',random_state=SEED,n_jobs=1)
        perplexity=min(30,max(5,(n-1)//3))
        with threadpool_limits(limits=1), warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            u=um.fit_transform(scores)
            ts=TSNE(n_components=2,perplexity=perplexity,early_exaggeration=12,learning_rate='auto',max_iter=1500,init='pca',metric='euclidean',random_state=SEED,method='barnes_hut',n_jobs=1)
            t=ts.fit_transform(scores)
            for w in caught: warnings_record.append(dict(scope=scope,message=str(w.message)))
        assert np.isfinite(u).all() and np.isfinite(t).all()
        labels=[label_info[s]['breed' if scope=='combined' else 'outcome'] for s in ids]
        save_table(folder/'coordinates.private.tsv', ['study_id','breed','outcome','sex','in_scheme_a',*[f'PC{k+1}' for k in range(20)],'UMAP1','UMAP2','TSNE1','TSNE2'],
                   [[s,label_info[s]['breed'],label_info[s]['outcome'],label_info[s]['sex'],label_info[s]['primary'],*scores[i].tolist(),*u[i].tolist(),*t[i].tolist()] for i,s in enumerate(ids)])
        plotdata[scope]=dict(ids=ids,labels=labels,PCA=scores[:,:2],UMAP=u,TSNE=t,explained=explained)
        sensitivity[scope]={}
        for method,coords in [('UMAP',u),('TSNE',t)]:
            metrics.append([scope,method,'primary',SEED,15 if method=='UMAP' else perplexity,trustworthiness(scores,coords,n_neighbors=5)])
        print(f'{scope}: checking nonlinear layouts with alternative settings',flush=True)
        with threadpool_limits(limits=1):
            for neighbors in [10,25]:
                coords=UMAP(n_components=2,n_neighbors=neighbors,min_dist=.1,metric='euclidean',random_state=42,n_jobs=1).fit_transform(scores)
                sensitivity[scope,f'UMAP n={neighbors}, seed=42']=coords
                metrics.append([scope,'UMAP','alternative',42,neighbors,trustworthiness(scores,coords,n_neighbors=5)])
            for per in [5, min(30,n-1)]:
                coords=TSNE(n_components=2,perplexity=per,max_iter=1500,early_exaggeration=12,learning_rate='auto',init='pca',metric='euclidean',random_state=42,method='barnes_hut',n_jobs=1).fit_transform(scores)
                sensitivity[scope,f't-SNE p={per}, seed=42']=coords
                metrics.append([scope,'TSNE','alternative',42,per,trustworthiness(scores,coords,n_neighbors=5)])
        if scope!='combined':
            ibd=folder/'04_relatedness'
            run(f'{scope} IBD',source,ibd,['--genome'])
            pairs=report(ibd.with_suffix('.genome'))
            assert len(pairs)==n*(n-1)//2
            pi=np.array([float(d['PI_HAT']) for d in pairs])
            assert np.isfinite(pi).all()
            plotdata[scope]['PI_HAT']=pi
        summary[scope]=dict(dogs=n,markers=len(rows(source.with_suffix('.bim'))),pca_first_two_percent=(explained[:2]*100).tolist(),pca_first_20_percent=float(explained[:20].sum()*100),
                            tsne_perplexity=perplexity,tsne_kl_divergence=float(ts.kl_divergence_),
                            relatedness_pairs_ge_0125=int((pi>=.125).sum()) if scope!='combined' else None,
                            relatedness_pairs_ge_025=int((pi>=.25).sum()) if scope!='combined' else None,
                            relatedness_pairs_ge_05=int((pi>=.5).sum()) if scope!='combined' else None)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white','axes.facecolor':'white'})
    def scatter(ax,scope,method,coords=None,title=None):
        d=plotdata[scope]; xy=d[method] if coords is None else coords
        labels=d['labels']; classes=['GSP','LR'] if scope=='combined' else ['Validated / Certified','Other elimination','Behavioral elimination']
        # Plot behavioral cases last so the smaller case group remains visible.
        for label in classes:
            mask=np.array([s==label for s in labels])
            ax.scatter(xy[mask,0],xy[mask,1],s=32,color=COLORS[label],alpha=.85,edgecolor='white',linewidth=.4,label=f'{label} (n={mask.sum()})')
        ax.set_title(title or ('t-SNE' if method=='TSNE' else method),loc='left',fontweight='bold')
        if method=='PCA':
            ax.set_xlabel(f"PC1 ({d['explained'][0]*100:.1f}% of variation)")
            ax.set_ylabel(f"PC2 ({d['explained'][1]*100:.1f}% of variation)")
        else:
            ax.set_xlabel(('t-SNE' if method=='TSNE' else method)+' dimension 1')
            ax.set_ylabel(('t-SNE' if method=='TSNE' else method)+' dimension 2')
            ax.set_xticks([]); ax.set_yticks([])
        ax.grid(alpha=.12)
    def save(fig,name):
        fig.savefig(figures/f'{name}.png',dpi=180,bbox_inches='tight')
        fig.savefig(figures/f'{name}.svg',bbox_inches='tight')
        plt.close(fig)
    fig,axes=plt.subplots(3,3,figsize=(15,12),layout='constrained')
    for row,scope in enumerate(SCOPES):
        for col,method in enumerate(['PCA','UMAP','TSNE']):
            scatter(axes[row,col],scope,method,title=f'{scope} · '+('t-SNE' if method=='TSNE' else method))
        axes[row,1].legend(loc='upper center',bbox_to_anchor=(.5,-.2),ncol=2 if scope=='combined' else 1,frameon=False,fontsize=9)
    fig.suptitle('Genetic structure of the 2024TSA study dogs',fontweight='bold',fontsize=20)
    save(fig,'structure_overview')
    for scope in SCOPES:
        fig,axs=plt.subplots(1,3,figsize=(15,4.5),layout='constrained')
        for ax,method in zip(axs,['PCA','UMAP','TSNE']): scatter(ax,scope,method)
        axs[1].legend(loc='upper center',bbox_to_anchor=(.5,-.22),frameon=False,fontsize=9,ncol=2 if scope=='combined' else 1)
        fig.suptitle(f'{scope}: {summary[scope]["dogs"]} dogs',fontsize=17,fontweight='bold')
        save(fig,f'{scope}_structure')
        fig,axs=plt.subplots(2,3,figsize=(14,8),layout='constrained')
        panels=[('UMAP',plotdata[scope]['UMAP'],'Primary UMAP (seed=20261002)'),
                *[('UMAP',sensitivity[scope,k],k) for k in [f'UMAP n={n}, seed=42' for n in [10,25]]],
                ('TSNE',plotdata[scope]['TSNE'],'Primary t-SNE (seed=20261002)'),
                *[('TSNE',sensitivity[scope,k],k) for k in [f't-SNE p={p}, seed=42' for p in [5,min(30,len(plotdata[scope]['ids'])-1)]]]]
        for ax,(method,coords,title) in zip(axs.flat,panels): scatter(ax,scope,method,coords=coords,title=title)
        fig.suptitle(f'{scope}: how the layouts change with settings',fontsize=17,fontweight='bold')
        save(fig,f'{scope}_layout_sensitivity')
    fig,axs=plt.subplots(1,3,figsize=(14,4),layout='constrained')
    for ax,scope in zip(axs,SCOPES):
        explained=plotdata[scope]['explained']; ax.plot(np.arange(1,21),explained[:20]*100,'o-',color='#246B8C',ms=4)
        ax.set_title(scope,fontweight='bold'); ax.set_xlabel('Principal component'); ax.set_ylabel('Percent of total variation'); ax.grid(alpha=.15)
    save(fig,'pca_variation')
    fig,axs=plt.subplots(1,2,figsize=(12,4),layout='constrained')
    for ax,scope in zip(axs,['GSP','LR']):
        ax.hist(plotdata[scope]['PI_HAT'],bins=np.linspace(0,1,26),color=COLORS[scope],edgecolor='white')
        ax.set_title(f'{scope}: pairwise genetic sharing',fontweight='bold'); ax.set_xlabel('PLINK estimated fraction shared by descent'); ax.set_ylabel('Number of dog pairs')
    save(fig,'within_breed_relatedness')
    save_table(out/'embedding_neighbor_preservation.tsv',['scope','method','layout','seed','neighborhood_setting','trustworthiness_k5'],metrics)
    versions={p:metadata.version(p) for p in ['numpy','scipy','scikit-learn','umap-learn','numba','matplotlib','nbformat']}
    for p,h in tracked.items(): assert sha(input_path(root,p))==h
    plink_warnings=[]
    for log in out.rglob('*.log'):
        plink_warnings.extend(dict(log=str(log.relative_to(root)),message=line) for line in log.read_text().splitlines() if 'Warning' in line)
    manifest=dict(status='complete',summary=summary,packages=versions,qc_manifest_sha256=sha(qm),script_sha256=sha(Path(__file__)),seed=SEED,
                  pca_components_for_nonlinear=20,pc_scaling='eigenvectors multiplied by sqrt(eigenvalues), no whitening',
                  ld_pruning='50-variant window, step 5, r2 0.2 within breed; pooled intersection of both prune lists and shared QC markers',
                  primary_umap=dict(n_neighbors=15,min_dist=.1,metric='euclidean',n_jobs=1),
                  primary_tsne=dict(perplexity='min(30, (n-1)//3), minimum 5',init='pca',learning_rate='auto',max_iter=1500),
                  python_warnings=warnings_record,plink_warnings=plink_warnings,source_qc_hashes_unchanged=True,
                  figures={str(p.relative_to(root)):sha(p) for p in figures.glob('*')},
                  interpretation='descriptive, no identity inference or outcome-driven exclusion; no model/covariate choice frozen')
    dump(out/'manifest.private.json',manifest)
    dump(out/'summary.private.json',summary)
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__': main()
