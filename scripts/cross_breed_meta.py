"""Separate GEMMA effect meta-analysis, GMMAT score meta-analysis and Fisher comparisons."""
from pathlib import Path
import argparse,csv,gzip,json
import numpy as np
from scipy.stats import norm,chi2
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from notebook_support import read_json,file_digest,study_path
from replication_models import rows
from summarize_gwas import validated_job

ROOT=Path(__file__).resolve().parents[1]
TESTS=[('GEMMA inverse variance','iv_p'),('GMMAT score meta','score_meta_p'),
       ('GEMMA Wald Fisher','GEMMA_wald_p_fisher'),('GEMMA likelihood ratio Fisher','GEMMA_lrt_p_fisher'),
       ('GEMMA score Fisher','GEMMA_p_fisher'),('GMMAT score Fisher','GMMAT_p_fisher')]

def write_table(path,records):
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'wt',encoding='utf-8',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(records[0]),delimiter='\t');w.writeheader();w.writerows(records)

def aligned_family(scenario):
    maps=[]
    for breed in ['LR','GSP']:
        fam=rows(ROOT/f'data/02-aligned/{scenario}/{breed}/2024TSA_{breed}.fam')
        bim=rows(ROOT/f'data/02-aligned/{scenario}/{breed}/2024TSA_{breed}.bim')
        assert len({r[1] for r in bim})==len(bim),'Duplicate marker IDs'
        maps.append({r[1]:(r[0],r[3],r[4],r[5]) for r in bim})
        if breed=='LR':ids={r[1] for r in fam}
        else:assert not ids.intersection(r[1] for r in fam),'Overlapping dog IDs'
    shared=set(maps[0])&set(maps[1]);mismatch=[s for s in shared if maps[0][s]!=maps[1][s]]
    assert not mismatch,'Shared marker coordinate/allele mismatch: '+str(mismatch[:10])
    return shared

def fixed_effect(beta,se):
    weights=1/np.square(se);total=weights.sum(axis=0)
    mean=(weights*beta).sum(axis=0)/total;error=np.sqrt(1/total)
    q=(weights*np.square(beta-mean)).sum(axis=0)
    return mean,error,2*norm.sf(np.abs(mean/error)),q,weights[0]/total

def figures(records,threshold,stem,output):
    c=np.array([int(r['chromosome']) for r in records]);p=np.array([int(r['position']) for r in records])
    offsets={};cursor=0
    for ch in sorted(set(c)):offsets[ch]=cursor;cursor+=p[c==ch].max()+1000000
    x=np.array([offsets[ch]+pos for ch,pos in zip(c,p)])/1e6
    fig,axes=plt.subplots(6,1,figsize=(14,18),sharex=True);qq,qaxes=plt.subplots(3,2,figsize=(11,13))
    for ax,qax,(label,key) in zip(axes,qaxes.flat,TESTS):
        prob=np.array([r[key] for r in records]);ok=np.isfinite(prob)&(prob>=0)&(prob<=1)
        for ch in offsets:
            mask=(c==ch)&ok;ax.scatter(x[mask],-np.log10(np.maximum(prob[mask],1e-300)),s=2,color='#246b8f' if ch%2 else '#72a5ba',rasterized=True)
        ax.axhline(-np.log10(threshold),ls='--',color='#b05c36');ax.set(title=label,ylabel='−log₁₀(p)')
        observed=np.sort(prob[ok]);expected=-np.log10((np.arange(len(observed))+.5)/len(observed))
        take=np.unique(np.r_[np.arange(min(100,len(observed))),np.linspace(0,len(observed)-1,min(1000,len(observed))).astype(int)])
        y=-np.log10(np.maximum(observed[take],1e-300));qax.scatter(expected[take],y,s=5)
        limit=max(expected.max(),y.max());qax.plot([0,limit],[0,limit],ls='--',color='.5');qax.set(title=label,xlabel='Expected −log₁₀(p)',ylabel='Observed −log₁₀(p)')
    for ax in axes:
        ax.set(xticks=[(x[c==ch].min()+x[c==ch].max())/2 for ch in offsets], xticklabels=list(offsets))
        ax.tick_params(axis='x', labelbottom=True, labelsize=8)
    axes[-1].set_xlabel('CanFam4 chromosome')
    fig.suptitle(stem+' — separate cross-breed analyses');fig.tight_layout();fig.savefig(output/(stem+'_manhattan.png'),dpi=130);plt.close(fig)
    qq.tight_layout();qq.savefig(output/(stem+'_qq.png'),dpi=130);plt.close(qq)
    best=sorted([r for r in records if np.isfinite(r['iv_p'])],key=lambda r:r['iv_p'])[:10]
    forest,ax=plt.subplots(figsize=(10,7));ys=np.arange(len(best))*3
    for offset,breed,col in [(-.6,'LR','#b05c36'),(0,'GSP','#246b8f'),(.6,'meta','#444444')]:
        means=[r['iv_beta' if breed=='meta' else breed+'_GEMMA_beta'] for r in best]
        errors=[1.96*r['iv_se' if breed=='meta' else breed+'_GEMMA_se'] for r in best]
        ax.errorbar(means,ys+offset,xerr=errors,fmt='o',label=breed,color=col,markersize=4)
    ax.axvline(0,color='.6',ls='--');ax.set(yticks=ys,yticklabels=[r['marker'] for r in best],xlabel='ALT effect on binary outcome scale (approximate 95% intervals)',title='Top inverse-variance markers — ranking for display, not confirmation');ax.legend();forest.tight_layout();forest.savefig(output/(stem+'_effects.png'),dpi=140);plt.close(forest)

def run(run_dir,output):
    cfg=read_json(ROOT/'config/cross_breed.json');plan=read_json(run_dir/'plan.private.json');manifest=read_json(run_dir/'manifest.private.json')
    assert manifest['status']=='paired_measured_scans_complete'
    assert all(file_digest(study_path(ROOT,p))==h for p,h in plan['sources'].items())
    output.mkdir(parents=True,exist_ok=True);summary=[];sources={};family_counts=[]
    for scenario in cfg['schemes']:
        common=aligned_family(scenario);maps=[]
        for breed in ['LR','GSP']:
            data=[]
            for job in manifest['jobs']:
                if job['breed']==breed and job['scenario']==scenario:
                    folder=study_path(ROOT,job['folder']);data.extend(validated_job(folder)[0]);sources[str((folder/'receipt.private.json').relative_to(ROOT))]=file_digest(folder/'receipt.private.json')
            maps.append({r['marker']:r for r in data})
        keys=sorted(common&set(maps[0])&set(maps[1]),key=lambda k:(int(maps[0][k]['chromosome']),int(maps[0][k]['position']),k))
        records=[]
        for key in keys:
            a,b=maps[0][key],maps[1][key]
            assert all(a[f]==b[f] for f in ['chromosome','position','ALT','REF'])
            row={f:a[f] for f in ['marker','chromosome','position','ALT','REF']}
            for breed,r in zip(['LR','GSP'],[a,b]):
                row.update({breed+'_'+f:r[f] for f in ['GEMMA_beta','GEMMA_se','GEMMA_wald_p','GEMMA_lrt_p','GEMMA_p','GMMAT_p','GMMAT_score','GMMAT_variance','AF','MAF','model_dogs','cases','called_dogs']})
            records.append(row)
        beta=np.array([[r[breed+'_GEMMA_beta'] for r in records] for breed in ['LR','GSP']]);se=np.array([[r[breed+'_GEMMA_se'] for r in records] for breed in ['LR','GSP']])
        ok=np.isfinite(beta).all(axis=0)&np.isfinite(se).all(axis=0)&(se>0).all(axis=0)
        with np.errstate(invalid='ignore',divide='ignore',over='ignore'):mean,error,p,q,fraction=fixed_effect(beta,se)
        u=np.array([[r[breed+'_GMMAT_score'] for r in records] for breed in ['LR','GSP']]);v=np.array([[r[breed+'_GMMAT_variance'] for r in records] for breed in ['LR','GSP']])
        score_ok=np.isfinite(u).all(axis=0)&np.isfinite(v).all(axis=0)&(v>0).all(axis=0)
        score_p=chi2.sf(np.square(u.sum(axis=0))/v.sum(axis=0),1)
        for i,r in enumerate(records):
            r.update(iv_beta=float(mean[i]) if ok[i] else np.nan,iv_se=float(error[i]) if ok[i] else np.nan,iv_p=float(p[i]) if ok[i] else np.nan,
                Q=float(q[i]) if ok[i] else np.nan,Q_p=float(chi2.sf(q[i],1)) if ok[i] else np.nan,
                I2=float(max(0,1-1/q[i])) if ok[i] and q[i]>0 else (0.0 if ok[i] else np.nan),LR_weight_fraction=float(fraction[i]) if ok[i] else np.nan,
                same_effect_direction=bool(beta[0,i]*beta[1,i]>0) if ok[i] else None,
                score_meta_score=float(u[:,i].sum()) if score_ok[i] else np.nan,score_meta_variance=float(v[:,i].sum()) if score_ok[i] else np.nan,
                score_meta_p=float(score_p[i]) if score_ok[i] else np.nan)
            for key in cfg['fisher_tests']:
                ps=np.array([r['LR_'+key],r['GSP_'+key]])
                r[key+'_fisher']=float(chi2.sf(-2*np.log(ps).sum(),4)) if np.isfinite(ps).all() and ((ps>0)&(ps<=1)).all() else np.nan
        threshold=cfg['alpha']/len(common);family_counts.append(dict(scenario=scenario,common_candidates=len(common),eligible_both=len(keys),threshold=threshold,invalid_effect_pairs=int((~ok).sum())))
        write_table(output/(scenario+'.private.tsv.gz'),records)
        # Native GMMAT validation inputs retain consistent ALT allele orientation.
        for breed in ['LR','GSP']:
            native=[dict(SNP=r['marker'],A1=r['ALT'],A2=r['REF'],N=r[breed+'_model_dogs'],AF=r[breed+'_AF'],SCORE=r[breed+'_GMMAT_score'],VAR=r[breed+'_GMMAT_variance'],PVAL=r[breed+'_GMMAT_p']) for r in records]
            write_table(output/(scenario+'_'+breed+'_score.private.tsv'),native)
        for label,key in TESTS:
            valid=[r for r in records if np.isfinite(r[key]) and 0<=r[key]<=1];ranked=sorted(valid,key=lambda r:r[key]);hits=[r for r in ranked if r[key]<threshold]
            write_table(output/(scenario+'_'+key+'_candidates.private.tsv'),hits if hits else [dict.fromkeys(records[0],'')])
            summary.append(dict(scenario=scenario,analysis=label,p_column=key,valid_tests=len(valid),invalid_tests=len(records)-len(valid),minimum_p=ranked[0][key],below_cutoff=len(hits),threshold=threshold,top=ranked[:10]))
        figures(records,threshold,scenario,output)
        print('Cross-breed meta:',family_counts[-1],flush=True)
    summary=dict(status='cross_breed_meta_complete',config=cfg,families=family_counts,analyses=summary,input_receipts=sources,config_sha256=file_digest(ROOT/'config/cross_breed.json'),script_sha256=file_digest(Path(__file__)),original_sources_unchanged=True)
    (output/'summary.private.json').write_text(json.dumps(summary,indent=2,allow_nan=True)+'\n',encoding='utf-8')
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run-dir',type=Path,default=ROOT/'data/07-gwas/measured-20261004');p.add_argument('--results-dir',type=Path,default=ROOT/'results/private/version2/08_cross_breed')
    args=p.parse_args();run(args.run_dir,args.results_dir)
