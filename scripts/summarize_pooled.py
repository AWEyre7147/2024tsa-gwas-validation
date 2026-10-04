"""Validate exploratory pooled fits and save clearly separated diagnostic figures."""
from pathlib import Path
import argparse,json
import numpy as np
from scipy.stats import chi2
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from notebook_support import read_json,file_digest,study_path
from summarize_gwas import validated_job
from cross_breed_meta import write_table

ROOT=Path(__file__).resolve().parents[1]

def run(run_dir,output):
    manifest=read_json(run_dir/'manifest.private.json');plan=read_json(run_dir/'plan.private.json');assert manifest['status']=='exploratory_pooled_scans_complete'
    assert all(file_digest(study_path(ROOT,p))==h for p,h in plan['sources'].items())
    output.mkdir(parents=True,exist_ok=True);analyses=[];diagnostics=[]
    for group in manifest['groups']:
        scenario=group['scenario'];records=[]
        for job in manifest['jobs']:
            if job['scenario']!=scenario:continue
            folder=study_path(ROOT,job['folder']);data,d=validated_job(folder);records.extend(data);info=read_json(folder/'input.private.json')
            diagnostics.append(dict(job,**d,minimum_eigenvalue=info['minimum_eigenvalue'],cross_breed_mean=info['cross_breed_mean']))
        assert len(records)==group['eligible_markers'];records.sort(key=lambda r:(int(r['chromosome']),int(r['position'])))
        write_table(output/(scenario+'.private.tsv.gz'),records);threshold=plan['config']['alpha']/group['candidate_markers']
        c=np.array([int(r['chromosome']) for r in records]);pos=np.array([int(r['position']) for r in records]);offsets={};cursor=0
        for ch in sorted(set(c)):offsets[ch]=cursor;cursor+=pos[c==ch].max()+1000000
        x=np.array([offsets[ch]+p for ch,p in zip(c,pos)])/1e6
        fig,axes=plt.subplots(4,1,figsize=(14,13),sharex=True);qq,qaxes=plt.subplots(2,2,figsize=(11,9))
        for ax,qax,(label,key) in zip(axes,qaxes.flat,[('GEMMA Wald','GEMMA_wald_p'),('GEMMA likelihood ratio','GEMMA_lrt_p'),('GEMMA score','GEMMA_p'),('GMMAT score','GMMAT_p')]):
            p=np.array([r[key] for r in records]);ok=np.isfinite(p)&(p>=0)&(p<=1);ranked=sorted([r for r in records if np.isfinite(r[key]) and 0<=r[key]<=1],key=lambda r:r[key]);hits=[r for r in ranked if r[key]<threshold]
            for ch in offsets:
                mask=(c==ch)&ok;ax.scatter(x[mask],-np.log10(np.maximum(p[mask],1e-300)),s=2,color='#246b8f' if ch%2 else '#72a5ba',rasterized=True)
            ax.axhline(-np.log10(threshold),ls='--',color='#b05c36');ax.set(title=label+' — exploratory pooled sensitivity',ylabel='−log₁₀(p)')
            ordered=np.sort(p[ok]);expected=-np.log10((np.arange(len(ordered))+.5)/len(ordered));take=np.unique(np.r_[np.arange(min(100,len(ordered))),np.linspace(0,len(ordered)-1,min(1000,len(ordered))).astype(int)])
            observed=-np.log10(np.maximum(ordered[take],1e-300));qax.scatter(expected[take],observed,s=5);limit=max(expected.max(),observed.max());qax.plot([0,limit],[0,limit],ls='--',color='.5');qax.set(title=label,xlabel='Expected −log₁₀(p)',ylabel='Observed −log₁₀(p)')
            # Empty candidate files contain a header only.
            import csv
            with (output/(scenario+'_'+key+'_candidates.private.tsv')).open('w',encoding='utf-8',newline='') as h:
                w=csv.DictWriter(h,fieldnames=list(records[0]),delimiter='\t');w.writeheader();w.writerows(hits)
            analyses.append(dict(scenario=scenario,analysis=label,valid_tests=int(ok.sum()),invalid_tests=int((~ok).sum()),minimum_p=ranked[0][key],below_cutoff=len(hits),threshold=threshold,median_chi_square_ratio=float(np.median(chi2.isf(ordered,1))/chi2.ppf(.5,1)),gemma_mle_boundary_count=sum(r['GEMMA_lambda_mle']>=99999 or r['GEMMA_lambda_mle']<=1.00001e-5 for r in records) if label.startswith('GEMMA') else None,top=ranked[:10]))
        axes[-1].set(xticks=[(x[c==ch].min()+x[c==ch].max())/2 for ch in offsets],xticklabels=list(offsets),xlabel='CanFam4 chromosome')
        fig.suptitle(scenario+' — pooled measured dogs, exploratory');fig.tight_layout();fig.savefig(output/(scenario+'_manhattan.png'),dpi=130);plt.close(fig);qq.tight_layout();qq.savefig(output/(scenario+'_qq.png'),dpi=130);plt.close(qq)
    summary=dict(status='exploratory_pooled_review_complete',groups=manifest['groups'],analyses=analyses,diagnostics=diagnostics,original_sources_unchanged=True,publication_inclusion_undecided=True,plan_sha256=file_digest(run_dir/'plan.private.json'),config_sha256=file_digest(ROOT/'config/cross_breed.json'),script_sha256=file_digest(Path(__file__)))
    (output/'summary.private.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print('Exploratory pooled summaries complete.',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run-dir',type=Path,default=ROOT/'data/gwas/pooled-measured-20261004');p.add_argument('--results-dir',type=Path,default=ROOT/'results/private/version2/08_cross_breed/pooled');a=p.parse_args();run(a.run_dir,a.results_dir)
