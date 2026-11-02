"""Validate paired measured GWAS outputs and create private diagnostic figures."""
from pathlib import Path
import argparse
import csv
import gzip
import json
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import chi2
from notebook_support import read_json,file_digest,study_path

ROOT=Path(__file__).resolve().parents[1]


def table(path):
    with path.open(encoding='utf-8',newline='') as handle:
        return list(csv.DictReader(handle,delimiter='\t'))


def number(value):
    try:return float(value)
    except (TypeError,ValueError):return float('nan')


def validated_job(folder):
    receipt=read_json(folder/'receipt.private.json')
    assert all(file_digest(folder/p)==h for p,h in receipt['files'].items()),'Completed output changed'
    markers=table(folder/'markers.private.tsv');meta=read_json(folder/'input.private.json')
    diagnostic=read_json(folder/'diagnostics.private.json')
    linear=table(folder/'output/gemma.assoc.txt')
    logistic=table(folder/'gmmat.score.txt') if diagnostic['gmmat_null_converged'] else None
    assert len(markers)==len(linear)==meta['markers']
    if logistic is not None:assert len(logistic)==len(markers)
    output=[]
    for i,(m,l) in enumerate(zip(markers,linear)):
        assert m['marker']==l['rs'] and m['ALT']==l['allele1'] and m['REF']==l['allele0']
        assert abs(number(m['AF'])-number(l['af']))<.0006
        g=logistic[i] if logistic is not None else {}
        if g:
            assert m['marker']==g['SNP'] and int(g['N'])==meta['dogs']
            assert abs(number(m['AF'])-number(g['AF']))<1e-6
        row=dict(m,GMMAT_p=number(g.get('PVAL')),GMMAT_score=number(g.get('SCORE')),
            GMMAT_variance=number(g.get('VAR')),GEMMA_p=number(l['p_score']),
            GEMMA_beta=number(l['beta']),GEMMA_se=number(l['se']),
            GEMMA_wald_p=number(l['p_wald']),GEMMA_lrt_p=number(l['p_lrt']),
            GEMMA_lambda_mle=number(l['l_mle']),model_dogs=meta['dogs'],cases=meta['cases'])
        output.append(row)
    return output,diagnostic


def summarize(run,output):
    manifest=read_json(run/'manifest.private.json');plan=read_json(run/'plan.private.json')
    assert manifest['status']=='paired_measured_scans_complete'
    assert all(file_digest(study_path(ROOT,p))==h for p,h in plan['sources'].items())
    output.mkdir(parents=True,exist_ok=True)
    reporting=read_json(ROOT/'config/gwas_reporting.json')
    thresholds={breed:reporting['family_alpha']/count for breed,count in reporting['candidate_markers'].items()}
    summaries=[];diagnostics=[]
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    for group in manifest['groups']:
        breed=group['breed'];scenario=group['scenario'];combined=[]
        threshold=thresholds[breed]
        assert reporting['candidate_markers'][breed]==group['candidate_markers']
        for job in manifest['jobs']:
            if job['breed']!=breed or job['scenario']!=scenario:continue
            records,d=validated_job(study_path(ROOT,job['folder']));combined.extend(records)
            diagnostics.append(dict(job,**d))
        combined.sort(key=lambda r:(int(r['chromosome']),int(r['position'])))
        if manifest['chromosomes']==list(range(1,39)):
            assert len(combined)==group['eligible_markers']
        stem=breed+'_'+scenario
        with gzip.open(output/(stem+'.private.tsv.gz'),'wt',encoding='utf-8',newline='') as h:
            writer=csv.DictWriter(h,fieldnames=list(combined[0]),delimiter='\t');writer.writeheader();writer.writerows(combined)
        chrom=np.array([int(r['chromosome']) for r in combined]);position=np.array([int(r['position']) for r in combined])
        offsets={};cursor=0
        for c in sorted(set(chrom)):
            offsets[c]=cursor;cursor+=position[chrom==c].max()+1_000_000
        x=np.array([offsets[c]+p for c,p in zip(chrom,position)])/1e6
        midpoints=[(x[chrom==c].min()+x[chrom==c].max())/2 for c in offsets]
        fig,axes=plt.subplots(4,1,figsize=(14,14),sharex=True)
        qq,qa=plt.subplots(2,2,figsize=(11,9))
        title=f'{breed}: '+('Scheme A, primary' if scenario==manifest['config']['primary_cohort'] else 'Scheme B, sensitivity')
        tests=[('GEMMA','Wald','GEMMA_wald_p'),('GEMMA','Likelihood ratio','GEMMA_lrt_p'),
               ('GEMMA','Score','GEMMA_p'),('GMMAT','Score','GMMAT_p')]
        for ax,qax,(method,test,key) in zip(axes,qa.flat,tests):
            p=np.array([r[key] for r in combined]);valid=np.isfinite(p)&(p>=0)&(p<=1)
            positive=p[valid];logp=-np.log10(np.maximum(p,1e-300))
            for c in offsets:
                mask=(chrom==c)&valid
                ax.scatter(x[mask],logp[mask],s=2,color='#246b8f' if c%2 else '#72a5ba',rasterized=True)
            ax.axhline(-math.log10(threshold),color='#b05c36',ls='--',lw=1,label=f'{breed} within-analysis cutoff: p = {threshold:.3g}')
            ax.legend(loc='upper right',fontsize=8)
            ax.set(ylabel='−log₁₀(p)',title=method+' — '+test+(' (method comparison)' if method=='GMMAT' else ''));ax.grid(axis='y',alpha=.12)
            ordered=np.sort(positive);expected=-np.log10((np.arange(len(ordered))+.5)/len(ordered))
            observed=-np.log10(np.maximum(ordered,1e-300));take=np.unique(np.r_[np.arange(min(100,len(ordered))),np.linspace(0,len(ordered)-1,min(1500,len(ordered))).astype(int)])
            qax.scatter(expected[take],observed[take],s=8,color='#246b8f');limit=max(expected.max(),observed.max())
            qax.plot([0,limit],[0,limit],color='.5',ls='--');qax.set(title=method+' — '+test,xlabel='Expected −log₁₀(p)',ylabel='Observed −log₁₀(p)')
            inflation=float(np.median(chi2.isf(positive,1))/chi2.ppf(.5,1))
            ranked=sorted([r for r in combined if math.isfinite(r[key]) and 0<=r[key]<=1],key=lambda r:r[key])
            hits=[r for r in ranked if r[key]<threshold]
            with (output/(stem+'_'+method+'_'+test.replace(' ','_')+'_candidates.private.tsv')).open('w',encoding='utf-8',newline='') as h:
                writer=csv.DictWriter(h,fieldnames=list(combined[0]),delimiter='\t');writer.writeheader();writer.writerows(hits)
            summaries.append(dict(breed=breed,scenario=scenario,method=method,test=test,p_value_column=key,
                reporting_role='laboratory workflow' if method=='GEMMA' else 'method comparison',dogs=group['dogs'],cases=group['cases'],
                candidates=group['candidate_markers'],eligible_markers=group['eligible_markers'],scanned_markers=len(combined),
                valid_tests=int(valid.sum()),invalid_tests=int((~valid).sum()),minimum_p=float(positive.min()),
                median_chi_square_ratio=inflation,threshold=threshold,below_analysis_threshold=len(hits),
                top_variants=ranked[:10],gemma_mle_boundary_count=sum(r['GEMMA_lambda_mle']>=99999 or r['GEMMA_lambda_mle']<=1.00001e-5 for r in combined) if method=='GEMMA' else None))
            if method=='GEMMA':
                summaries[-1].update(wald_invalid=sum(not math.isfinite(r['GEMMA_wald_p']) for r in combined),
                    lrt_invalid=sum(not math.isfinite(r['GEMMA_lrt_p']) for r in combined),
                    effect_invalid=sum(not math.isfinite(r['GEMMA_beta']) or not math.isfinite(r['GEMMA_se']) for r in combined))
        for ax in axes:
            ax.set_xticks(midpoints, labels=list(offsets))
            ax.tick_params(axis='x', labelbottom=True, labelsize=8)
        axes[-1].set_xlabel('CanFam4 autosome')
        axes[-1].tick_params(axis='x',labelsize=8,labelrotation=90)
        fig.suptitle(title+' — measured-genotype GWAS');fig.tight_layout();fig.savefig(output/(stem+'_manhattan.png'),dpi=160);plt.close(fig)
        qq.suptitle(title+' — GEMMA tests and GMMAT comparison');qq.tight_layout();qq.savefig(output/(stem+'_qq.png'),dpi=160);plt.close(qq)
    summary=dict(status='measured_gwas_summarized',full_autosomes=manifest['chromosomes']==list(range(1,39)),
        thresholds=thresholds,reporting_policy=reporting,reporting_policy_sha256=file_digest(ROOT/'config/gwas_reporting.json'),
        groups=manifest['groups'],tests=summaries,
        chromosome_null_models=len(diagnostics),gmmat_null_models_converged=sum(d['gmmat_null_converged'] for d in diagnostics),
        discoveries_validated=False,imputed_variants_tested=False)
    (output/'summary.private.json').write_text(json.dumps(summary,indent=2,allow_nan=True)+'\n',encoding='utf-8')
    (output/'diagnostics.private.json').write_text(json.dumps(diagnostics,indent=2)+'\n',encoding='utf-8')
    print('Measured scan summaries:',[(r['breed'],r['scenario'],r['method'],r['test'],r['below_analysis_threshold']) for r in summaries])
    return summary


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir',type=Path,default=ROOT/'data/07-gwas/measured-20261004')
    parser.add_argument('--results-dir',type=Path,default=ROOT/'results/private/version2/07_breed_specific_gwas')
    args=parser.parse_args();summarize(args.run_dir,args.results_dir)
