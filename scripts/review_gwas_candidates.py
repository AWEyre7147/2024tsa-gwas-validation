"""Review measured candidate regions, genotypes, LD and fit diagnostics without new GWAS tests."""
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
from notebook_support import read_json,file_digest,study_path
from breed_gwas import decode

ROOT=Path(__file__).resolve().parents[1]


def read_table(path):
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt',encoding='utf-8',newline='') as h:return list(csv.DictReader(h,delimiter='\t'))


def write_table(path,records):
    with path.open('w',encoding='utf-8',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(records[0]),delimiter='\t');w.writeheader();w.writerows(records)


def r_squared(a,b):
    ok=np.isfinite(a)&np.isfinite(b)
    if ok.sum()<3 or np.std(a[ok])==0 or np.std(b[ok])==0:return None,int(ok.sum())
    return float(np.corrcoef(a[ok],b[ok])[0,1]**2),int(ok.sum())


def pair_summary(K,selected):
    sub=K[np.ix_(selected,selected)];values=sub[np.triu_indices(len(selected),1)]
    return dict(pairs=len(values),mean=float(np.mean(values)) if len(values) else None,
        median=float(np.median(values)) if len(values) else None,maximum=float(np.max(values)) if len(values) else None)


def review(run,scan_results,output):
    cfg=read_json(ROOT/'config/candidate_review.json');policy=read_json(ROOT/'config/gwas_reporting.json')
    manifest=read_json(run/'manifest.private.json');plan=read_json(run/'plan.private.json')
    scan_summary=read_json(scan_results/'summary.private.json')
    assert scan_summary['reporting_policy_sha256']==file_digest(ROOT/'config/gwas_reporting.json')
    assert all(file_digest(study_path(ROOT,p))==h for p,h in plan['sources'].items())
    output.mkdir(parents=True,exist_ok=True);regions=[];flow=[];all_ld=[];all_counts=[];all_models=[]
    for breed in ['LR','GSP']:
        primary=read_table(scan_results/(breed+'_success-vs-behelim.private.tsv.gz'))
        sensitivity={r['marker']:r for r in read_table(scan_results/(breed+'_not-behelim-vs-behelim.private.tsv.gz'))}
        threshold=policy['family_alpha']/policy['candidate_markers'][breed]
        candidates=[r for r in primary if any(math.isfinite(float(r[k])) and float(r[k])<threshold for k in ['GEMMA_wald_p','GEMMA_lrt_p'])]
        candidates.sort(key=lambda r:(int(r['chromosome']),int(r['position'])))
        windows=[]
        for row in candidates:
            c=int(row['chromosome']);p=int(row['position']);start=max(1,p-cfg['window_bp']);end=p+cfg['window_bp']
            if windows and windows[-1]['chromosome']==c and start<=windows[-1]['end']:
                windows[-1]['end']=max(end,windows[-1]['end']);windows[-1]['selected'].append(row)
            else:windows.append(dict(chromosome=c,start=start,end=end,selected=[row]))
        prefix=ROOT/'data/aligned/success-vs-behelim'/breed/('2024TSA_'+breed)
        fam,bim,g=decode(prefix);marker_index={r[1]:i for i,r in enumerate(bim)}
        case=np.array([r[5]=='2' for r in fam]);male=np.array([r[4]=='1' for r in fam])
        assert len({r[1] for r in fam})==len(fam)
        flow.append(dict(breed=breed,selected_markers=len(candidates),display_regions=len(windows),primary_dogs=len(fam),cases=int(case.sum())))
        for number,window in enumerate(windows,1):
            region_id=f'{breed}_{number:02d}_chr{window["chromosome"]}'
            selected=window['selected'];anchor=min(selected,key=lambda r:float(r['GEMMA_wald_p']))
            nearby=[r for r in primary if int(r['chromosome'])==window['chromosome'] and window['start']<=int(r['position'])<=window['end']]
            a=g[marker_index[anchor['marker']]]
            assert int(np.isfinite(a[case]).sum())==int(anchor['case_called'])
            assert int(np.nansum(a[case]))==int(anchor['case_alt_copies'])
            ld=[]
            for row in nearby:
                value,n=r_squared(a,g[marker_index[row['marker']]])
                ld.append(dict(region=region_id,marker=row['marker'],chromosome=window['chromosome'],position=int(row['position']),r_squared=value,paired_dogs=n,selected_candidate=row['marker'] in {s['marker'] for s in selected}))
            all_ld.extend(ld)
            af=float(np.nanmean(a)/2);minor_alt=af<=.5;minor=a if minor_alt else 2-a;carriers=np.flatnonzero(np.isfinite(minor)&(minor>0))
            folder=run/breed/'success-vs-behelim'/f'chr{window["chromosome"]:02d}'
            K=np.loadtxt(folder/'kinship.txt');dogs=read_table(folder/'dogs.private.tsv')
            assert [d['id'] for d in dogs]==[r[1] for r in fam]
            counts=[]
            for mask,label in [(case,'case'),(~case,'control')]:
                for sex_mask,sex in [(np.ones(len(fam),bool),'all'),(male,'male'),(~male,'female')]:
                    values=a[mask&sex_mask];m=minor[mask&sex_mask]
                    counts.append(dict(region=region_id,marker=anchor['marker'],group=label,sex=sex,dogs=len(values),called=int(np.isfinite(values).sum()),missing=int(np.isnan(values).sum()),ALT_0=int((values==0).sum()),ALT_1=int((values==1).sum()),ALT_2=int((values==2).sum()),minor_allele_carriers=int((m>0).sum())))
            all_counts.extend(counts)
            models=[]
            for candidate in selected:
                for scenario,r in [('A',candidate),('B',sensitivity.get(candidate['marker']))]:
                    if r is None:continue
                    native_folder=run/breed/('success-vs-behelim' if scenario=='A' else 'not-behelim-vs-behelim')/f'chr{window["chromosome"]:02d}'
                    receipt=read_json(native_folder/'receipt.private.json')
                    native_path=native_folder/'output/gemma.assoc.txt'
                    assert file_digest(native_path)==receipt['files']['output/gemma.assoc.txt']
                    native=next(row for row in read_table(native_path) if row['rs']==r['marker'])
                    remle=float(native['l_remle'])
                    models.append(dict(region=region_id,marker=r['marker'],scheme=scenario,
                        GEMMA_beta=float(r['GEMMA_beta']),GEMMA_SE=float(r['GEMMA_se']),GEMMA_wald_p=float(r['GEMMA_wald_p']),
                        GEMMA_lrt_p=float(r['GEMMA_lrt_p']),GEMMA_score_p=float(r['GEMMA_p']),GMMAT_score_p=float(r['GMMAT_p']),
                        GEMMA_lambda_mle=float(r['GEMMA_lambda_mle']),GEMMA_mle_boundary=float(r['GEMMA_lambda_mle'])>=99999 or float(r['GEMMA_lambda_mle'])<=1.00001e-5,
                        GEMMA_lambda_remle=remle,GEMMA_remle_boundary=remle>=99999 or remle<=1.00001e-5,
                        observed_dogs=int(r['called_dogs']),model_dogs=int(r['model_dogs'])))
            all_models.extend(models)
            ctrl_af=float(np.nanmean(a[~case])/2)
            leave_one=[]
            for index in np.flatnonzero(case):
                remaining=case.copy();remaining[index]=False
                leave_one.append(float(np.nanmean(a[remaining])/2-ctrl_af))
            entry=dict(region=region_id,breed=breed,chromosome=window['chromosome'],start=window['start'],end=window['end'],
                selected_markers=[r['marker'] for r in selected],anchor=anchor['marker'],anchor_position=int(anchor['position']),
                ALT=anchor['ALT'],REF=anchor['REF'],carrier_allele=anchor['ALT'] if minor_alt else anchor['REF'],threshold=threshold,
                measured_neighbors=len(nearby),nonanchor_neighbors_r2_ge_08=sum(r['marker']!=anchor['marker'] and r['r_squared'] is not None and r['r_squared']>=.8 for r in ld),
                case_carriers=int(np.sum((minor>0)&case)),control_carriers=int(np.sum((minor>0)&~case)),
                carrier_pair_relationship=pair_summary(K,carriers),all_pair_relationship=pair_summary(K,np.arange(len(fam))),
                case_carrier_pair_relationship=pair_summary(K,np.flatnonzero((minor>0)&case)),
                case_control_ALT_frequency_difference=float(np.nanmean(a[case])/2-ctrl_af),
                leave_one_case_out_frequency_difference_min=min(leave_one),leave_one_case_out_frequency_difference_max=max(leave_one),
                independence_established=False,causal_variant_identified=False)
            regions.append(entry)
            fig,axes=plt.subplots(2,2,figsize=(12,8));positions=np.array([int(r['position']) for r in nearby])/1e6
            color=np.array([r['r_squared'] if r['r_squared'] is not None else np.nan for r in ld])
            for ax,key,title in [(axes[0,0],'GEMMA_wald_p','GEMMA Wald'),(axes[0,1],'GEMMA_lrt_p','GEMMA likelihood ratio')]:
                p=np.array([float(r[key]) for r in nearby]);valid=np.isfinite(p)&(p>=0)&(p<=1)
                scatter=ax.scatter(positions[valid],-np.log10(np.maximum(p[valid],1e-300)),c=color[valid],vmin=0,vmax=1,cmap='viridis',s=18,plotnonfinite=True)
                ax.axhline(-math.log10(threshold),color='#b05c36',ls='--');ax.axvline(entry['anchor_position']/1e6,color='.5',ls=':')
                ax.set(title=title,xlabel=f'CanFam4 chr{window["chromosome"]} position (Mb)',ylabel='−log₁₀(p)')
            for mask,label,offset,color_bar in [(case,'Cases',-.18,'#b05c36'),(~case,'Controls',.18,'#246b8f')]:
                values=a[mask];fractions=[np.sum(values==k)/np.isfinite(values).sum() for k in range(3)]
                axes[1,0].bar(np.arange(3)+offset,fractions,width=.35,label=label,color=color_bar)
            axes[1,0].set(xticks=np.arange(3),xticklabels=[f'{anchor["REF"]}/{anchor["REF"]}',f'{anchor["REF"]}/{anchor["ALT"]}',f'{anchor["ALT"]}/{anchor["ALT"]}'],ylabel='Fraction of observed genotypes',title='Display-anchor genotypes');axes[1,0].legend()
            all_pairs=K[np.triu_indices(len(K),1)];cp=K[np.ix_(carriers,carriers)];carrier_pairs=cp[np.triu_indices(len(carriers),1)]
            axes[1,1].hist(all_pairs,bins=25,density=True,alpha=.5,label='All dog pairs',color='#246b8f')
            if len(carrier_pairs):axes[1,1].hist(carrier_pairs,bins=15,density=True,alpha=.5,label='Minor-allele carrier pairs',color='#b05c36')
            axes[1,1].set(title='Descriptive genetic similarity',xlabel='Chromosome-excluded relationship value',ylabel='Density');axes[1,1].legend(fontsize=8)
            fig.suptitle(region_id+' — '+anchor['marker']+' (Scheme A)')
            fig.subplots_adjust(top=.9,right=.86,hspace=.4,wspace=.4)
            fig.colorbar(scatter,cax=fig.add_axes([.89,.56,.015,.32]),label='Measured r² with display anchor')
            fig.savefig(output/(region_id+'.png'),dpi=150,bbox_inches='tight');plt.close(fig)
    write_table(output/'neighbor_ld.private.tsv',all_ld);write_table(output/'anchor_genotypes.private.tsv',all_counts);write_table(output/'candidate_models.private.tsv',all_models)
    summary=dict(status='candidate_regions_reviewed',config=cfg,flow=flow,regions=regions,
        reporting_policy_sha256=file_digest(ROOT/'config/gwas_reporting.json'),review_config_sha256=file_digest(ROOT/'config/candidate_review.json'),
        script_sha256=file_digest(Path(__file__)),scan_summary_sha256=file_digest(scan_results/'summary.private.json'),
        original_sources_unchanged=all(file_digest(study_path(ROOT,p))==h for p,h in plan['sources'].items()))
    (output/'summary.private.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print('Candidate review:',flow,flush=True)
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-dir',type=Path,default=ROOT/'data/gwas/measured-20261004')
    p.add_argument('--scan-results',type=Path,default=ROOT/'results/private/version2/07_breed_specific_gwas')
    p.add_argument('--results-dir',type=Path,default=ROOT/'results/private/version2/09_candidate_review')
    args=p.parse_args();review(args.run_dir,args.scan_results,args.results_dir)
