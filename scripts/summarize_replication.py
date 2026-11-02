"""Summarize paired targeted models without declaring directional replication."""
from pathlib import Path
import argparse
import csv
import json
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from notebook_support import read_json, file_digest, study_path

ROOT = Path(__file__).resolve().parents[1]


def wilson(successes, total):
    """Binomial uncertainty for independent null draws at a single target."""
    if not total:
        return None, None
    z = 1.959963984540054
    p = successes / total
    denominator = 1 + z*z/total
    center = (p + z*z/(2*total)) / denominator
    half = z * math.sqrt(p*(1-p)/total + z*z/(4*total*total)) / denominator
    return max(0,center-half), min(1,center+half)


def one(value):
    return value[0] if isinstance(value, list) else value


def summarize(run, output):
    manifest = read_json(run / 'manifest.private.json')
    assert manifest['status'] == 'paired_models_complete'
    plan = read_json(run / 'plan.private.json')
    assert all(file_digest(study_path(ROOT, p)) == h for p, h in plan['sources'].items())
    rows = []
    for job in manifest['jobs']:
        path = study_path(ROOT, job['result'])
        assert file_digest(path) == read_json(path.parent/'result-receipt.private.json')['sha256']
        result = read_json(path)
        meta = read_json(path.parent/'input.private.json')
        score = one(result['GMMAT_score']); linear = one(result['GEMMA'])
        wald = one(result['GMMAT_wald'])
        assert result['null_converged']
        assert score['SNP'] == linear['rs'] == job['marker']
        assert linear['allele1'] == meta['alternate_allele']
        assert linear['allele0'] == meta['reference_allele']
        with (path.parent/'dogs.private.tsv').open(encoding='utf-8', newline='') as handle:
            dogs = list(csv.DictReader(handle, delimiter='\t'))
        measured_af = sum(float(d['genotype']) for d in dogs)/(2*len(dogs))
        assert score['N'] == len(dogs) == meta['complete_dogs']
        assert abs(score['AF']-measured_af) < 1e-6
        assert abs(linear['af']-measured_af) < 0.0006  # GEMMA prints three decimals.
        simulation = result['calibration']
        assert simulation['completed']+simulation['failed'] == simulation['requested']
        for method, p in [('GMMAT', score['PVAL']), ('GEMMA', linear['p_score'])]:
            assert p is not None and math.isfinite(p) and 0 <= p <= 1
            ps = [v[method] for v in simulation['values']]
            assert len(ps) == simulation['completed']
            hits = sum(v < .05 for v in ps)
            low, high = wilson(hits, len(ps))
            converged = bool(wald.get('converged', False)) if method == 'GMMAT' else None
            beta = wald.get('BETA') if method == 'GMMAT' else linear['beta']
            se = wald.get('SE') if method == 'GMMAT' else linear['se']
            valid_effect = beta is not None and se is not None and math.isfinite(beta) and math.isfinite(se) and se > 0
            if method == 'GMMAT':
                valid_effect = valid_effect and converged
            rows.append(dict(scenario=job['scenario'], marker=job['marker'],
                historical_target=meta['historical_target'], counted_allele=meta['alternate_allele'],
                reference_allele=meta['reference_allele'], dogs=len(dogs), cases=meta['case_count'],
                background_markers=meta['background_markers'], method=method, score_p=p,
                beta=beta, standard_error=se, effect_scale='log odds' if method=='GMMAT' else '0/1 outcome',
                effect_fit_converged=converged, effect_available=valid_effect,
                wald_p=wald.get('PVAL') if method=='GMMAT' else linear['p_wald'],
                gemma_lrt_p=linear.get('p_lrt') if method=='GEMMA' else None,
                gemma_lambda_mle=linear.get('l_mle') if method=='GEMMA' else None,
                gemma_lambda_at_search_boundary=(linear.get('l_mle',1)<=1.00001e-5 or linear.get('l_mle',1)>=99999) if method=='GEMMA' else None,
                null_replicates_requested=simulation['requested'], null_replicates_valid=len(ps),
                null_replicates_failed=simulation['failed'], null_p_below_05=hits,
                null_rate_below_05=hits/len(ps) if ps else None,
                null_rate_ci_low=low, null_rate_ci_high=high,
                null_rate_failure_lower=hits/simulation['requested'] if simulation['requested'] else None,
                null_rate_failure_upper=(hits+simulation['failed'])/simulation['requested'] if simulation['requested'] else None))
    output.mkdir(parents=True, exist_ok=True)
    with (output/'paired_models.private.tsv').open('w', encoding='utf-8', newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t')
        writer.writeheader();writer.writerows(rows)
    threshold = manifest['config']['family_alpha']/(manifest['config']['historical_target_count']*2)
    summary = dict(status='paired_models_summarized', directional_replication_assessed=False,
        family_comparison_threshold=threshold, tested_targets_per_cohort=12,
        observed_null_models_converged=24, null_draws_requested=sum(r['null_replicates_requested'] for r in rows[::2]),
        paired_null_draws_valid=sum(r['null_replicates_valid'] for r in rows[::2]),
        paired_null_draws_failed=sum(r['null_replicates_failed'] for r in rows[::2]),
        gmmat_effect_fits_converged=sum(r['effect_fit_converged'] for r in rows if r['method']=='GMMAT'),
        gemma_alternative_mle_at_search_boundary=sum(r['gemma_lambda_at_search_boundary'] for r in rows if r['method']=='GEMMA'),
        primary={method:dict(minimum_score_p=min(r['score_p'] for r in rows if r['method']==method and r['scenario']=='success-vs-behelim'),
            nominal_below_05=sum(r['score_p']<.05 for r in rows if r['method']==method and r['scenario']=='success-vs-behelim'),
            below_family_threshold=sum(r['score_p']<threshold for r in rows if r['method']==method and r['scenario']=='success-vs-behelim')) for method in ['GMMAT','GEMMA']},
        result_hashes={job['result']:file_digest(study_path(ROOT,job['result'])) for job in manifest['jobs']})
    (output/'summary.private.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(13,6),sharey=True)
    for ax,scenario,title in zip(axes,['success-vs-behelim','not-behelim-vs-behelim'],['Scheme A: primary','Scheme B: sensitivity']):
        selected=[r for r in rows if r['scenario']==scenario]
        targets=[r['historical_target'] for r in selected if r['method']=='GMMAT']
        for method,offset,color in [('GMMAT',-.10,'#246b8f'),('GEMMA',.10,'#c27c2c')]:
            values=[r for r in selected if r['method']==method]
            ax.scatter([-math.log10(max(r['score_p'],1e-300)) for r in values],np.arange(12)+offset,label=method,color=color,s=35)
        ax.axvline(-math.log10(.05),color='.5',ls=':',label='Nominal 0.05')
        ax.axvline(-math.log10(threshold),color='.4',ls='--',label='18 targets × 2 methods')
        ax.set(title=title,xlabel='−log₁₀(score-test p-value)',yticks=np.arange(12),yticklabels=targets)
        ax.grid(axis='y',alpha=.15)
    axes[0].invert_yaxis();axes[1].legend(loc='lower right',fontsize=8)
    fig.suptitle('Measured Labrador targets: paired mixed-model comparison')
    fig.tight_layout();fig.savefig(output/'paired_score_tests.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(13,6),sharey=True)
    for ax,scenario,title in zip(axes,['success-vs-behelim','not-behelim-vs-behelim'],['Scheme A: primary','Scheme B: sensitivity']):
        selected=[r for r in rows if r['scenario']==scenario]
        for method,offset,color in [('GMMAT',-.10,'#246b8f'),('GEMMA',.10,'#c27c2c')]:
            values=[r for r in selected if r['method']==method]
            rates=np.array([r['null_rate_below_05'] for r in values])
            lows=np.array([r['null_rate_ci_low'] for r in values]);highs=np.array([r['null_rate_ci_high'] for r in values])
            ax.errorbar(rates,np.arange(12)+offset,xerr=np.array([rates-lows,highs-rates]),fmt='o',color=color,label=method,markersize=4,capsize=2)
        ax.axvline(.05,color='.4',ls='--');ax.set(title=title,xlabel='Fraction of valid null draws with p < 0.05',yticks=np.arange(12),yticklabels=[r['historical_target'] for r in selected if r['method']=='GMMAT'])
        ax.grid(axis='y',alpha=.15)
    axes[0].set_yticks(np.arange(12), labels=[r['historical_target'] for r in rows if r['scenario']=='success-vs-behelim' and r['method']=='GMMAT'])
    axes[0].tick_params(axis='y', labelleft=True)
    axes[0].invert_yaxis();axes[1].legend(loc='lower right')
    fig.suptitle('Initial null screen: 95% binomial intervals, converged paired draws only')
    fig.tight_layout();fig.savefig(output/'null_screen.png',dpi=180);plt.close(fig)
    print(json.dumps({k:v for k,v in summary.items() if k!='result_hashes'},indent=2))
    return summary


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir',type=Path,default=ROOT/'data/06-replication/models-aligned-20261004')
    parser.add_argument('--results-dir',type=Path,default=ROOT/'results/private/version2/06_2013tsa_lr_replication')
    args=parser.parse_args();summarize(args.run_dir,args.results_dir)
