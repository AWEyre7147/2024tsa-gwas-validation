"""Analytical benchmark and native-model, finite-cohort power/null simulations.

Default pilot checks the framework with explicit Monte Carlo uncertainty. Full
runs are opt-in and use a new output directory. No whole-genome inference or
existing association receipts/configs are changed.
"""
from pathlib import Path
import argparse
import hashlib
import math
import shutil
import subprocess
import numpy as np
from scipy.stats import norm
from breed_gwas import decode
from extension_support import ROOT, save, sources, start, finish, rows, tsv
from notebook_support import read_json,file_digest,study_path

TESTS=['GMMAT_score','GEMMA_wald','GEMMA_lrt','GEMMA_score']


def benchmark_power(frequency, OR, cases, controls, alpha):
    # HWE controls with rare-outcome approximation: retrospective case genotype
    # frequencies are proportional to OR^dosage times control probabilities.
    control=np.array([(1-frequency)**2,2*frequency*(1-frequency),frequency**2])
    case=control*np.power(OR,np.arange(3));case/=case.sum()
    pcase=float(case@np.arange(3)/2)
    se=math.sqrt(pcase*(1-pcase)/(2*cases)+frequency*(1-frequency)/(2*controls))
    delta=(pcase-frequency)/se
    critical=norm.isf(alpha/2)
    return float(norm.sf(critical-delta)+norm.cdf(-critical-delta)),pcase


def wilson(k,n):
    if not n:return None,None
    z=norm.ppf(.975);p=k/n;denom=1+z*z/n
    center=(p+z*z/(2*n))/denom
    half=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/denom
    return float(max(0,center-half)),float(min(1,center+half))


def verified_files(folder,names,expected):
    for name in names:
        if file_digest(folder/name)!=expected[name]:raise ValueError(f'Completed native source changed: {name}')


def collect_contexts(config):
    records=[];paths=[]
    manifest=read_json(ROOT/'data/06-replication/models-aligned-20261004/manifest.private.json')
    for job in manifest['jobs']:
        folder=(ROOT/job['result']).parent
        meta=read_json(folder/'input.private.json')
        records.append(dict(context='target_LR_'+job['scenario']+'_'+job['marker'],family='historical',breed='LR',
            scheme=job['scenario'],marker=job['marker'],native_folder=folder,
            ALT=meta['alternate_allele'],REF=meta['reference_allele']))
    alignment=read_json(ROOT/'data/02-aligned/manifest.private.json')
    for breed in ['LR','GSP']:
        for scheme in config['schemes']:
            prefix=ROOT/'data/02-aligned'/scheme/breed/('2024TSA_'+breed)
            dataset=next(d for d in alignment['datasets'] if d['breed']==breed and d['scenario']==scheme)
            for p,digest in dataset['files'].items():
                if Path(p).suffix in {'.bed','.bim','.fam'}:
                    path=study_path(ROOT,p);paths.append(path)
                    if file_digest(path)!=digest:raise ValueError('Aligned input hash mismatch')
            dogs,bim,g=decode(prefix);maf=np.minimum(np.nanmean(g,axis=1)/2,1-np.nanmean(g,axis=1)/2)
            for lo,hi in config['representative_maf_bins']:
                eligible=[i for i in range(len(bim)) if lo<=maf[i]<hi]
                if not eligible:continue
                selected=min(eligible,key=lambda i:hashlib.sha256((str(config['seed'])+breed+bim[i][1]).encode()).hexdigest())
                row=bim[selected];folder=ROOT/'data/07-gwas/measured-20261004'/breed/scheme/f'chr{int(row[0]):02d}'
                receipt=read_json(folder/'receipt.private.json')
                verified_files(folder,['dogs.private.tsv','kinship.txt','gmmat_null.private.rds','markers.private.tsv'],receipt['files'])
                observed=next(r for r in rows(folder/'markers.private.tsv') if r['marker']==row[1])
                native_dogs=rows(folder/'dogs.private.tsv')
                if [r['id'] for r in native_dogs]!=[r[1] for r in dogs]:raise ValueError('Native GWAS sample order differs')
                filled=np.where(np.isnan(g[selected]),np.nanmean(g[selected]),g[selected])
                records.append(dict(context='GWAS_'+breed+'_'+scheme+'_'+row[1],family='GWAS',breed=breed,scheme=scheme,
                    marker=row[1],native_folder=folder,ALT=row[4],REF=row[5],g=filled,
                    MAF=float(maf[selected]),observed_called=int(np.isfinite(g[selected]).sum()),
                    maf_bin=[lo,hi]))
    return records,paths


def summarize(folder,contexts,policy):
    summaries=[]
    for context in contexts:
        path=folder/context['context']/'draws.private.tsv'
        draws=rows(path)
        cells=sorted({(r['OR'],r['variance_multiplier']) for r in draws},key=lambda x:(float(x[1]),float(x[0])))
        family_alpha=.05/36 if context['family']=='historical' else policy['family_alpha']/policy['candidate_markers'][context['breed']]
        for OR,variance in cells:
            subset=[r for r in draws if (r['OR'],r['variance_multiplier'])==(OR,variance)]
            for test in TESTS:
                values=[]
                for r in subset:
                    try:p=float(r[test])
                    except ValueError:continue
                    if math.isfinite(p) and 0<=p<=1:values.append(p)
                for label,alpha in [('nominal',.05),('analysis_family',family_alpha)]:
                    successes=sum(p<alpha for p in values);valid=len(values);total=len(subset)
                    low,high=wilson(successes,valid)
                    summaries.append(dict(context=context['context'],family=context['family'],breed=context['breed'],scheme=context['scheme'],
                        marker=context['marker'],generating_OR=float(OR),variance_multiplier=float(variance),test=test,
                        threshold_label=label,alpha=alpha,requested=total,valid=valid,failed=total-valid,
                        rejections=successes,rate_valid=successes/valid if valid else None,Wilson_low=low,Wilson_high=high,
                        failure_lower=successes/total,failure_upper=(successes+total-valid)/total,
                        null_scenario=float(OR)==1))
    tsv(folder/'rates.private.tsv',summaries)
    return summaries


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile',choices=['pilot','full'],default='pilot')
    parser.add_argument('--replicates',type=int)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--tools',type=Path,default=Path.home()/'anaconda3/envs/tsa2024-association/bin')
    args=parser.parse_args();config=read_json(ROOT/'config/power.json');policy=read_json(ROOT/'config/gwas_reporting.json')
    contexts,paths=collect_contexts(config)
    if args.profile=='pilot':
        chosen=[next(r for r in contexts if r['family']=='historical' and r['scheme']=='success-vs-behelim')]
        chosen += [next(r for r in contexts if r['family']=='GWAS' and r['breed']==b and r['scheme']=='success-vs-behelim') for b in ['LR','GSP']]
        contexts=chosen
    replicates=args.replicates or config[args.profile]['replicates']
    if replicates<1:raise ValueError('Replicates must be positive')
    odds=config['pilot']['odds_ratios'] if args.profile=='pilot' else config['odds_ratios']
    variances=config['pilot']['variance_multipliers'] if args.profile=='pilot' else config['variance_multipliers']
    if args.profile=='full' and args.output is None:raise ValueError('Supply a new --output directory for a full run')
    paths += [ROOT/'config/power.json',ROOT/'config/gwas_reporting.json',ROOT/'config/replication.json',Path(__file__),
              ROOT/'scripts/power_simulation.R',ROOT/'scripts/breed_gwas.py',ROOT/'scripts/replication_models.py',
              ROOT/'scripts/extension_support.py',ROOT/'scripts/notebook_support.py']
    for r in contexts:
        receipt_name='result-receipt.private.json' if r['family']=='historical' else 'receipt.private.json'
        paths += [r['native_folder']/name for name in ['input.private.json','dogs.private.tsv','kinship.txt','gmmat_null.private.rds',receipt_name]]
        if r['family']=='historical':
            receipt=read_json(r['native_folder']/receipt_name)
            result_path=r['native_folder']/'result.private.json'
            if file_digest(result_path)!=receipt['sha256']:raise ValueError('Historical native result receipt differs')
            paths.append(result_path)
            paths += [r['native_folder']/name for name in ['genotypes.txt','gemma_genotypes.txt']]
    versions=dict(R=subprocess.check_output([str(args.tools/'Rscript'),'--version'],text=True).strip(),
        GMMAT=subprocess.check_output([str(args.tools/'Rscript'),'-e','cat(as.character(packageVersion("GMMAT")))'],text=True),
        GEMMA_sha256=file_digest(args.tools/'gemma'))
    plan=dict(config=config,profile=args.profile,replicates=replicates,OR=odds,variance_multipliers=variances,
        context_ids=[r['context'] for r in contexts],sources=sources(paths),versions=versions,
        numpy=np.__version__,scipy=__import__('scipy').__version__)
    folder=args.output or ROOT/config['output_directory']
    if not start(folder,plan):return
    benchmarks=[]
    for breed,scheme,cases,controls in [('LR','A',7,39),('LR','B',7,41),('GSP','A',9,96),('GSP','B',9,107)]:
        thresholds=[('nominal',.05),('GWAS',policy['family_alpha']/policy['candidate_markers'][breed])]
        if breed=='LR':thresholds.append(('historical',.05/36))
        for p in config['analytical_benchmark']['maf']:
            for OR in config['odds_ratios']:
                for label,alpha in thresholds:
                    power,case_af=benchmark_power(p,OR,cases,controls,alpha)
                    benchmarks.append(dict(breed=breed,scheme=scheme,cases=cases,controls=controls,control_frequency=p,
                        OR=OR,case_frequency=case_af,threshold_label=label,alpha=alpha,power=power))
    tsv(folder/'analytical_benchmark.tsv',benchmarks)
    public_contexts=[]
    for r in contexts:
        dest=folder/r['context'];dest.mkdir();native=r['native_folder']
        dogs=rows(native/'dogs.private.tsv');n=len(dogs);cases=sum(int(d['y']) for d in dogs)
        for name in ['dogs.private.tsv','kinship.txt']:
            shutil.copyfile(native/name,dest/name)
        if r['family']=='historical':
            for name in ['genotypes.txt','gemma_genotypes.txt']:
                shutil.copyfile(native/name,dest/name)
            g=np.array([float(x) for x in (dest/'genotypes.txt').read_text().split()[1:]])
            r['MAF']=float(min(g.mean()/2,1-g.mean()/2));r['observed_called']=n
        else:
            g=r['g'];values=[format(x,'.16g') for x in g]
            (dest/'genotypes.txt').write_text(r['marker']+'\t'+'\t'.join(values)+'\n')
            (dest/'gemma_genotypes.txt').write_text(','.join([r['marker'],r['ALT'],r['REF'],*values])+'\n')
        K=np.loadtxt(dest/'kinship.txt')
        if K.shape!=(n,n) or np.linalg.eigvalsh(K).min() < -1e-8:raise ValueError('Invalid matched LOCO matrix')
        np.savetxt(dest/'covariates.txt',np.column_stack([np.ones(n),[int(d['male']) for d in dogs]]),fmt='%d',delimiter='\t')
        seed=int(hashlib.sha256((str(config['seed'])+r['context']).encode()).hexdigest()[:7],16)
        meta=dict(context=r['context'],cases=cases,seed=seed,replicates=replicates,odds_ratios=odds,variance_multipliers=variances,
                  maximum_attempts=config['maximum_draw_attempts_per_replicate'],ALT=r['ALT'],REF=r['REF'],
                  null_model=str(native/'gmmat_null.private.rds'),gemma=str(args.tools/'gemma'))
        save(dest/'input.private.json',meta)
        print('Simulating',r['context'],flush=True)
        with (dest/'simulation.console.txt').open('w') as log:
            subprocess.run([str(args.tools/'Rscript'),str(ROOT/'scripts/power_simulation.R'),str(dest)],stdout=log,stderr=subprocess.STDOUT,
                check=True,env={**__import__('os').environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'})
        public_contexts.append({k:v for k,v in r.items() if k not in {'g','native_folder'}})
        print('Completed',r['context'],flush=True)
    save(folder/'contexts.private.json',public_contexts)
    rates=summarize(folder,contexts,policy)
    finish(folder,plan,dict(status='power_framework_pilot_complete' if args.profile=='pilot' else 'power_framework_run_complete',
        contexts=len(contexts),replicates_per_cell=replicates,total_simulated_outcomes=len(contexts)*len(odds)*len(variances)*replicates,
        rate_rows=len(rates),sources_unchanged=True,extreme_tail_validated=False,
        limitations=['Conditional on observed cohorts, genotype patterns and stated logistic generating model.',
                     'Illustrative OR grid is not a historical effect estimate or observed-effect power calculation.',
                     'Analytical benchmark assumes independent HWE controls and a rare-outcome approximation.',
                     'Finite pilot and ordinary null replicate counts do not validate genome-wide tail probabilities.',
                     'Full profile may be computationally expensive; it is not launched by notebook review.']))


if __name__=='__main__':main()
