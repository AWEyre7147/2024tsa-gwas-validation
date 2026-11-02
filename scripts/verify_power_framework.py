"""Verify simulation journals, sample matching and independently recompute rates/CIs."""
from pathlib import Path
import argparse
import math
from scipy.stats import norm
from extension_support import ROOT,rows,save
from notebook_support import read_json,file_digest


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run',type=Path)
    args=parser.parse_args();config=read_json(ROOT/'config/power.json');folder=args.run or ROOT/config['output_directory']
    receipt=read_json(folder/'receipt.private.json')
    for name,digest in receipt['outputs'].items():
        if file_digest(folder/name)!=digest:raise ValueError('Power output receipt mismatch')
    contexts=read_json(folder/'contexts.private.json');groups={};outcomes=0;failed_generated=0;boundary_count=0
    for context in contexts:
        dest=folder/context['context'];meta=read_json(dest/'input.private.json');dogs=rows(dest/'dogs.private.tsv')
        native=Path(meta['null_model']).parent
        for name in ['dogs.private.tsv','kinship.txt']:
            if file_digest(dest/name)!=file_digest(native/name):raise ValueError('Power sample/matrix differs from original native input')
        draws=rows(dest/'draws.private.tsv')
        for r in draws:
            key=(r['context'],float(r['OR']),float(r['variance_multiplier']))
            groups.setdefault(key,[]).append(r)
            if r['generated']=='TRUE':
                if len(r['outcome_bits'])!=len(dogs) or set(r['outcome_bits'])-{'0','1'} or r['outcome_bits'].count('1')!=meta['cases']:
                    raise ValueError('Generated phenotype journal does not preserve case count/order')
                if abs(float(r['expected_cases'])-meta['cases'])>1e-6:raise ValueError('Generating-model expected count differs')
                outcomes+=1
            else:failed_generated+=1
            boundary_count+=r['GEMMA_mle_boundary']=='TRUE'
    checked=0;total_failures={}
    for r in rows(folder/'rates.private.tsv'):
        subset=groups[r['context'],float(r['generating_OR']),float(r['variance_multiplier'])]
        values=[]
        for x in subset:
            try:p=float(x[r['test']])
            except ValueError:continue
            if math.isfinite(p) and 0<=p<=1:values.append(p)
        n=len(values);total=len(subset);k=sum(p<float(r['alpha']) for p in values)
        if (n,total,k)!=(int(r['valid']),int(r['requested']),int(r['rejections'])):raise ValueError('Rejection counts disagree')
        if n:
            # Independent quadratic inversion of binomial score inequality.
            z2=float(norm.ppf(.975)**2);a=n+z2;b=-(2*k+z2);c=k*k/n
            disc=max(0,b*b-4*a*c);lo=(-b-math.sqrt(disc))/(2*a);hi=(-b+math.sqrt(disc))/(2*a)
            if abs(lo-float(r['Wilson_low']))>1e-10 or abs(hi-float(r['Wilson_high']))>1e-10:
                raise ValueError('Independent binomial interval calculation disagrees')
        if r['threshold_label']=='nominal':total_failures[r['test']]=total_failures.get(r['test'],0)+total-n
        checked+=1
    summary=dict(status='passed',generated_outcomes_checked=outcomes,sampler_failures=failed_generated,
        rate_rows_independently_checked=checked,failed_tests_by_method=total_failures,
        GEMMA_boundary_draws=boundary_count,extreme_tail_validated=False,
        receipt_sha256=file_digest(folder/'receipt.private.json'),verifier_sha256=file_digest(Path(__file__)))
    save(folder/'verification.private.json',summary);print(summary)


if __name__=='__main__':main()
