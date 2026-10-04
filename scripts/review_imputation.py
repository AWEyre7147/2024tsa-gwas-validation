"""Review held-out marker accuracy and directly measured replication targets.

Run with the tsa2024 Python environment in WSL. Original inference files are
read only; review outputs contain private genotype-derived information.
"""
from pathlib import Path
import csv
import argparse
import json
import subprocess
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / 'data/imputation/autosomal-20261003'
OUT = RUN / 'review'
BCF = Path.home() / 'anaconda3/envs/tsa2024-imputation/bin/bcftools'


def query(path, positions, label):
    bed = OUT / (label + '.bed')
    bed.write_text(''.join(f'{chrom}\t{pos-1}\t{pos}\n' for chrom, pos in positions), encoding='utf-8')
    output = subprocess.check_output([str(BCF), 'query', '-R', str(bed), '-f',
        '%CHROM\t%POS\t%REF\t%ALT\t%INFO/DR2[\t%GT:%DS]\n', str(path)], text=True)
    return [line.split('\t') for line in output.splitlines()]


def main():
    global RUN, OUT, BCF
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir',type=Path,default=RUN)
    parser.add_argument('--bcftools',type=Path,default=BCF)
    args=parser.parse_args()
    RUN=args.run_dir.resolve(); OUT=RUN/'review'
    BCF=args.bcftools
    OUT.mkdir(exist_ok=True)
    metrics = []
    for breed in ['GSP', 'LR']:
        for chrom in [1, 3, 13, 20, 36]:
            region = f'chr{chrom}'
            truth = json.loads((RUN / 'targets' / breed / region / 'truth.private.json').read_text(encoding='utf-8'))
            stored = json.loads((RUN / 'validation' / breed / region / 'masked_variant_metrics.private.json').read_text(encoding='utf-8'))
            bypos = {str(row['position']): row for row in stored}
            predictions = query(RUN / 'validation' / breed / region / 'masked.vcf.gz',
                [(region, int(pos)) for pos in truth], f'{breed}_{region}_hidden')
            assert len(predictions) == len(truth)
            for row in predictions:
                pos = row[1]; t = truth[pos]
                assert row[2:4] == [t['ref'], t['alt']]
                assert len(row[5:]) == len(t['gt'])
                pairs = [(sum(map(int, gt.replace('|', '/').split('/'))), float(pred.split(':')[1]))
                    for gt, pred in zip(t['gt'], row[5:]) if '.' not in gt]
                obs, ds = np.asarray(pairs, dtype=float).T
                corr = float(np.corrcoef(obs, ds)[0, 1]) if np.std(obs) > 0 and np.std(ds) > 0 else None
                m = dict(bypos[pos], breed=breed, chromosome=region,
                    signed_correlation=corr, dosage_rmse=float(np.sqrt(np.mean((obs-ds)**2))))
                if corr is not None: assert abs(corr*corr-m['dosage_r2']) < 1e-7
                metrics.append(m)
    summaries = []
    for breed in ['GSP', 'LR']:
        for cutoff in [0, .8, .9, .95, .99]:
            selected = [m for m in metrics if m['breed']==breed and m['reported_DR2'] is not None and m['reported_DR2']>=cutoff]
            defined = [m for m in selected if m['dosage_r2'] is not None]
            summaries.append(dict(breed=breed, cutoff=cutoff, retained=len(selected), measurable=len(defined),
                median_r2=float(np.median([m['dosage_r2'] for m in defined])),
                below_08=sum(m['dosage_r2']<.8 for m in defined),
                negative_correlations=sum(m['signed_correlation']<0 for m in defined)))
    targets = list(csv.DictReader((ROOT/'data/aligned/reports/replication_target_alignment.private.tsv').open(encoding='utf-8'), delimiter='\t'))
    for target in targets:
        target['review_status'] = 'Excluded during QC; no imputation rescue assumed'
        if target['QC_status'] != 'retained': continue
        breed = target['breed']; chrom = target['historical_target'].split('.')[0]; pos = int(target['canFam4_position'])
        rows = query(RUN/'not-behelim-vs-behelim'/breed/chrom/'imputed.vcf.gz', [(chrom,pos)], f'{breed}_{target["input_marker"]}')
        assert len(rows)==1 and int(rows[0][1])==pos
        target['review_status'] = 'Directly measured target present in completed output'
        target['reference_allele'], target['alternate_allele'] = rows[0][2:4]
    result = dict(confidence_thresholds=summaries, heldout_marker_count=len(metrics),
        negative_correlations=sum(m['signed_correlation'] is not None and m['signed_correlation']<0 for m in metrics),
        targets=targets, recommendation='Use QC-passed observed target genotypes for primary replication. DR2 alone does not establish imputed variant accuracy. Effect alleles and directions require historical source confirmation.')
    (OUT/'review.private.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    (OUT/'heldout_markers.private.json').write_text(json.dumps(metrics,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
