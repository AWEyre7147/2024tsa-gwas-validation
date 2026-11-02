"""Prepare measured LR historical-target genotypes and descriptive allele counts.

No association tests, effect-direction claims, or model selection are performed.
"""
from pathlib import Path
from collections import Counter
import argparse
import csv
import json
import subprocess
from notebook_support import read_json, study_path, verify_datasets, file_digest, plink_environment

ROOT=Path(__file__).resolve().parents[1]


def table(path):
    return [line.split() for line in path.read_text(encoding='utf-8').splitlines() if line.strip()]


def write_table(path, rows):
    with path.open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t')
        writer.writeheader();writer.writerows(rows)


def allele_counts(dosages):
    """Count ALT copies among observed diploid genotypes; leave missing calls out."""
    called=[value for value in dosages if value is not None]
    alternate=sum(called)
    return dict(called_dogs=len(called),missing_dogs=len(dosages)-len(called),
                alt_copies=alternate,total_alleles=2*len(called),
                alt_frequency=alternate/(2*len(called)) if called else None,
                carriers=sum(value>0 for value in called))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'data/06-replication/target-preparation-20261004')
    args=parser.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    alignment=read_json(ROOT/'data/02-aligned/manifest.private.json')
    assert alignment['status']=='complete';verify_datasets(ROOT,alignment)
    with (ROOT/'data/02-aligned/reports/replication_target_alignment.private.tsv').open(encoding='utf-8',newline='') as handle:
        targets=[row for row in csv.DictReader(handle,delimiter='\t') if row['breed']=='LR']
    assert len(targets)==18 and len({row['historical_target'] for row in targets})==18
    available=[row for row in targets if row['QC_status']=='retained']
    assert len(available)==12
    (out/'measured_targets.txt').write_text(''.join(row['input_marker']+'\n' for row in available),encoding='utf-8')
    environment=read_json(plink_environment(ROOT));exe=environment['executable']
    counts=[];sample_summaries=[];commands=[];matrices={}
    for scenario in ['success-vs-behelim','not-behelim-vs-behelim']:
        dataset=next(d for d in alignment['datasets'] if d['breed']=='LR' and d['scenario']==scenario)
        prefix=next(study_path(ROOT,p) for p in dataset['files'] if p.endswith('.fam')).with_suffix('')
        fam=table(prefix.with_suffix('.fam'));bim=table(prefix.with_suffix('.bim'))
        positions={row[1]:i for i,row in enumerate(bim)}
        assert len(positions)==len(bim)
        raw=prefix.with_suffix('.bed').read_bytes();stride=(len(fam)+3)//4
        assert raw[:3]==bytes([108,27,1]) and len(raw)==3+len(bim)*stride
        sex=Counter((row[5],row[4]) for row in fam)
        sample_summaries.append(dict(scenario=scenario,dogs=len(fam),cases=sum(r[5]=='2' for r in fam),
            controls=sum(r[5]=='1' for r in fam),case_males=sex['2','1'],case_females=sex['2','2'],
            control_males=sex['1','1'],control_females=sex['1','2'],unknown_sex=sum(r[4]=='0' for r in fam)))
        matrix=[dict(study_id=row[1],phenotype=int(row[5]),sex=int(row[4])) for row in fam]
        for target in available:
            marker=target['input_marker'];index=positions[marker];variant=bim[index]
            assert int(variant[3])==int(target['canFam4_position'])
            calls=[{0:2,1:None,2:1,3:0}[(raw[3+index*stride+j//4]>>(2*(j%4)))&3] for j in range(len(fam))]
            for row,value in zip(matrix,calls):row[marker]=value
            groups={name:allele_counts([value for r,value in zip(fam,calls) if r[5]==code])
                    for name,code in [('case','2'),('control','1')]}
            record=dict(scenario=scenario,historical_target=target['historical_target'],marker=marker,
                        canFam4_chromosome='chr'+variant[0],canFam4_position=int(variant[3]),
                        reference_allele=variant[5],alternate_allele=variant[4])
            for group,values in groups.items():record.update({group+'_'+key:value for key,value in values.items()})
            record['case_minus_control_alt_frequency']=groups['case']['alt_frequency']-groups['control']['alt_frequency']
            counts.append(record)
        write_table(out/(scenario+'_dosages.private.tsv'),matrix)
        matrices[scenario]={row['study_id']:row for row in matrix}
        export=out/(scenario+'_plink_check')
        argv=[exe,'--bfile',str(prefix),'--dog','--keep-allele-order','--extract',str(out/'measured_targets.txt'),
              '--recode','A','--out',str(export)]
        proc=subprocess.run(argv,capture_output=True,text=True)
        (out/(scenario+'_plink_check.console.txt')).write_text(proc.stdout+proc.stderr,encoding='utf-8')
        proc.check_returncode();commands.append(dict(argv=argv,returncode=proc.returncode))
        exported=table(export.with_suffix('.raw'));headers=exported[0]
        assert [row[1] for row in exported[1:]]==[row[1] for row in fam]
        for target in available:
            marker=target['input_marker'];allele=bim[positions[marker]][4]
            column=headers.index(marker+'_'+allele)
            observed=[None if row[column]=='NA' else int(row[column]) for row in exported[1:]]
            assert observed==[row[marker] for row in matrix],marker
    primary=matrices['success-vs-behelim'];broader=matrices['not-behelim-vs-behelim']
    assert primary.keys()<=broader.keys() and all(row==broader[identifier] for identifier,row in primary.items())
    assert {key for key,row in primary.items() if row['phenotype']==2}=={key for key,row in broader.items() if row['phenotype']==2}
    verify_datasets(ROOT,alignment)
    write_table(out/'target_availability.private.tsv',targets)
    write_table(out/'allele_counts.private.tsv',counts)
    result=dict(status='measured_target_inputs_prepared',breed='LR',historical_targets=18,measured_targets=12,
        excluded_targets=6,cohorts=sample_summaries,allele_count_rows=len(counts),
        observed_dosages_independently_verified_with_plink=True,association_tests_performed=False,
        historical_direction_status='Requires allele-specific source review; not inferred from ALT labels',
        model_status='awaiting candidate-model choice and calibration',
        source_alignment_manifest_sha256=file_digest(ROOT/'data/02-aligned/manifest.private.json'),commands=commands)
    (out/'manifest.private.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='commands'},indent=2))


if __name__=='__main__':main()
