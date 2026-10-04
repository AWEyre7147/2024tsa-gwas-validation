"""Independently export aligned BEDs with PLINK and compare every VCF call."""
from pathlib import Path
import gzip
import json
import subprocess
from reference_alignment import ROOT, OUT, sha, dump


def records(f):
    for line in f:
        if line.startswith('#CHROM'): yield line.rstrip().split('\t')
        elif not line.startswith('#'): yield line.rstrip().split('\t')


def main():
    m=json.loads((OUT/'manifest.private.json').read_text())
    env=json.loads((ROOT/'data/processed/private/plink_environment.private.json').read_text())
    exe=Path(env['executable']);assert sha(exe)==env['executable_sha256']
    results=[]
    for ds in m['datasets']:
        assert all(sha(ROOT/p)==h for p,h in ds['files'].items())
        prefix=next(ROOT/p for p in ds['files'] if p.endswith('.bed')).with_suffix('')
        dest=OUT/'reports'/('plink_check_'+ds['scenario']+'_'+ds['breed'])
        args=[str(exe),'--dog','--bfile',str(prefix),'--keep-allele-order','--real-ref-alleles','--recode','vcf-iid','--out',str(dest)]
        result=subprocess.run(args,capture_output=True,text=True);result.check_returncode()
        count=0
        with gzip.open(prefix.with_suffix('.vcf.gz'),'rt') as source, dest.with_suffix('.vcf').open() as check:
            for a,b in zip(records(source),records(check),strict=True):
                if a[0]=='#CHROM': assert a[9:]==b[9:];continue
                assert a[0].removeprefix('chr')==b[0] and a[1:5]==b[1:5]
                assert [sorted(x.split(':')[0].replace('|','/').split('/')) for x in a[9:]]==[sorted(x.split(':')[0].replace('|','/').split('/')) for x in b[9:]]
                count+=1
        assert count==ds['aligned_variants']
        results.append(dict(scenario=ds['scenario'],breed=ds['breed'],markers=count,all_calls_agree=True,command=args,export_sha256=sha(dest.with_suffix('.vcf'))))
        # Temporary independent exports are reproducible; keep the small logs.
        dest.with_suffix('.vcf').unlink()
        print(ds['breed'],ds['scenario'],'all PLINK-exported calls agree',flush=True)
    dump(OUT/'reports/plink_verification.private.json',dict(status='passed',plink=env['version'],datasets=results))


if __name__=='__main__':main()
