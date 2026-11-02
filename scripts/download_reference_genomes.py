"""Cache the public CanFam3.1 and CanFam4 genomes with published MD5 checks."""
from pathlib import Path
import concurrent.futures
import hashlib
import json
import shutil
import subprocess

ROOT=Path(__file__).resolve().parents[1]
DEST=ROOT/'data/02-reference/genomes'


def fetch(assembly):
    curl=shutil.which('curl.exe') or shutil.which('curl')
    base=f'https://hgdownload.soe.ucsc.edu/goldenPath/{assembly}/bigZips/'
    for name in (assembly+'.2bit',assembly+'.md5sum.txt'):
        p=DEST/name
        if not p.exists():
            temp=p.with_name(p.name+'.partial')
            url=base+('md5sum.txt' if name.endswith('.txt') else name)
            subprocess.run([curl,'--fail','--location','--silent','--show-error','--retry','2',url,'--output',str(temp)],check=True)
            temp.replace(p)
    expected=next(line.split()[0] for line in (DEST/(assembly+'.md5sum.txt')).read_text().splitlines() if line.split()[-1].lstrip('*')==assembly+'.2bit')
    md5=hashlib.md5();sha=hashlib.sha256();p=DEST/(assembly+'.2bit')
    with p.open('rb') as f:
        while chunk:=f.read(1048576):md5.update(chunk);sha.update(chunk)
    assert md5.hexdigest()==expected,'Reference checksum mismatch'
    return assembly,dict(url=base+assembly+'.2bit',checksum_url=base+'md5sum.txt',bytes=p.stat().st_size,md5=expected,sha256=sha.hexdigest(),published_checksum_verified=True)


if __name__=='__main__':
    DEST.mkdir(parents=True,exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool: report=dict(pool.map(fetch,['canFam3','canFam4']))
    (DEST/'assets.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Both public reference genomes match published checksums.')
