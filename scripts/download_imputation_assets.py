"""Download public assets for a bounded Dog10K pilot; no genotype upload."""
from pathlib import Path
import concurrent.futures
import hashlib
import json
import subprocess
import shutil
import urllib.request

ROOT=Path(__file__).resolve().parents[1]
DEST=ROOT/'data/02-reference/imputation-pilot'
ASSETS={
 'canFam3ToCanFam4.over.chain.gz':'https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/canFam3ToCanFam4.over.chain.gz',
 'canFam3ToROS.over.chain.gz':'https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/canFam3ToGCF_014441545.1.over.chain.gz',
 'beagle.27Feb25.75f.jar':'https://faculty.washington.edu/browning/beagle/beagle.27Feb25.75f.jar',
 'campbell_sex_average_canFam4.tar.gz':'https://zenodo.org/records/17095604/files/campbell_sex_average_canFam4.tar.gz?download=1',
 'samples.txt':'https://zenodo.org/records/8084059/files/dog10K-alignment-sample-table.2022-02-23-v8.txt.UPDATESRA.txt?download=1',
 'panel.md5sums.txt':'https://kiddlabshare.med.umich.edu/dog10K/phased-imputation-panel/md5sums.txt',
 'AutoAndXPAR.Dog10K.phased.bcf.csi':'https://kiddlabshare.med.umich.edu/dog10K/phased-imputation-panel/AutoAndXPAR.Dog10K.phased.bcf.csi',
 'twoBitToFa':'https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa',
}

def fetch(item):
 name,url=item;p=DEST/name
 if not p.exists():
  tmp=p.with_name(p.name+'.partial')
  downloader=shutil.which('curl.exe') or shutil.which('curl')
  if not downloader:raise RuntimeError('curl is required for certificate-verified downloads')
  subprocess.run([downloader,'--fail','--location','--silent','--show-error','--retry','2','--max-time','180',url,'--output',str(tmp)],check=True)
  tmp.replace(p)
 return name,dict(url=url,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())

if __name__=='__main__':
 DEST.mkdir(parents=True,exist_ok=True)
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
  entries=dict(pool.map(fetch,ASSETS.items()))
 assert hashlib.md5((DEST/'AutoAndXPAR.Dog10K.phased.bcf.csi').read_bytes()).hexdigest()=='78992e81d64fe2cf3dc911b07ccadf95'
 (DEST/'assets.json').write_text(json.dumps(entries,indent=2)+'\n')
 print('Public assets downloaded and hashed; panel index matches published MD5.')
