"""Cache the complete public Dog10K panel and verify its published MD5."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess

ROOT=Path(__file__).resolve().parents[1]
URL='https://kiddlabshare.med.umich.edu/dog10K/phased-imputation-panel/AutoAndXPAR.Dog10K.phased.bcf'
EXPECTED_MD5='b88ac98f30c32b938b54410633fe94e8'


def checks(path):
    md5=hashlib.md5();sha=hashlib.sha256()
    with path.open('rb') as f:
        while b:=f.read(1048576):md5.update(b);sha.update(b)
    assert md5.hexdigest()==EXPECTED_MD5,'Published reference-panel MD5 mismatch'
    return dict(md5=md5.hexdigest(),sha256=sha.hexdigest(),bytes=path.stat().st_size)


def ensure_panel():
    folder=ROOT/'data/reference/dog10k';folder.mkdir(parents=True,exist_ok=True)
    path=folder/'AutoAndXPAR.Dog10K.phased.bcf'
    if not path.exists():
        partial=path.with_name(path.name+'.partial');curl=shutil.which('curl.exe') or shutil.which('curl')
        assert curl,'curl is required to download the public reference panel'
        subprocess.run([curl,'--fail','--location','--silent','--show-error','--retry','2','--continue-at','-',URL,'--output',str(partial)],check=True)
        result=checks(partial);partial.replace(path)
    else:result=checks(path)
    index=path.with_name(path.name+'.csi')
    original=ROOT/'data/reference/imputation-pilot/AutoAndXPAR.Dog10K.phased.bcf.csi'
    assert hashlib.md5(original.read_bytes()).hexdigest()=='78992e81d64fe2cf3dc911b07ccadf95'
    if not index.exists():shutil.copyfile(original,index)
    assert index.read_bytes()==original.read_bytes()
    result.update(url=URL,published_md5_verified=True)
    (folder/'panel_assets.json').write_text(json.dumps(result,indent=2)+'\n')
    return path


if __name__=='__main__':
    print(ensure_panel())
