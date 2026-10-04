"""Cache the study's pinned PLINK 1.9 executable for Windows or Linux."""
from pathlib import Path
import hashlib
import json
import platform
import subprocess
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parents[1]
BUILD = '20260927'


def main():
    system = platform.system()
    if system == 'Linux' and platform.machine().lower() in {'x86_64', 'amd64'}:
        name, asset, executable = 'linux64', 'linux_x86_64', 'plink'
    elif system == 'Windows' and platform.machine().lower() in {'amd64', 'x86_64'}:
        name, asset, executable = 'win64', 'win64', 'plink.exe'
    else:
        raise RuntimeError('This setup supports x86-64 Windows and Linux/WSL.')
    folder = ROOT / 'data/tools' / f'plink-1.9.0-{BUILD}-{name}'
    folder.mkdir(parents=True, exist_ok=True)
    archive = folder / 'download.zip'
    url = f'https://s3.amazonaws.com/plink1-assets/plink_{asset}_{BUILD}.zip'
    target = folder / executable
    receipt = folder / 'installation.private.json'
    previous = json.loads(receipt.read_text(encoding='utf-8')) if receipt.exists() else None
    if not target.exists():
        print('Downloading pinned PLINK 1.9 release...', flush=True)
        urllib.request.urlretrieve(url, archive)
        with zipfile.ZipFile(archive) as bundle:
            # Extract only the known executable, without arbitrary archive paths.
            target.write_bytes(bundle.read(executable))
    target.chmod(target.stat().st_mode | 0o111)
    digest = hashlib.sha256(target.read_bytes()).hexdigest()
    if previous:
        assert digest == previous['executable_sha256'], 'Cached executable differs from recorded installation'
    version = subprocess.check_output([str(target), '--version'], text=True).strip()
    assert 'PLINK v1.9.0' in version and '27 Sep 2026' in version
    details = dict(executable=str(target), executable_sha256=digest,
                   archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest() if archive.exists() else None,
                   source_url=url, version=version)
    receipt.write_text(json.dumps(details, indent=2)+'\n', encoding='utf-8')
    print(version)
    print('Installed:', target)


if __name__ == '__main__':
    main()
