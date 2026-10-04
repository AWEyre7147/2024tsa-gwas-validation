"""Small, dependency-free helpers shared by the executable publication chapters."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys


def study_path(root, relative):
    return root / str(relative).replace('\\', '/')


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8'))


def file_digest(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def verify_datasets(root, manifest):
    for dataset in manifest['datasets']:
        for relative, expected in dataset['files'].items():
            path = study_path(root, relative)
            if file_digest(path) != expected:
                raise ValueError(f'Input file differs from its recorded version: {relative}')
    return len(manifest['datasets'])


def run_script(root, name, *arguments, linux=False):
    if linux and os.name == 'nt':
        raise RuntimeError('Run this analysis cell with the WSL/Linux project Python kernel.')
    subprocess.run([sys.executable, str(root/'scripts'/name), *map(str, arguments)],
                   cwd=root, check=True)


def plink_environment(root):
    """Locate the pinned native binary without importing analysis dependencies."""
    folder = 'plink-1.9.0-20260927-win64' if os.name == 'nt' else 'plink-1.9.0-20260927-linux64'
    binary = root/'data/tools'/folder/('plink.exe' if os.name == 'nt' else 'plink')
    if not binary.exists():
        run_script(root, 'setup_plink.py')
    version = subprocess.check_output([str(binary), '--version'], text=True)
    assert 'PLINK v1.9.0' in version
    path = root/'data/processed/private'/f'notebook-plink-{os.name}.private.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(dict(executable=str(binary), executable_sha256=file_digest(binary)), indent=2), encoding='utf-8')
    return path
