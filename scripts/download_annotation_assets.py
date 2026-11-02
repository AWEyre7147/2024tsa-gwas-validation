"""Cache release-pinned ROS annotations and UCSC conversion assets without replacement.

Run with the Linux notebook environment. Published MD5 (NCBI/UCSC) and
BSD sum (Ensembl) checks are verified; SHA256 locks every cached resource.
"""
from pathlib import Path
import argparse
import concurrent.futures
import datetime
import gzip
import hashlib
import json
import shutil
import subprocess
from notebook_support import file_digest, read_json

ROOT = Path(__file__).resolve().parents[1]
HUB = 'https://hgdownload.soe.ucsc.edu/hubs/GCF/014/441/545/GCF_014441545.1/'
NCBI = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9615/106/GCF_014441545.1_ROS_Cfam_1.0/'
ENSEMBL = 'https://ftp.ensembl.org/pub/release-116/gtf/canis_lupus_familiaris/'
ASSETS = {
    'ros.2bit': HUB + 'GCF_014441545.1.2bit',
    'ucsc.md5sum.txt': HUB + 'md5sum.txt',
    'refseq.gff.gz': NCBI + 'GCF_014441545.1_ROS_Cfam_1.0_genomic.gff.gz',
    'assembly_report.txt': NCBI + 'GCF_014441545.1_ROS_Cfam_1.0_assembly_report.txt',
    'ncbi.md5checksums.txt': NCBI + 'md5checksums.txt',
    'ensembl.gtf.gz': ENSEMBL + 'Canis_lupus_familiaris.ROS_Cfam_1.0.116.gtf.gz',
    'ensembl.CHECKSUMS': ENSEMBL + 'CHECKSUMS',
    'canFam3ToROS.chain.gz': 'https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/canFam3ToGCF_014441545.1.over.chain.gz',
    'canFam4ToROS.chain.gz': 'https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/liftOver/canFam4ToGCF_014441545.1.over.chain.gz',
    'ROSToCanFam3.chain.gz': HUB + 'liftOver/GCF_014441545.1ToCanFam3.over.chain.gz',
    'ROSToCanFam4.chain.gz': HUB + 'liftOver/GCF_014441545.1ToCanFam4.over.chain.gz',
}


def published_md5(path, listing, upstream_name):
    matches = [line.split()[0] for line in listing.read_text().splitlines()
               if line.split() and Path(line.split()[-1].lstrip('*')).name == upstream_name]
    if len(matches) != 1:
        raise ValueError(f'Missing/ambiguous published checksum for {upstream_name}')
    with path.open('rb') as handle:
        actual = hashlib.file_digest(handle, 'md5').hexdigest()
    if actual != matches[0]:
        raise ValueError(f'Published checksum mismatch: {path.name}')


def download(destination):
    destination.mkdir(parents=True, exist_ok=True)
    lock = destination / 'assets.json'
    if lock.exists():
        saved = read_json(lock)
        if saved['urls'] != ASSETS:
            raise ValueError('Resource URLs changed; use a new cache directory')
        for name, record in saved['files'].items():
            if file_digest(destination / name) != record['sha256']:
                raise ValueError(f'Cached resource differs: {name}')
        print('Verified completed ROS resource cache.', flush=True)
        return saved
    curl = shutil.which('curl') or shutil.which('curl.exe')
    if not curl:
        raise RuntimeError('curl is required')

    def fetch(item):
        name, url = item
        path = destination / name
        if not path.exists():
            partial = path.with_suffix(path.suffix + '.partial')
            subprocess.run([curl, '--fail', '--location', '--silent', '--show-error',
                            '--retry', '3', '--connect-timeout', '30', '--max-time', '900',
                            url, '--output', str(partial)], check=True)
            partial.replace(path)
        if name.endswith('.gz'):
            with gzip.open(path, 'rb') as handle:
                while handle.read(1024 * 1024):
                    pass
        print(f'Cached {name}', flush=True)
        return name, {'sha256': file_digest(path), 'bytes': path.stat().st_size}

    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        files = dict(pool.map(fetch, ASSETS.items()))
    published_md5(destination / 'ros.2bit', destination / 'ucsc.md5sum.txt', 'GCF_014441545.1.2bit')
    for local, upstream in [('refseq.gff.gz', 'GCF_014441545.1_ROS_Cfam_1.0_genomic.gff.gz'),
                            ('assembly_report.txt', 'GCF_014441545.1_ROS_Cfam_1.0_assembly_report.txt')]:
        published_md5(destination / local, destination / 'ncbi.md5checksums.txt', upstream)
    line = next(line for line in (destination / 'ensembl.CHECKSUMS').read_text().splitlines()
                if line.split()[-1] == 'Canis_lupus_familiaris.ROS_Cfam_1.0.116.gtf.gz')
    expected = [int(x) for x in line.split()[:2]]
    actual = [int(x) for x in subprocess.check_output(['sum', '-r', str(destination / 'ensembl.gtf.gz')], text=True).split()[:2]]
    if actual != expected:
        raise ValueError('Ensembl published BSD checksum mismatch')
    report = {'status': 'complete', 'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'urls': ASSETS, 'files': files,
              'published_checks': ['UCSC genome MD5', 'NCBI GFF and assembly-report MD5', 'Ensembl GTF BSD sum'],
              'chain_checks': 'gzip integrity and SHA256; no published chain MD5 claimed'}
    lock.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path)
    args = parser.parse_args()
    config = read_json(ROOT / 'config/annotation.json')
    download(args.directory or ROOT / config['resource_directory'])
