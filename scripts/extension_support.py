"""Provenance helpers for new study extensions; existing frozen modules stay unchanged."""
from pathlib import Path
import csv
import datetime
import gzip
import json
import platform
import subprocess
from notebook_support import file_digest, read_json, study_path

ROOT = Path(__file__).resolve().parents[1]


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def rows(path, whitespace=False):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8', newline='') as handle:
        if whitespace:
            data = [line.split() for line in handle if line.strip()]
            return [dict(zip(data[0], r)) for r in data[1:]]
        return list(csv.DictReader(handle, delimiter='\t'))


def tsv(path, records, fields=None):
    fields = fields or list(records[0])
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'wt', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader(); writer.writerows(records)


def sources(paths):
    return {str(p.resolve().relative_to(ROOT)): file_digest(p) for p in paths}


def start(folder, plan):
    receipt = folder / 'receipt.private.json'
    if receipt.exists():
        saved = read_json(receipt)
        if saved['plan'] != plan:
            raise ValueError('Completed extension specification changed; choose a new output directory')
        for name, digest in saved['outputs'].items():
            if file_digest(folder / name) != digest:
                raise ValueError(f'Completed output differs: {name}')
        print('Verified completed extension:', folder.name, flush=True)
        return False
    if folder.exists():
        raise FileExistsError('Incomplete/existing output folder; choose a new directory')
    folder.mkdir(parents=True)
    save(folder / 'plan.private.json', plan)
    return True


def finish(folder, plan, summary):
    for relative, digest in plan['sources'].items():
        if file_digest(study_path(ROOT, relative)) != digest:
            raise ValueError('Source changed during extension execution')
    save(folder / 'summary.private.json', summary)
    outputs = {str(p.relative_to(folder)): file_digest(p) for p in folder.rglob('*') if p.is_file()}
    save(folder / 'receipt.private.json', dict(plan=plan, outputs=outputs,
        completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), python=platform.python_version()))


def command(argv, log):
    proc = subprocess.run(list(map(str, argv)), capture_output=True, text=True)
    log.write_text(proc.stdout + proc.stderr, encoding='utf-8')
    proc.check_returncode()
    return dict(argv=list(map(str, argv)), returncode=proc.returncode)
