"""Run reproducible, breed-specific PLINK 1.9 array QC; retain private reports.

Python requires numpy. Supply a recorded PLINK environment with --environment.
Inputs are selected exclusively from data/01-input/manifest.private.json.
Existing QC outputs are never overwritten.
"""
from pathlib import Path
from collections import Counter
import argparse
import csv
import hashlib
import json
import re
import shutil
import subprocess
import numpy as np

SCENARIOS = [('a', 'success-vs-behelim'), ('b', 'not-behelim-vs-behelim')]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def study_path(root, relative):
    return root / relative.replace('\\', '/')

def lines(path):
    return [s.split() for s in path.read_text().splitlines() if s.strip()]

def table(path):
    rows = lines(path)
    return [dict(zip(rows[0], r)) for r in rows[1:]]

def tsv(path, headers, rows):
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f, delimiter='\t'); w.writerow(headers); w.writerows(rows)

def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2) + '\n', encoding='utf-8')

def verify_calls(source, output):
    """Check every retained call and allele against its prepared input."""
    sf, of = lines(source.with_suffix('.fam')), lines(output.with_suffix('.fam'))
    sb, ob = lines(source.with_suffix('.bim')), lines(output.with_suffix('.bim'))
    si = {r[1]: i for i, r in enumerate(sf)}
    vi = {r[1]: i for i, r in enumerate(sb)}
    assert len(si) == len(sf) and len(vi) == len(sb)
    inds = np.array([si[r[1]] for r in of])
    vind = np.array([vi[r[1]] for r in ob])
    assert all(r == sf[si[r[1]]] for r in of)
    assert all(r == sb[vi[r[1]]] for r in ob)
    def packed(prefix, samples, variants):
        data = prefix.with_suffix('.bed').read_bytes()
        stride = (samples + 3) // 4
        assert data[:3] == bytes([108, 27, 1])
        assert len(data) == 3 + variants * stride
        return np.frombuffer(data, dtype=np.uint8, offset=3).reshape(variants, stride)
    original = packed(source, len(sf), len(sb))[vind]
    expected = (original[:, inds // 4] >> (2 * (inds % 4)).astype(np.uint8)) & 3
    stored = packed(output, len(of), len(ob))
    js = np.arange(len(of))
    actual = (stored[:, js // 4] >> (2 * (js % 4)).astype(np.uint8)) & 3
    assert np.array_equal(expected, actual)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--environment', type=Path, required=True)
    parser.add_argument('--run-id', default='20261002-plink19-autosomal')
    parser.add_argument('--output',type=Path,help='Separate QC output directory inside the repository')
    args = parser.parse_args()
    root = args.root.resolve()
    if not re.fullmatch(r'[A-Za-z0-9_-]+', args.run_id):
        raise ValueError('Run ID must be a simple directory name')
    env = json.loads(args.environment.read_text())
    exe = Path(env['executable'])
    assert sha(exe) == env['executable_sha256']
    version = subprocess.run([str(exe), '--version'], check=True, capture_output=True, text=True).stdout.strip()
    assert 'v1.9.0' in version and '27 Sep 2026' in version
    input_manifest = root / 'data/01-input/manifest.private.json'
    inputs = json.loads(input_manifest.read_text())
    prefixes = {}
    tracked_hashes = dict(inputs['source_hashes'])
    for ds in inputs['datasets']:
        for p, h in ds['files'].items():
            assert sha(study_path(root,p)) == h
            tracked_hashes[p] = h
        fam_path = next(study_path(root,p) for p in ds['files'] if p.endswith('.fam'))
        prefixes[ds['scenario'], ds['breed']] = fam_path.with_suffix('')
    for p, h in tracked_hashes.items(): assert sha(study_path(root,p)) == h
    outroot = (args.output or root / 'data/02-qc').resolve()
    assert outroot.is_relative_to(root),'Keep QC outputs inside the repository'
    if outroot.exists(): raise FileExistsError('QC output exists; use a new directory for a new run')
    reports = outroot / 'reports' / args.run_id
    reports.mkdir(parents=True)
    policy = dict(plink_version=version, input_manifest_sha256=sha(input_manifest),
                  sample_missingness_max=0.10, marker_missingness_max=0.10, maf_min=0.025,
                  hwe_p_min=1e-5, hwe_test='standard exact, no mid-p',
                  missingness_and_maf_population='retained Scheme B dogs within breed',
                  hwe_population='retained Scheme A controls within breed',
                  structural_scope='autosomes 1-38, positive positions, unique IDs/coordinates, A/C/G/T SNP alleles',
                  order=['structural eligibility', 'sample missingness', 'marker missingness', 'MAF', 'HWE'],
                  comparison_policy='same breed-specific marker set; A subset of B',
                  assembly='array assembly not independently established; no coordinate conversion',
                  script_sha256=sha(Path(__file__)))
    write_json(reports/'policy.private.json', policy)
    commands = []
    def run(label, source, output, flags):
        cmd = [str(exe), '--bfile', str(source), '--dog', '--keep-allele-order',
               '--threads', '1', '--memory', '2048', *map(str,flags), '--out', str(output)]
        completed = subprocess.run(cmd, capture_output=True, text=True)
        output.with_suffix('.console.txt').write_text(completed.stdout + completed.stderr, encoding='utf-8')
        commands.append(dict(step=label, argv=cmd, returncode=completed.returncode))
        write_json(reports/'commands.private.json', commands)
        if completed.returncode: raise RuntimeError(f'PLINK failed: {label}; see private reports')
    datasets = []
    flows = []
    targets = []
    target_text = (root/'docs/references/2013TSA Significant Hits.md').read_text()
    target_coords = list(dict.fromkeys(re.findall(r'chr(\d+)\.(\d+)', target_text)))
    warnings = []
    for breed in ['GSP', 'LR']:
        print(f'{breed}: checking markers and missing data', flush=True)
        folder = reports/breed
        folder.mkdir()
        a = prefixes[SCENARIOS[0][1], breed]
        b = prefixes[SCENARIOS[1][1], breed]
        bim = lines(b.with_suffix('.bim'))
        assert lines(a.with_suffix('.bim')) == bim
        fam_a, fam_b = lines(a.with_suffix('.fam')), lines(b.with_suffix('.fam'))
        assert {r[1] for r in fam_a} <= {r[1] for r in fam_b}
        assert {r[1] for r in fam_a if r[5]=='2'} == {r[1] for r in fam_b if r[5]=='2'}
        ids = Counter(r[1] for r in bim)
        coords = Counter((r[0], r[3]) for r in bim if int(r[3]) > 0)
        exclusions = []
        eligible = []
        for row in bim:
            reasons = []
            if not (row[0].isdigit() and 1 <= int(row[0]) <= 38): reasons.append('nonautosomal_or_unplaced')
            if int(row[3]) <= 0: reasons.append('nonpositive_position')
            if ids[row[1]] > 1: reasons.append('duplicate_id')
            if int(row[3]) > 0 and coords[row[0],row[3]] > 1: reasons.append('repeated_coordinate')
            if any(x not in ['A','C','G','T'] for x in row[4:6]): reasons.append('not_ACGT_SNP')
            if reasons: exclusions.append([row[1],row[0],row[3],';'.join(reasons)])
            else: eligible.append(row[1])
        tsv(folder/'structural_exclusions.tsv', ['SNP','CHR','BP','reasons'], exclusions)
        extract = folder/'eligible_markers.txt'
        extract.write_text('\n'.join(eligible)+'\n')
        structural = folder/'01_structural'
        run(f'{breed}: structural', b, structural, ['--extract',extract,'--make-bed','--missing'])
        imiss = table(structural.with_suffix('.imiss'))
        tsv(folder/'sample_missingness_exclusions.tsv', ['FID','IID','F_MISS'],
            [[d['FID'],d['IID'],d['F_MISS']] for d in imiss if float(d['F_MISS'])>0.10])
        samples = folder/'02_samples'
        run(f'{breed}: sample missingness', structural,samples,['--mind','0.10','--make-bed','--missing'])
        markers = folder/'03_marker_missingness'
        run(f'{breed}: marker missingness', samples,markers,['--geno','0.10','--make-bed','--missing'])
        maf = folder/'04_maf'
        run(f'{breed}: MAF',markers,maf,['--maf','0.025','--make-bed','--freq'])
        retained_samples = {r[1] for r in lines(maf.with_suffix('.fam'))}
        controls = [r[:2] for r in fam_a if r[5]=='1' and r[1] in retained_samples]
        assert controls
        keep_controls = folder/'primary_controls.keep'
        keep_controls.write_text(''.join(' '.join(r)+'\n' for r in controls))
        hardy = folder/'05_primary_control_hwe'
        run(f'{breed}: primary-control HWE',maf,hardy,['--keep',keep_controls,'--hardy'])
        hwe_rows = [d for d in table(hardy.with_suffix('.hwe')) if d['TEST']=='ALL']
        assert len(hwe_rows) == len(lines(maf.with_suffix('.bim')))
        assert all(d['P'] not in ['NA','nan'] for d in hwe_rows)
        hwe_fail = [d['SNP'] for d in hwe_rows if float(d['P']) < 1e-5]
        fail_path = folder/'hwe_exclusions.txt'; fail_path.write_text('\n'.join(hwe_fail)+ ('\n' if hwe_fail else ''))
        final_b = folder/'06_final_b'
        run(f'{breed}: final B',maf,final_b,(['--exclude',fail_path] if hwe_fail else []) + ['--make-bed','--freq','--missing'])
        kept_b = {r[1] for r in lines(final_b.with_suffix('.fam'))}
        keep_a = folder/'primary_samples.keep'
        keep_a.write_text(''.join(' '.join(r[:2])+'\n' for r in fam_a if r[1] in kept_b))
        final_a = folder/'07_final_a'
        run(f'{breed}: final A', final_b,final_a,['--keep',keep_a,'--make-bed','--missing','--freq'])
        assert lines(final_a.with_suffix('.bim')) == lines(final_b.with_suffix('.bim'))
        assert all(float(d['F_MISS']) <= 0.10 for d in table(final_b.with_suffix('.lmiss')))
        assert all(float(d['F_MISS']) <= 0.10 for d in table(final_b.with_suffix('.imiss')))
        assert all(float(d['MAF']) >= 0.025 for d in table(final_b.with_suffix('.frq')))
        final_ids = {r[1] for r in lines(final_b.with_suffix('.bim'))}
        assert not final_ids.intersection(hwe_fail)
        source_by_coordinate = {(r[0],r[3]):r[1] for r in bim}
        for chrom,pos in target_coords:
            snp = source_by_coordinate.get((chrom,pos))
            targets.append([breed,f'chr{chrom}.{pos}',snp or '', 'retained' if snp in final_ids else 'excluded_by_QC' if snp else 'no_input_coordinate_match'])
        stages = [('input',b),('structural',structural),('sample_missingness',samples),('marker_missingness',markers),('MAF',maf),('HWE',final_b)]
        prior_ids = {r[1] for r in bim}
        for stage,prefix in stages:
            rows = lines(prefix.with_suffix('.fam')); variant_ids = {r[1] for r in lines(prefix.with_suffix('.bim'))}
            flows.append([breed,stage,len(rows),sum(r[5]=='1' for r in rows),sum(r[5]=='2' for r in rows),len(variant_ids)])
            if stage!='input':
                (folder/f'excluded_at_{stage}.txt').write_text('\n'.join(sorted(prior_ids-variant_ids))+'\n')
            prior_ids = variant_ids
        for (scheme,scenario),source in zip(SCENARIOS,[final_a,final_b]):
            output = outroot/scenario/breed/f'2024TSA_{breed}'
            output.parent.mkdir(parents=True)
            for ext in ['bed','bim','fam']: shutil.copy2(source.with_suffix('.'+ext),output.with_suffix('.'+ext))
            verify_calls(prefixes[scenario,breed],output)
            rows = lines(output.with_suffix('.fam'))
            datasets.append(dict(scenario=scenario,breed=breed,dogs=len(rows),controls=sum(r[5]=='1' for r in rows),cases=sum(r[5]=='2' for r in rows),variants=len(final_ids),all_retained_calls_verified=True,
                files={str(output.with_suffix('.'+e).relative_to(root)):sha(output.with_suffix('.'+e)) for e in ['bed','bim','fam']}))
        initial_max = max(float(d['F_MISS']) for d in imiss)
        print(f'{breed}: {len(final_ids)} markers retained; maximum initial dog missingness {initial_max:.4%}',flush=True)
    for p,h in tracked_hashes.items(): assert sha(study_path(root,p))==h
    assert sha(input_manifest)==policy['input_manifest_sha256']
    for log in reports.rglob('*.log'):
        for line in log.read_text().splitlines():
            if 'Warning' in line: warnings.append(dict(log=str(log.relative_to(root)),message=line))
    tsv(reports/'sample_variant_flow.tsv',['breed','stage','dogs','controls','cases','variants'],flows)
    tsv(reports/'replication_target_coordinate_review.tsv',['breed','historical_target','input_marker','QC_status'],targets)
    manifest = dict(status='complete',run_id=args.run_id,policy=policy,datasets=datasets,source_and_input_hashes_unchanged=True,
                    report_directory=str(reports.relative_to(root)),warnings=warnings,
                    limitation='autosomal SNP baseline; no sex inference, PCA, kinship, imputation, or association; historical target coordinates require assembly confirmation')
    write_json(outroot/'manifest.private.json',manifest)
    print(json.dumps(dict(datasets=[{k:d[k] for k in ['scenario','breed','dogs','controls','cases','variants']} for d in datasets],warnings=warnings),indent=2))

if __name__=='__main__': main()
