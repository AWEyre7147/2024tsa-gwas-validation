"""Native PLINK heterozygosity, non-PAR X sex and case/control missingness diagnostics.

Produces review flags only. Uses broader within-breed frequency estimates and
pruned markers; no impute-sex, filtering of study dogs, or association tests.
"""
from pathlib import Path
from collections import Counter
import argparse
import csv
import math
import numpy as np
from scipy.stats import fisher_exact
from breed_gwas import decode
from notebook_support import read_json, file_digest, study_path
from extension_support import ROOT, save, sources, start, finish, command, rows, tsv


def get_prefix(manifest, breed, scenario):
    ds = next(d for d in manifest['datasets'] if d['breed'] == breed and d['scenario'] == scenario)
    for relative, expected in ds['files'].items():
        if Path(relative).suffix in {'.bed', '.bim', '.fam'} and file_digest(study_path(ROOT, relative)) != expected:
            raise ValueError('Prepared/aligned input differs from accepted manifest')
    return next(study_path(ROOT, p) for p in ds['files'] if p.endswith('.fam')).with_suffix('')


def eligible(bim):
    counts = Counter((r[0], r[3]) for r in bim)
    ids = Counter(r[1] for r in bim)
    return [i for i, r in enumerate(bim) if r[0].isdigit() and 1 <= int(r[0]) <= 38 and int(r[3]) > 0
            and counts[r[0], r[3]] == 1 and ids[r[1]] == 1 and r[4] in 'ACGT' and r[5] in 'ACGT'
            and len(r[4]) == len(r[5]) == 1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--plink', type=Path, default=ROOT / 'data/02-tools/plink-1.9.0-20260927-linux64/plink')
    args = parser.parse_args(); config = read_json(ROOT / 'config/diagnostics.json')
    prepared = read_json(ROOT / 'data/01-input/manifest.private.json')
    aligned = read_json(ROOT / 'data/02-aligned/manifest.private.json')
    prefixes = {(kind, breed, scheme): get_prefix(manifest, breed, scheme)
                for kind, manifest in [('prepared', prepared), ('aligned', aligned)]
                for breed in ['LR', 'GSP'] for scheme in ['success-vs-behelim', 'not-behelim-vs-behelim']}
    backgrounds = {b: ROOT / 'data/03-structure/20261002-pca-umap-tsne' / b / '02_pruned_genotypes' for b in ['LR', 'GSP']}
    paths = [ROOT / 'config/diagnostics.json', Path(__file__), ROOT / 'scripts/extension_support.py',
             ROOT / 'scripts/breed_gwas.py', ROOT / 'scripts/replication_models.py', ROOT / 'scripts/notebook_support.py',
             ROOT / 'data/01-input/manifest.private.json', ROOT / 'data/02-aligned/manifest.private.json']
    paths += [p.with_suffix(s) for p in list(prefixes.values()) + list(backgrounds.values()) for s in ['.bed', '.bim', '.fam']]
    version = __import__('subprocess').check_output([str(args.plink), '--version'], text=True).strip()
    if '27 Sep 2026' not in version:
        raise ValueError('Use the pinned PLINK release')
    plan = dict(config=config, sources=sources(paths), PLINK=version, plink_sha256=file_digest(args.plink),
                numpy=np.__version__, scipy=__import__('scipy').__version__)
    output = args.output or ROOT / config['output_directory']
    if not start(output, plan):
        return
    commands = []; heter_records = []; sex_records = []; missing_summary = []
    for breed in ['LR', 'GSP']:
        folder = output / breed; folder.mkdir()
        bg = backgrounds[breed]
        base = [args.plink, '--bfile', bg, '--dog', '--keep-allele-order']
        for name, extra in [('frequencies', ['--freq']), ('heterozygosity', ['--read-freq', folder / 'frequencies.frq', '--het'])]:
            commands.append(command(base + extra + ['--out', folder / name], folder / (name + '.console.txt')))
        fam, variants, g = decode(bg)
        broadfam = [line.split() for line in prefixes['prepared', breed, 'not-behelim-vs-behelim'].with_suffix('.fam').read_text().splitlines()]
        if fam != broadfam:
            raise ValueError('Pruned background dog/metadata order differs from broader input')
        freq = {r['SNP']: float(r['MAF']) for r in rows(folder / 'frequencies.frq', whitespace=True)}
        expected_het = np.array([2*freq[r[1]]*(1-freq[r[1]]) for r in variants])
        native = rows(folder / 'heterozygosity.het', whitespace=True)
        rates = np.nanmean(g == 1, axis=0) * len(g) / np.isfinite(g).sum(axis=0)
        mean, sd = float(rates.mean()), float(rates.std(ddof=1))
        primary_ids = {line.split()[1] for line in prefixes['prepared', breed, 'success-vs-behelim'].with_suffix('.fam').read_text().splitlines()}
        for j, (dog, r) in enumerate(zip(fam, native)):
            observed = int(np.sum(g[:, j] == 1)); called = int(np.isfinite(g[:, j]).sum())
            eh = float(expected_het[np.isfinite(g[:, j])].sum()); F = 1-observed/eh
            if r['IID'] != dog[1] or int(r['N(NM)']) != called or int(r['O(HOM)']) != called-observed:
                raise ValueError('Independent heterozygosity count check failed')
            if abs(float(r['F'])-F) > .0001:
                raise ValueError('Independent F-coefficient check failed')
            z = (rates[j]-mean)/sd if sd else 0.
            heter_records.append(dict(breed=breed, study_id=dog[1], sex=int(dog[4]), phenotype=int(dog[5]),
                in_scheme_A=dog[1] in primary_ids, called_markers=called, observed_heterozygotes=observed,
                heterozygosity_rate=float(rates[j]), F=F, z_within_breed=float(z),
                review_flag=bool(abs(z)>config['heterozygosity_z_flag'])))
        prefix = prefixes['prepared', breed, 'not-behelim-vs-behelim']
        xf, xb, xg = decode(prefix)
        par = [r for r in xb if r[0] == config['PAR_code']]
        nonpar = [r for r in xb if r[0] == config['X_nonPAR_code']]
        if not par or not nonpar or max(int(r[3]) for r in par) >= min(int(r[3]) for r in nonpar):
            raise ValueError('Prepared canine X/PAR encoding is not separated as expected')
        xmarkers = [r[1] for r in nonpar if int(r[3]) > config['conservative_X_min_bp']
                    and len(r[4]) == len(r[5]) == 1 and r[4] in 'ACGT' and r[5] in 'ACGT']
        (folder / 'X_nonPAR_markers.txt').write_text('\n'.join(xmarkers)+'\n')
        xbase = [args.plink, '--bfile', prefix, '--dog', '--keep-allele-order', '--extract', folder / 'X_nonPAR_markers.txt',
                 '--maf', config['X_maf_min'], '--geno', config['X_missingness_max']]
        commands.append(command(xbase + ['--indep-pairwise', *config['X_pruning'], '--out', folder / 'X_prune'], folder / 'X_prune.console.txt'))
        prune = folder / 'X_prune.prune.in'
        if not prune.read_text().strip():
            raise ValueError('No usable pruned non-PAR X markers')
        xbase = [args.plink, '--bfile', prefix, '--dog', '--keep-allele-order', '--extract', prune]
        commands.append(command(xbase + ['--freq', '--out', folder / 'X_freq'], folder / 'X_freq.console.txt'))
        commands.append(command(xbase + ['--read-freq', folder / 'X_freq.frq', '--check-sex', config['female_max_F'],
                         config['male_min_F'], '--out', folder / 'sex_check'], folder / 'sex_check.console.txt'))
        sx = {r['IID']: r for r in rows(folder / 'sex_check.sexcheck', whitespace=True)}
        selected = {line.strip() for line in prune.read_text().splitlines()}
        ix = [i for i, r in enumerate(xb) if r[1] in selected]
        for j, dog in enumerate(xf):
            r = sx[dog[1]]
            if int(r['PEDSEX']) != int(dog[4]):
                raise ValueError('Sex diagnostic changed the recorded sex')
            sex_records.append(dict(breed=breed, study_id=dog[1], recorded_sex=int(dog[4]),
                in_scheme_A=dog[1] in primary_ids, X_F=float(r['F']), screening_sex=int(r['SNPSEX']),
                PLINK_status=r['STATUS'], X_markers=len(ix), X_called=int(np.isfinite(xg[ix,j]).sum()),
                X_heterozygous_calls=int(np.sum(xg[ix,j] == 1)),
                review_flag=r['STATUS'] != 'OK', interpretation='Conventional screen; no sex reassignment'))
        for scheme in ['success-vs-behelim', 'not-behelim-vs-behelim']:
            for label, kind in [('prepared_structurally_eligible', 'prepared'), ('aligned_analysis', 'aligned')]:
                prefix = prefixes[kind, breed, scheme]
                dogs, markers, g = decode(prefix); ix = eligible(markers)
                case = np.array([r[5] == '2' for r in dogs]); nc, nu = int(case.sum()), int((~case).sum())
                tag = scheme + '_' + label
                extract = folder / (tag + '.markers.txt')
                extract.write_text(''.join(markers[i][1]+'\n' for i in ix))
                native_prefix = folder / tag
                commands.append(command([args.plink, '--bfile', prefix, '--dog', '--keep-allele-order', '--extract', extract,
                    '--test-missing', '--out', native_prefix], folder / (tag + '.console.txt')))
                native = {r['SNP']: r for r in rows(native_prefix.with_suffix('.missing'), whitespace=True)}
                counts = Counter(); minimum = 1.; records = []; checked = 0
                cutoff = config['missingness_alpha']/len(ix)
                for i in ix:
                    missing = np.isnan(g[i]); a, u = int(missing[case].sum()), int(missing[~case].sum())
                    key = (a,u)
                    if key not in counts:
                        counts[key] = float(fisher_exact([[a,nc-a],[u,nu-u]]).pvalue)
                    p = counts[key] if a+u < len(dogs) else None
                    if p is not None:minimum = min(minimum,p)
                    r = native.get(markers[i][1])
                    if r is not None:
                        checked += 1
                        # PLINK prints fractions to four significant digits;
                        # compare at the precision available in the report.
                        if abs(float(r['P'])-p) > max(1e-6,p*.001) or abs(float(r['F_MISS_A'])-a/nc) > .000051 or abs(float(r['F_MISS_U'])-u/nu) > .000051:
                            raise ValueError(f'Independent missingness/Fisher check failed: {markers[i][1]}, {r}, {a}/{nc}, {u}/{nu}, p={p}')
                    elif (a or u) and a+u < len(dogs):
                        raise ValueError('Missingness native output absent for a marker with missing calls')
                    records.append(dict(marker=markers[i][1], chromosome=markers[i][0], position=int(markers[i][3]),
                        case_missing=a, case_called=nc-a, control_missing=u, control_called=nu-u,
                        case_missing_fraction=a/nc, control_missing_fraction=u/nu, fisher_p=p,
                        diagnostic_bonferroni=cutoff, review_flag=p is not None and p<cutoff,
                        test_status='all_calls_missing_untestable' if p is None else 'no_missing_calls' if a+u==0 else 'tested'))
                tsv(folder / (tag + '.private.tsv.gz'), records)
                flagged = sum(r['review_flag'] for r in records)
                missing_summary.append(dict(breed=breed, scheme=scheme, family=label, markers=len(ix), cases=nc, controls=nu,
                    markers_with_missing=sum(r['case_missing']+r['control_missing']>0 for r in records),
                    minimum_p=minimum, cutoff=cutoff, flagged=flagged, native_rows_checked=checked,
                    all_missing_untestable=sum(r['fisher_p'] is None for r in records)))
                print(breed, scheme, label, 'diagnostic flags:', flagged, flush=True)
    tsv(output / 'heterozygosity.private.tsv', heter_records)
    tsv(output / 'sex_checks.private.tsv', sex_records)
    tsv(output / 'missingness_summary.private.tsv', missing_summary)
    save(output / 'commands.private.json', commands)
    finish(output, plan, dict(status='diagnostics_complete', heterozygosity_dogs=len(heter_records),
        heterozygosity_flags=sum(r['review_flag'] for r in heter_records), sex_dogs=len(sex_records),
        sex_flags=sum(r['review_flag'] for r in sex_records), missingness=missing_summary,
        sources_unchanged=True, exclusions_applied=False,
        limitations=['F thresholds are conventional screening thresholds, not validated sex assignments.',
                     'Within-breed z flags are descriptive; relatives and inbreeding can affect heterozygosity.',
                     'Fisher missingness tests assume independent dogs; relatedness and sex can affect interpretation.',
                     'No data changes or automatic exclusions.']))


if __name__ == '__main__':
    main()
