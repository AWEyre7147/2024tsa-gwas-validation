"""Map the agreed measured/historical inventory to ROS and annotate both providers.

This is a descriptive, SNP-only interval annotation workflow. It preserves source
coordinates, alleles and group membership. It does not run association tests,
infer historical effects, predict amino-acid substitutions, or alter genotypes.
"""
from pathlib import Path
from collections import Counter, defaultdict
import argparse
import csv
import datetime
import gzip
import json
import math
import re
import shutil
import subprocess
from urllib.parse import unquote
from notebook_support import file_digest, read_json, study_path
from reference_alignment import TwoBit, context, COMP
from investigate_array_manifests import parse as parse_manifest

ROOT = Path(__file__).resolve().parents[1]


def table(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_table(path, records, fields=None):
    fields = fields or list(records[0])
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(records)


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def bim(path):
    rows = [line.split() for line in path.read_text().splitlines() if line.strip()]
    result = {r[1]: r for r in rows}
    if len(rows) != len(result):
        raise ValueError(f'Duplicate BIM markers: {path}')
    return result


def verified_bim(manifest, breed, sources):
    dataset = next(d for d in manifest['datasets'] if d['breed'] == breed and d['scenario'] == 'success-vs-behelim')
    relative = next(p for p in dataset['files'] if p.endswith('.bim'))
    path = study_path(ROOT, relative)
    digest = file_digest(path)
    if digest != dataset['files'][relative]:
        raise ValueError('BIM differs from completed upstream manifest')
    sources[str(path.relative_to(ROOT))] = digest
    return bim(path)


def inventory(config, sources):
    def read_source(relative):
        path = ROOT / relative
        sources[relative] = file_digest(path)
        return path
    aligned_manifest = read_json(read_source('data/02-aligned/manifest.private.json'))
    original_manifest = read_json(read_source('data/01-input/manifest.private.json'))
    if aligned_manifest['target_assembly'] != 'CanFam4 / UU_Cfam_GSD_1.0':
        raise ValueError('Unexpected analysis assembly')
    aligned = {breed: verified_bim(aligned_manifest, breed, sources) for breed in ['LR', 'GSP']}
    original = verified_bim(original_manifest, 'LR', sources)
    historical = re.findall(r'\|\s*(chr\d+\.\d+)\s*\|', read_source(config['historical_source']).read_text())
    available = {r['historical_target']: r for r in table(read_source(config['historical_availability'])) if r['breed'] == 'LR'}
    if len(historical) != 18 or len(set(historical)) != 18 or set(historical) != set(available):
        raise ValueError('Historical inventory differs from the agreed 18 targets')
    records = []
    for target in historical:
        a = available[target]
        chromosome, position = target[3:].split('.')
        marker = a['input_marker']
        source = original[marker]
        if (source[0], source[3]) != (chromosome, position):
            raise ValueError('Historical coordinate disagrees with prepared input')
        current = aligned['LR'].get(marker)
        records.append(dict(record_id='historical:' + target, evidence_group='historical_target', breed='LR',
                            marker=marker, historical_target=target, QC_status=a['QC_status'],
                            current_tested=a['QC_status'] == 'retained', source_assembly='CanFam3.1',
                            source_chromosome='chr' + chromosome, source_position=int(position),
                            source_A1=source[4], source_A2=source[5],
                            canFam4_chromosome='chr' + current[0] if current else '',
                            canFam4_position=int(current[3]) if current else '',
                            canFam4_REF=current[5] if current else '', canFam4_ALT=current[4] if current else '',
                            selection_tests='', display_region=''))
    review = read_json(read_source(config['candidate_review']))
    policy = read_json(read_source('config/gwas_reporting.json'))
    summary = read_json(read_source(config['scan_results'] + '/summary.private.json'))
    if summary['reporting_policy_sha256'] != sources['config/gwas_reporting.json']:
        raise ValueError('GWAS summary uses a different reporting policy')
    model_path = config['candidate_review'].rsplit('/', 1)[0] + '/candidate_models.private.tsv'
    models = {(r['region'], r['marker']): r for r in table(read_source(model_path)) if r['scheme'] == 'A'}
    for breed in ['LR', 'GSP']:
        scan = table(read_source(config['scan_results'] + '/' + breed + '_success-vs-behelim.private.tsv.gz'))
        threshold = policy['family_alpha'] / policy['candidate_markers'][breed]
        selected = {r['marker']: r for r in scan if any(math.isfinite(float(r[k])) and float(r[k]) < threshold
                                                      for k in ['GEMMA_wald_p', 'GEMMA_lrt_p'])}
        reviewed = {marker: region['region'] for region in review['regions'] if region['breed'] == breed
                    for marker in region['selected_markers']}
        if set(selected) != set(reviewed) or len(selected) != {'LR': 3, 'GSP': 21}[breed]:
            raise ValueError('Primary candidate inventory differs from Chapter 09')
        for marker, row in sorted(selected.items()):
            current = aligned[breed][marker]
            if (int(current[0]), int(current[3]), current[5], current[4]) != (
                    int(row['chromosome']), int(row['position']), row['REF'], row['ALT']):
                raise ValueError('Candidate representation disagrees with aligned BIM')
            if (reviewed[marker], marker) not in models:
                raise ValueError('Candidate missing from model review')
            tests = [label for label, key in [('Wald', 'GEMMA_wald_p'), ('likelihood_ratio', 'GEMMA_lrt_p')]
                     if math.isfinite(float(row[key])) and float(row[key]) < threshold]
            records.append(dict(record_id=breed + ':' + marker, evidence_group='breed_specific_candidate', breed=breed,
                                marker=marker, historical_target='', QC_status='retained', current_tested=True,
                                source_assembly='CanFam4', source_chromosome='chr' + current[0], source_position=int(current[3]),
                                source_A1=current[4], source_A2=current[5], canFam4_chromosome='chr' + current[0],
                                canFam4_position=int(current[3]), canFam4_REF=current[5], canFam4_ALT=current[4],
                                selection_tests=';'.join(tests), display_region=reviewed[marker]))
    if len(records) != 42:
        raise ValueError('Expected 18 historical and 24 breed-specific evidence records')
    return records


def liftover(sites, chain, folder, name, executable, config, commands):
    source = folder / (name + '.input.bed')
    mapped = folder / (name + '.mapped.bed')
    unmapped = folder / (name + '.unmapped.bed')
    source.write_text(''.join(f'{c}\t{p-1}\t{p}\t{key}\t0\t+\n' for key, c, p in sites), encoding='utf-8')
    command = [executable, '-multiple', '-noSerial', f'-minMatch={config["liftover_min_match"]}',
               str(source), str(chain), str(mapped), str(unmapped)]
    proc = subprocess.run(command, capture_output=True, text=True)
    commands.append(dict(argv=command, returncode=proc.returncode, stdout=proc.stdout, stderr=proc.stderr))
    proc.check_returncode()
    hits = defaultdict(list)
    for line in mapped.read_text().splitlines():
        r = line.split()
        hits[r[3]].append((r[0], int(r[1]) + 1, int(r[2]), r[5]))
    return dict(hits)


def assembly_aliases(path, genome):
    aliases = {}
    evidence = []
    for line in path.read_text().splitlines():
        if line.startswith('#') or not line.strip():
            continue
        r = line.split('\t')
        if r[1] != 'assembled-molecule' or r[2] not in {str(i) for i in range(1, 39)}:
            continue
        if r[5] != '=':
            raise ValueError('GenBank and RefSeq chromosome sequences are not identical')
        choices = [r[6], r[4], r[0], 'chr' + r[2], r[9]]
        keys = [key for key in choices if key in genome.index]
        if len(keys) != 1 or genome.index[keys[0]][0] != int(r[8]):
            raise ValueError('ROS genome chromosome/length disagrees with assembly report')
        for name in set(choices) - {'na'}:
            aliases[name] = keys[0]
        evidence.append(dict(chromosome=r[2], genbank=r[4], refseq=r[6], relationship=r[5],
                             length=int(r[8]), genome_name=keys[0]))
    if len(evidence) != 38:
        raise ValueError('Expected all 38 identical GenBank/RefSeq autosomes')
    return aliases, evidence


def flank_check(source, target, chrom, pos, target_chrom, target_pos, strand, config):
    n = config['context_flank_bp']
    s = source.get(chrom, pos-1-n, pos+n)
    t = target.get(target_chrom, target_pos-1-n, target_pos+n)
    if len(s) != 2*n+1 or len(t) != 2*n+1:
        return 0, 0, False
    if strand == '-':
        s = s.translate(COMP)[::-1]
    pairs = [(a, b) for i, (a, b) in enumerate(zip(s, t)) if i != n and a in 'ACGT' and b in 'ACGT']
    mismatches = sum(a != b for a, b in pairs)
    return len(pairs), mismatches, len(pairs) >= config['minimum_context_bases'] and mismatches <= config['maximum_context_mismatches']


def map_records(records, resources, output, executable, config):
    commands = []
    reports = output / 'mapping'
    reports.mkdir()
    cf3 = TwoBit(ROOT / 'data/02-reference/genomes/canFam3.2bit')
    cf4 = TwoBit(ROOT / 'data/02-reference/genomes/canFam4.2bit')
    ros = TwoBit(resources / 'ros.2bit')
    aliases, assembly_evidence = assembly_aliases(resources / 'assembly_report.txt', ros)
    save_json(output / 'assembly_compatibility.private.json', assembly_evidence)
    _, _, manifest, _ = parse_manifest(ROOT / 'data/02-reference/array-manifests/caninehd_b.csv')
    direct = {}
    for assembly in ['CanFam3.1', 'CanFam4']:
        label = 'canFam3' if assembly == 'CanFam3.1' else 'canFam4'
        sites = [(r['record_id'], r['source_chromosome'], r['source_position']) for r in records if r['source_assembly'] == assembly]
        direct.update(liftover(sites, resources / (label + 'ToROS.chain.gz'), reports, label, executable, config, commands))
    alternate_sites = [(r['record_id'], r['canFam4_chromosome'], r['canFam4_position']) for r in records
                       if r['evidence_group'] == 'historical_target' and r['canFam4_position']]
    alternate = liftover(alternate_sites, resources / 'canFam4ToROS.chain.gz', reports, 'historical_via_canFam4', executable, config, commands)
    reverse = {}
    for assembly in ['CanFam3.1', 'CanFam4']:
        sites = [(r['record_id'], direct[r['record_id']][0][0], direct[r['record_id']][0][1]) for r in records
                 if r['source_assembly'] == assembly and len(direct.get(r['record_id'], [])) == 1]
        label = 'canFam3' if assembly == 'CanFam3.1' else 'canFam4'
        reverse.update(liftover(sites, resources / ('ROSTo' + label[0].upper() + label[1:] + '.chain.gz'), reports,
                                'reverse_' + label, executable, config, commands))
    for r in records:
        key = r['record_id']
        hits = direct.get(key, [])
        r.update(mapping_count=len(hits), mapping_status='unmapped', ROS_chromosome='', ROS_position='',
                 chain_strand='', source_reference='', ROS_reference='', ROS_other_allele='', allele_status='not_checked',
                 context_bases=0, context_mismatches=0, reverse_status='not_checked', historical_route_status='not_applicable')
        if len(hits) != 1:
            r['mapping_status'] = 'multiple_mappings' if hits else 'unmapped'
            continue
        c, p, end, strand = hits[0]
        r.update(ROS_chromosome=c, ROS_position=p, chain_strand=strand)
        if c not in set(aliases.values()) or end != p:
            r['mapping_status'] = 'unsupported_destination'
            continue
        rev = reverse.get(key, [])
        expected = (r['source_chromosome'], r['source_position'], r['source_position'], strand)
        r['reverse_status'] = 'unique_exact_return' if rev == [expected] else 'reverse_mismatch_or_ambiguity'
        genome = cf3 if r['source_assembly'] == 'CanFam3.1' else cf4
        bases, mismatches, passed = flank_check(genome, ros, r['source_chromosome'], r['source_position'], c, p, strand, config)
        r.update(context_bases=bases, context_mismatches=mismatches)
        r['source_reference'] = genome.get(r['source_chromosome'], r['source_position']-1, r['source_position'])
        r['ROS_reference'] = ros.get(c, p-1, p)
        if r['evidence_group'] == 'historical_target' and r['canFam4_position']:
            other = alternate.get(key, [])
            r['historical_route_status'] = 'same_ROS_coordinate' if len(other) == 1 and other[0][:3] == hits[0][:3] else 'route_disagreement_or_ambiguity'
        if r['evidence_group'] == 'historical_target' and not r['canFam4_position']:
            r['historical_route_status'] = 'no_retained_canFam4_representation'
        if r['reverse_status'] != 'unique_exact_return':
            r['mapping_status'] = 'reverse_check_failed'
        elif not passed:
            r['mapping_status'] = 'sequence_context_failed'
        elif r['historical_route_status'] == 'route_disagreement_or_ambiguity':
            r['mapping_status'] = 'historical_route_failed'
        else:
            r['mapping_status'] = 'validated_coordinate'
        pair = None
        if r['source_assembly'] == 'CanFam4':
            pair = {r['canFam4_REF'], r['canFam4_ALT']}
            if r['source_reference'] != r['canFam4_REF']:
                raise ValueError('Analytical REF disagrees with CanFam4 sequence')
        else:
            d = manifest.get(r['marker'])
            orientation, _ = context(genome, r['source_chromosome'], r['source_position'], d['TopGenomicSeq']) if d else (None, {})
            if orientation:
                pair = {r['source_A1'], r['source_A2']}
                if orientation == '-':
                    pair = {a.translate(COMP) for a in pair}
        if pair is None or len(pair) != 2 or not pair <= set('ACGT') or r['source_reference'] not in pair:
            r['allele_status'] = 'source_allele_orientation_unresolved'
        else:
            converted = {a.translate(COMP) for a in pair} if strand == '-' else pair
            if r['ROS_reference'] not in converted:
                r['allele_status'] = 'destination_reference_not_in_allele_pair'
                r['mapping_status'] = 'allele_check_failed'
            else:
                r['ROS_other_allele'] = next(iter(converted - {r['ROS_reference']}))
                oriented_ref = r['source_reference'].translate(COMP) if strand == '-' else r['source_reference']
                r['allele_status'] = 'reference_preserved' if oriented_ref == r['ROS_reference'] else 'reference_switched'
    save_json(output / 'liftover_commands.private.json', commands)
    return aliases


def attributes(text, provider):
    if provider == 'Ensembl':
        return dict(re.findall(r'(\S+) "([^"]*)";', text))
    return {key: unquote(value) for item in text.split(';') if '=' in item for key, value in [item.split('=', 1)]}


def feature_rows(path, provider, aliases):
    with gzip.open(path, 'rt', encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            row = line.rstrip('\n').split('\t')
            if len(row) != 9 or row[0] not in aliases:
                continue
            yield aliases[row[0]], row[2], int(row[3]), int(row[4]), row[6], attributes(row[8], provider)


def gene_inventory(path, provider, aliases, config):
    genes = {}
    for chrom, kind, start, end, strand, a in feature_rows(path, provider, aliases):
        if kind not in {'gene', 'pseudogene'}:
            continue
        key = a['gene_id'] if provider == 'Ensembl' else a['ID']
        dbxref = a.get('Dbxref', '')
        entrez = next((x.split(':', 1)[1] for x in dbxref.split(',') if x.startswith('GeneID:')), '')
        genes[key] = dict(provider=provider, release=config['ensembl_release'] if provider == 'Ensembl' else config['refseq_release'],
                          assembly_accession=config['ensembl_assembly'] if provider == 'Ensembl' else config['refseq_assembly'],
                          gene_id=key, gene_version=a.get('gene_version', ''), entrezgene=entrez,
                          symbol=a.get('gene_name', a.get('gene', a.get('Name', ''))),
                          biotype=a.get('gene_biotype', 'pseudogene' if kind == 'pseudogene' else ''),
                          chromosome=chrom, start=start, end=end, strand=strand)
    return genes


def distance(pos, gene):
    return max(gene['start'] - pos, pos - gene['end'], 0)


def transcript_classification(pos, transcript):
    """Classify closed genomic intervals without asserting a functional effect."""
    cds = any(a <= pos <= b for a, b in transcript['CDS'])
    exon = any(a <= pos <= b for a, b in transcript['exon'])
    utr = any(a <= pos <= b for a, b in transcript['UTR'])
    if cds:
        return 'CDS_overlap'
    if utr or (exon and transcript['CDS']):
        return 'UTR_overlap'
    if exon:
        return 'noncoding_transcript_exon' if not transcript['CDS'] else 'exon_overlap'
    exons = sorted(set(transcript['exon']))
    for (_, left_end), (right_start, _) in zip(exons, exons[1:]):
        if left_end < pos < right_start and (pos <= left_end+2 or pos >= right_start-2):
            return 'splice_boundary_2bp'
    return 'intron_overlap'


def annotate(records, resources, aliases, config, output):
    accepted = [r for r in records if r['mapping_status'] == 'validated_coordinate']
    gene_links = []
    transcripts_out = []
    gene_totals = {}
    gene_stats = {}
    for provider, filename in [('Ensembl', 'ensembl.gtf.gz'), ('RefSeq', 'refseq.gff.gz')]:
        path = resources / filename
        genes = gene_inventory(path, provider, aliases, config)
        gene_totals[provider] = len(genes)
        gene_stats[provider] = dict(Counter(g['biotype'] for g in genes.values()))
        by_chrom = defaultdict(list)
        for key, gene in genes.items():
            by_chrom[gene['chromosome']].append((key, gene))
        selected_genes = set()
        for r in accepted:
            pos = r['ROS_position']
            ranked = sorted(by_chrom[r['ROS_chromosome']], key=lambda item: (distance(pos, item[1]), item[1]['gene_id']))
            nonoverlap = [distance(pos, gene) for _, gene in ranked if distance(pos, gene) > 0]
            nearest_limit = nonoverlap[min(config['nearest_genes_per_provider'], len(nonoverlap))-1] if nonoverlap else -1
            for key, gene in ranked:
                d = distance(pos, gene)
                # The neighborhood is a descriptive ROS window, not a lifted LD region.
                if d > 500000 and d > nearest_limit:
                    continue
                relation = 'overlap' if d == 0 else 'within_500kb' if d <= 500000 else 'nearest_outside_window'
                gene_links.append(dict(record_id=r['record_id'], evidence_group=r['evidence_group'], breed=r['breed'],
                                       marker=r['marker'], ROS_position=pos, relation=relation, distance_bp=d,
                                       nearest_nonoverlap=d > 0 and d <= nearest_limit, **gene))
                if d == 0:
                    selected_genes.add(key)
        transcripts = {}
        parent_genes = {}
        for chrom, kind, start, end, strand, a in feature_rows(path, provider, aliases):
            if provider == 'Ensembl':
                gid = a.get('gene_id')
                tid = a.get('transcript_id')
                if gid not in selected_genes or not tid:
                    continue
                if tid not in transcripts:
                    transcripts[tid] = dict(gene_id=gid, chromosome=chrom, start=start, end=end, strand=strand,
                                            biotype=a.get('transcript_biotype', ''), version=a.get('transcript_version', ''),
                                            CDS=[], exon=[], UTR=[])
                t = transcripts[tid]
                t['start'] = min(t['start'], start)
                t['end'] = max(t['end'], end)
                if kind in {'CDS', 'exon'}:
                    t[kind].append((start, end))
                if kind.endswith('utr') or kind == 'UTR':
                    t['UTR'].append((start, end))
            else:
                parents = a.get('Parent', '').split(',')
                parent = next((p for p in parents if p in selected_genes), None)
                if parent and kind not in {'exon', 'CDS'} and a.get('ID'):
                    tid = a['ID']
                    parent_genes[tid] = parent
                    transcripts[tid] = dict(gene_id=parent, chromosome=chrom, start=start, end=end, strand=strand,
                                            biotype=kind, version='', CDS=[], exon=[], UTR=[])
                for tid in parents:
                    if tid in transcripts and kind in {'CDS', 'exon'}:
                        transcripts[tid][kind].append((start, end))
        for r in accepted:
            for tid, t in transcripts.items():
                if t['chromosome'] == r['ROS_chromosome'] and t['start'] <= r['ROS_position'] <= t['end']:
                    g = genes[t['gene_id']]
                    transcripts_out.append(dict(record_id=r['record_id'], provider=provider, release=g['release'],
                                                assembly_accession=g['assembly_accession'], gene_id=t['gene_id'], symbol=g['symbol'],
                                                gene_biotype=g['biotype'], transcript_id=tid, transcript_version=t['version'],
                                                transcript_biotype=t['biotype'], strand=t['strand'],
                                                region_annotation=transcript_classification(r['ROS_position'], t)))
    write_table(output / 'gene_annotations.private.tsv', gene_links,
                ['record_id', 'evidence_group', 'breed', 'marker', 'ROS_position', 'relation', 'distance_bp', 'nearest_nonoverlap',
                 'provider', 'release', 'assembly_accession', 'gene_id', 'gene_version', 'entrezgene', 'symbol', 'biotype',
                 'chromosome', 'start', 'end', 'strand'])
    write_table(output / 'transcript_annotations.private.tsv', transcripts_out,
                ['record_id', 'provider', 'release', 'assembly_accession', 'gene_id', 'symbol', 'gene_biotype',
                 'transcript_id', 'transcript_version', 'transcript_biotype', 'strand', 'region_annotation'])
    provider_comparison = []
    for r in records:
        overlaps = {provider: [g for g in gene_links if g['record_id'] == r['record_id'] and g['provider'] == provider and g['relation'] == 'overlap']
                    for provider in ['Ensembl', 'RefSeq']}
        symbols = {provider: {g['symbol'] for g in values if g['symbol']} for provider, values in overlaps.items()}
        provider_comparison.append(dict(record_id=r['record_id'], mapping_status=r['mapping_status'],
                                         Ensembl_overlap_ids=';'.join(g['gene_id'] for g in overlaps['Ensembl']),
                                         RefSeq_overlap_ids=';'.join(g['gene_id'] for g in overlaps['RefSeq']),
                                         shared_overlap_symbols=';'.join(sorted(symbols['Ensembl'] & symbols['RefSeq'])),
                                         Ensembl_only_overlap_symbols=';'.join(sorted(symbols['Ensembl'] - symbols['RefSeq'])),
                                         RefSeq_only_overlap_symbols=';'.join(sorted(symbols['RefSeq'] - symbols['Ensembl']))))
    write_table(output / 'provider_comparison.private.tsv', provider_comparison)
    return dict(autosomal_gene_inventory=gene_totals, autosomal_gene_biotypes=gene_stats,
                gene_links=len(gene_links), transcript_links=len(transcripts_out),
                interval_annotations=dict(Counter(t['region_annotation'] for t in transcripts_out)),
                noncoding_overlap_links=sum(g['relation'] == 'overlap' and g['biotype'] != 'protein_coding' for g in gene_links))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--liftover', type=Path)
    args = parser.parse_args()
    config = read_json(ROOT / 'config/annotation.json')
    resources = ROOT / config['resource_directory']
    output = (args.output or ROOT / config['output_directory']).resolve()
    sources = {}
    for relative in ['config/annotation.json', 'scripts/annotate_candidates.py', 'scripts/download_annotation_assets.py',
                     'scripts/reference_alignment.py', 'scripts/investigate_array_manifests.py', 'scripts/notebook_support.py',
                     'data/02-reference/genomes/canFam3.2bit', 'data/02-reference/genomes/canFam4.2bit',
                     'data/02-reference/array-manifests/caninehd_b.csv']:
        sources[relative] = file_digest(ROOT / relative)
    assets_path = resources / 'assets.json'
    assets = read_json(assets_path)
    sources[str(assets_path.relative_to(ROOT))] = file_digest(assets_path)
    for name, evidence in assets['files'].items():
        relative = str((resources / name).relative_to(ROOT))
        sources[relative] = file_digest(resources / name)
        if sources[relative] != evidence['sha256']:
            raise ValueError(f'Resource differs from cache receipt: {name}')
    records = inventory(config, sources)
    executable = str(args.liftover or shutil.which('liftOver') or Path.home() / 'anaconda3/envs/tsa2024-imputation/bin/liftOver')
    executable_hash = file_digest(Path(executable))
    receipt = output / 'receipt.private.json'
    if receipt.exists():
        saved = read_json(receipt)
        if saved['sources'] != sources or saved['liftover_sha256'] != executable_hash:
            raise ValueError('Completed annotation inputs/code changed; choose a new output directory')
        for name, digest in saved['outputs'].items():
            if file_digest(output / name) != digest:
                raise ValueError(f'Completed annotation output changed: {name}')
        print('Verified completed ROS annotation without recomputing.', flush=True)
        return
    if output.exists():
        raise FileExistsError('Incomplete/existing annotation output; choose a new directory')
    output.mkdir(parents=True)
    save_json(output / 'plan.private.json', dict(config=config, sources=sources, liftover_sha256=executable_hash,
                                               created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    aliases = map_records(records, resources, output, executable, config)
    write_table(output / 'variant_crosswalk.private.tsv', records)
    annotation = annotate(records, resources, aliases, config, output)
    summary = dict(status='descriptive_ROS_annotation_complete', evidence_records=len(records),
                   groups=dict(Counter(r['evidence_group'] for r in records)),
                   mapping_status=dict(Counter(r['mapping_status'] for r in records)),
                   allele_status=dict(Counter(r['allele_status'] for r in records)),
                   historical_route_status=dict(Counter(r['historical_route_status'] for r in records if r['historical_target'])),
                   annotation=annotation,
                   limitations=['Interval annotations are not amino-acid or regulatory-effect predictions.',
                                'Nearest genes and shared symbols do not prove a functional target or cross-provider equivalence.',
                                'Historical effect contrasts remain unreconciled; excluded targets are annotated without restoration.',
                                'No imputed or cross-breed candidates, new association tests or causal claims.',
                                'MyGene functional enrichment and comparison with the prior publication are separate follow-up work.'])
    save_json(output / 'summary.private.json', summary)
    outputs = {str(p.relative_to(output)): file_digest(p) for p in output.rglob('*') if p.is_file()}
    for relative, digest in sources.items():
        if file_digest(study_path(ROOT, relative)) != digest:
            raise ValueError('Source changed during annotation')
    save_json(receipt, dict(status='complete', sources=sources, liftover_sha256=executable_hash, outputs=outputs,
                           completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
