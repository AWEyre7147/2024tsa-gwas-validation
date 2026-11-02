"""Plot descriptive provider gene overlap; do not rank association or causal evidence."""
from pathlib import Path
from collections import defaultdict
import argparse
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from notebook_support import read_json, file_digest

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--annotation', type=Path)
    args = parser.parse_args()
    config = read_json(ROOT / 'config/annotation.json')
    folder = args.annotation or ROOT / config['output_directory']
    receipt = read_json(folder / 'receipt.private.json')
    rows = {}
    for name in ['variant_crosswalk.private.tsv', 'gene_annotations.private.tsv']:
        path = folder / name
        if file_digest(path) != receipt['outputs'][name]:
            raise ValueError('Annotation output changed')
        with path.open(newline='') as handle:
            rows[name] = list(csv.DictReader(handle, delimiter='\t'))
    genes = defaultdict(list)
    for r in rows['gene_annotations.private.tsv']:
        if r['relation'] == 'overlap':
            genes[r['record_id'], r['provider']].append(r['biotype'])
    labels = ['Coding only', 'Noncoding / pseudogene only', 'Both', 'No gene overlap', 'Mapping unresolved']
    colors = ['#33658a', '#7a5195', '#2a9d8f', '#c8d1da', '#e9a23b']
    categories = ['historical_target', 'breed_specific_candidate']
    fig, ax = plt.subplots(figsize=(10, 4.7))
    bars = []
    for group in categories:
        variants = [r for r in rows['variant_crosswalk.private.tsv'] if r['evidence_group'] == group]
        for provider in ['Ensembl', 'RefSeq']:
            counts = [0] * len(labels)
            for r in variants:
                types = genes[r['record_id'], provider]
                if r['mapping_status'] != 'validated_coordinate':
                    index = 4
                elif not types:
                    index = 3
                elif 'protein_coding' in types and any(t != 'protein_coding' for t in types):
                    index = 2
                else:
                    index = 0 if 'protein_coding' in types else 1
                counts[index] += 1
            bars.append((('Historical targets' if group == 'historical_target' else 'Breed-specific candidates') + '\n' + provider, counts))
    for i, (label, color) in enumerate(zip(labels, colors)):
        left = [sum(counts[:i]) for _, counts in bars]
        values = [counts[i] for _, counts in bars]
        ax.barh(range(len(bars)), values, left=left, color=color, label=label, edgecolor='white', linewidth=.7)
        for y, width, offset in zip(range(len(bars)), values, left):
            if width:
                ax.text(offset+width/2, y, str(width), ha='center', va='center', color='white' if i < 3 else '#243447', fontsize=10)
    ax.set_yticks(range(len(bars)), [label for label, _ in bars])
    ax.invert_yaxis()
    ax.set_xlabel('Evidence records (one count per variant and provider)')
    ax.set_title('ROS gene overlap differs between annotation providers', loc='left', pad=16)
    ax.set_xlim(0, 25)
    ax.set_xticks(range(0, 26, 5))
    ax.spines[['top', 'right', 'left']].set_visible(False)
    ax.legend(loc='upper center', bbox_to_anchor=(.5, -.20), ncol=3, frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(folder / 'gene_overlap_summary.png', dpi=180, bbox_inches='tight')
    plt.close(fig)
    print('Saved descriptive annotation figure.')


if __name__ == '__main__':
    main()
