"""Plot reported imputation confidence against hidden-genotype recovery."""
from pathlib import Path
import argparse
import json
import os

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLBACKEND', 'Agg')
os.environ.setdefault('MPLCONFIGDIR', str(ROOT/'data/01-processed/private/matplotlib-cache'))
import matplotlib.pyplot as plt


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--figures-only',action='store_true',help='Generate figures without writing any notebook')
    parser.add_argument('--results-dir',type=Path,default=ROOT/'results/private/04_imputation')
    parser.add_argument('--run-dir',type=Path,default=ROOT/'data/04-imputation/autosomal-20261003')
    args=parser.parse_args()
    review = args.run_dir/'review'
    summary = json.loads((review/'review.private.json').read_text(encoding='utf-8'))
    metrics = json.loads((review/'heldout_markers.private.json').read_text(encoding='utf-8'))
    folder = args.results_dir/'figures'
    folder.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), sharex=True, sharey=True)
    for ax, breed in zip(axes, ['GSP', 'LR']):
        for chrom, color in zip([1,3,13,20,36], ['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00']):
            rows = [m for m in metrics if m['breed']==breed and m['chromosome']==f'chr{chrom}' and m['dosage_r2'] is not None]
            ax.scatter([m['reported_DR2'] for m in rows], [m['dosage_r2'] for m in rows], s=19, alpha=.6, color=color, label=f'Chr {chrom}')
        ax.plot([0,1],[0,1], color='gray', linewidth=1, linestyle=':')
        ax.set(title='German Shorthaired Pointer' if breed=='GSP' else 'Labrador Retriever', xlabel='Beagle estimated accuracy (DR2)', xlim=(-.02,1.02), ylim=(-.02,1.02))
        ax.grid(alpha=.15); ax.spines[['top','right']].set_visible(False)
    axes[0].set_ylabel('Measured dosage accuracy (correlation squared)')
    axes[1].legend(frameon=False, loc='upper left', bbox_to_anchor=(1.02, 1))
    fig.suptitle('Does reported confidence match recovery of hidden genotypes?')
    fig.tight_layout()
    png = folder/'imputation_confidence.png'
    fig.savefig(png,dpi=180); fig.savefig(folder/'imputation_confidence.svg'); plt.close(fig)
    print('Saved confidence figures; notebook sources unchanged.')


if __name__ == '__main__':
    main()
