"""Make validation figures and embed them in a readable local review notebook."""
from pathlib import Path
import argparse
import base64
import json
import os
import nbformat
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLBACKEND','Agg')
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'data/processed/private/matplotlib-cache'))
import matplotlib.pyplot as plt


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run-dir',type=Path,default=ROOT/'data/imputation/autosomal-20261003')
    ap.add_argument('--results-dir',type=Path,default=ROOT/'results/private/04_imputation')
    ap.add_argument('--allow-subset',action='store_true',help='Only for an archived implementation check')
    ap.add_argument('--figures-only',action='store_true',help='Generate figures without writing any notebook')
    args=ap.parse_args();m=json.loads((args.run_dir/'manifest.private.json').read_text())
    assert m['status']=='complete','Finish the selected run before summarizing it'
    if not args.allow_subset:assert set(m['jobs'])=={'chr'+str(c) for c in range(1,39)},'Not a completed whole-autosome run'
    dest=args.results_dir;figures=dest/'figures';figures.mkdir(parents=True,exist_ok=True)
    data={};table=[]
    for region,job in m['jobs'].items():
        key=region.replace(':','_').replace('-','_')
        for breed,record in job['breeds'].items():
            p=args.run_dir/'validation'/breed/key/'masked_variant_metrics.private.json'
            if not p.exists():continue
            data[breed,region]=json.loads(p.read_text());v=record['verification']
            table.append([breed,region,v['masked_markers'],v['called_truth_genotypes'],f"{100*v['hard_call_concordance']:.2f}%",f"{v['pooled_dosage_r2']:.3f}"])
    assert data,'No held-out genotype results found'
    regions=sorted({r for b,r in data},key=lambda r:int(r.split(':')[0][3:]))
    palette=['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00']
    colors={r:palette[i%len(palette)] for i,r in enumerate(regions)}
    fig,axes=plt.subplots(1,2,figsize=(11,4.6),sharex=True,sharey=True)
    for ax,breed in zip(axes,['GSP','LR']):
        for region in regions:
            rows=[r for r in data.get((breed,region),[]) if r['dosage_r2'] is not None]
            ax.scatter([r['maf'] for r in rows],[r['dosage_r2'] for r in rows],s=19,alpha=.65,color=colors[region],label=region.split(':')[0],edgecolors='none')
        ax.set(title='German Shorthaired Pointer' if breed=='GSP' else 'Labrador Retriever',xlabel='Measured minor allele frequency',xlim=(0,.51),ylim=(-.02,1.03))
        ax.grid(alpha=.15);ax.spines[['top','right']].set_visible(False)
    axes[0].set_ylabel('Measured dosage accuracy (correlation squared)')
    axes[1].legend(title='Validation region',frameon=False,loc='lower right')
    fig.suptitle('Recovery of array genotypes hidden before phasing')
    fig.text(.5,.015,'Each point is a held-out marker. Undefined correlations are omitted. These results assess array genotype recovery.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.06,1,.95))
    png=figures/'imputation_validation.png';fig.savefig(png,dpi=180);fig.savefig(figures/'imputation_validation.svg');plt.close(fig)
    if args.figures_only:
        print('Saved validation figures; notebook sources unchanged.')
        return
    text='## Recovery of hidden measurements\n\n| Breed | Chromosome / region | Hidden markers | Called measurements | Genotype agreement | Pooled dosage correlation squared |\n| --- | --- | ---: | ---: | ---: | ---: |\n'
    text+=''.join('| '+' | '.join(map(str,row))+' |\n' for row in table)
    text+='\nThese comparisons concern the predefined array-marker validation regions. They do not independently establish phase accuracy or reliability at every imputed sequence variant.\n'
    source=ROOT/'notebooks/05_imputation_validation.ipynb';n=nbformat.read(source,as_version=4)
    n.cells[-1].source=text
    cell=nbformat.v4.new_markdown_cell('![Recovery of hidden array genotypes](attachment:imputation_validation.png)')
    cell.attachments={'imputation_validation.png':{'image/png':base64.b64encode(png.read_bytes()).decode()}}
    n.cells.append(cell);n.metadata.setdefault('study',{})['status']='validation_results_available_association_review_pending'
    if args.allow_subset:
        n.cells.insert(0,nbformat.v4.new_markdown_cell('# Archived implementation check\n\nThis selected-region test is not the whole-genome study result.'))
    nbformat.validate(n);nbformat.write(n,dest/'05_imputation_validation.ipynb')
    (dest/'validation_summary.json').write_text(json.dumps(dict(regions=regions,table=table,subset=args.allow_subset),indent=2)+'\n')
    print('Saved validation figures and a notebook with the figure embedded.')


if __name__=='__main__':main()
