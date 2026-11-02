"""Render QC and pilot power figures from recorded extension results."""
from pathlib import Path
import argparse
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from extension_support import ROOT,rows
from notebook_support import read_json


def plot_qc():
    config=read_json(ROOT/'config/diagnostics.json');qc=ROOT/config['output_directory']
    sex=rows(qc/'sex_support.private.tsv');heter=rows(qc/'heterozygosity.private.tsv')
    fig,axes=plt.subplots(1,2,figsize=(11,4.5))
    for ax,breed in zip(axes,['LR','GSP']):
        for code,label,color in [('1','Recorded male','#33658a'),('2','Recorded female','#a65f8e')]:
            values=[float(r['X_F_independent']) for r in sex if r['breed']==breed and r['recorded_sex']==code]
            order=np.arange(len(values))
            ax.scatter(values,order,s=22,label=label,color=color,alpha=.85)
        ax.axvspan(.2,.8,color='#e9a23b',alpha=.12)
        ax.axvline(.2,color='#999999',linestyle='--',linewidth=1)
        ax.axvline(.8,color='#999999',linestyle='--',linewidth=1)
        ax.set_title(breed+' non-PAR X sex screen');ax.set_xlabel('X inbreeding coefficient F');ax.set_ylabel('Index within recorded-sex group')
        ax.legend(frameon=False,fontsize=8);ax.spines[['top','right']].set_visible(False)
    fig.tight_layout();fig.savefig(qc/'sex_screen.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(9,4))
    for ax,breed in zip(axes,['LR','GSP']):
        subset=[r for r in heter if r['breed']==breed];z=[float(r['z_within_breed']) for r in subset]
        ax.hist(z,bins=12,color='#33658a',edgecolor='white');ax.axvline(-3,color='#a65f8e',linestyle='--');ax.axvline(3,color='#a65f8e',linestyle='--')
        ax.set_title(breed);ax.set_xlabel('Within-breed heterozygosity z score');ax.set_ylabel('Dogs')
        ax.spines[['top','right']].set_visible(False)
    fig.tight_layout();fig.savefig(qc/'heterozygosity.png',dpi=180);plt.close(fig)


def plot_power():
    power=ROOT/read_json(ROOT/'config/power.json')['output_directory'];rates=rows(power/'rates.private.tsv')
    contexts=read_json(power/'contexts.private.json')
    fig,axes=plt.subplots(1,len(contexts),figsize=(13,4.8),sharey=True,squeeze=False)
    colors={'GMMAT_score':'#33658a','GEMMA_score':'#2a9d8f','GEMMA_wald':'#a65f8e','GEMMA_lrt':'#e9a23b'}
    for ax,context in zip(axes[0],contexts):
        for method,color in colors.items():
            subset=sorted([r for r in rates if r['context']==context['context'] and r['test']==method
                and r['variance_multiplier']=='1.0' and r['threshold_label']=='analysis_family' and r['valid']!='0'],key=lambda r:float(r['generating_OR']))
            x=[float(r['generating_OR']) for r in subset];y=[float(r['rate_valid']) for r in subset]
            lower=[v-float(r['Wilson_low']) for r,v in zip(subset,y)];upper=[float(r['Wilson_high'])-v for r,v in zip(subset,y)]
            ax.errorbar(x,y,yerr=[lower,upper],label=method.replace('_',' '),color=color,marker='o',markersize=4,capsize=2,linewidth=1.2)
        ax.set_title(context['breed']+' '+context['family']+'\nMAF '+format(context['MAF'],'.3f'),fontsize=10)
        ax.set_xlabel('Generating per-ALT odds ratio');ax.spines[['top','right']].set_visible(False);ax.set_ylim(-.03,1.03)
    axes[0,0].set_ylabel('Rejection fraction among valid fits')
    handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=4,frameon=False,fontsize=9)
    fig.suptitle('Initial finite-cohort simulation pilot: 40 replicates per scenario\nFamily thresholds; fitted background variance; Wilson 95% intervals',fontsize=11)
    fig.tight_layout(rect=(0,.08,1,.90));fig.savefig(power/'pilot_power.png',dpi=180);plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=['qc','power','all'],default='all')
    args=parser.parse_args()
    if args.stage in {'qc','all'}:plot_qc()
    if args.stage in {'power','all'}:plot_power()
    print('Saved requested study-extension figures.')


if __name__=='__main__':main()
