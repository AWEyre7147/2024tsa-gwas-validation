"""Independently reconstruct non-PAR X heterozygosity and add descriptive Y counts.

Conventional X sex-screen discordance is distinguished from ambiguous calls.
Y array call rates are supporting observations, not independently validated sex
assignments. The completed diagnostic receipt is not changed.
"""
from pathlib import Path
import argparse
import numpy as np
from extension_support import ROOT, rows, tsv, save
from notebook_support import read_json,file_digest


def raw_calls(prefix):
    fam=[r.split() for r in prefix.with_suffix('.fam').read_text().splitlines()]
    bim=[r.split() for r in prefix.with_suffix('.bim').read_text().splitlines()]
    binary=prefix.with_suffix('.bed').read_bytes();stride=(len(fam)+3)//4
    if binary[:3]!=bytes([108,27,1]) or len(binary)!=3+len(bim)*stride:raise ValueError('Unexpected PLINK format')
    return fam,bim,np.frombuffer(binary,dtype=np.uint8,offset=3).reshape(len(bim),stride)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--diagnostics',type=Path)
    args=parser.parse_args();config=read_json(ROOT/'config/diagnostics.json');folder=args.diagnostics or ROOT/config['output_directory']
    receipt=read_json(folder/'receipt.private.json')
    for name,digest in receipt['outputs'].items():
        if file_digest(folder/name)!=digest:raise ValueError('Diagnostic output receipt differs')
    expected={r['study_id']:r for r in rows(folder/'sex_checks.private.tsv')};records=[]
    for breed in ['LR','GSP']:
        fam,bim,packed=raw_calls(ROOT/'data/01-input/not-behelim-vs-behelim'/breed/('2024TSA_'+breed))
        marker_ids={r.strip() for r in (folder/breed/'X_prune.prune.in').read_text().splitlines()}
        x=[i for i,r in enumerate(bim) if r[1] in marker_ids]
        y=[i for i,r in enumerate(bim) if r[0]=='40' and int(r[3])>0]
        freq={r['SNP']:float(r['MAF']) for r in rows(folder/breed/'X_freq.frq',whitespace=True)}
        eh=np.array([2*freq[bim[i][1]]*(1-freq[bim[i][1]]) for i in x])
        for j,dog in enumerate(fam):
            codes=(packed[x,j//4]>>np.uint8(2*(j%4)))&3
            called=codes!=1;het=codes==2
            F=float(1-het.sum()/eh[called].sum())
            native=expected[dog[1]]
            if int(called.sum())!=int(native['X_called']) or int(het.sum())!=int(native['X_heterozygous_calls']) or abs(F-float(native['X_F']))>.001:
                raise ValueError('Independent X heterozygosity/F check disagrees with native result')
            screen=int(native['screening_sex']);recorded=int(dog[4])
            status='ambiguous' if screen==0 else 'concordant' if screen==recorded else 'discordant_screen'
            yc=(packed[y,j//4]>>np.uint8(2*(j%4)))&3
            records.append(dict(breed=breed,study_id=dog[1],recorded_sex=recorded,screening_sex=screen,
                screening_status=status,X_F_independent=F,X_heterozygosity_rate=float(het.sum()/called.sum()),
                Y_markers=len(y),Y_called=int(np.sum(yc!=1)),Y_call_fraction=float(np.mean(yc!=1)),
                Y_heterozygous_calls=int(np.sum(yc==2)),interpretation='Array Y call counts are descriptive; no reassignment'))
    tsv(folder/'sex_support.private.tsv',records)
    summary=dict(status='passed',X_F_checks=len(records),sex_screen_categories={name:sum(r['screening_status']==name for r in records)
        for name in ['concordant','ambiguous','discordant_screen']},
        diagnostic_receipt_sha256=file_digest(folder/'receipt.private.json'),verifier_sha256=file_digest(Path(__file__)),
        no_sex_or_outcome_changes=True)
    save(folder/'verification.private.json',summary)
    print(summary)


if __name__=='__main__':main()
