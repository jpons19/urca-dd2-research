#!/usr/bin/env python3
"""Reproducible DD2 fixed-coupling pion-screening scan and reference figures.

From repository root: python -m calculations.pion_softening.scan_dd2
Outputs in figures/pion_softening; never alters baseline CSV/plots.
"""
import argparse
import csv
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import qmc

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'calculations'))
import compare_compose_four as base
from .bottaro_softened import softened_rates


def run(npoints=36, power=14, ratios=(1.0,.75,.5,.25), output=None):
    if npoints<2 or power<8 or power>20:
        raise ValueError('npoints>=2, 8<=power<=20 are required')
    ratios=tuple(float(r) for r in ratios)
    output=Path(output) if output else ROOT/'figures'/'pion_softening'
    output.mkdir(parents=True,exist_ok=True)
    samples=qmc.Sobol(d=4,scramble=True,seed=42).random_base2(power)
    _,cold,micro,*_=base.setup()
    rows=[]
    for nb in np.linspace(.16,.98,npoints):
        d=base.make_data(float(nb),cold,micro)
        vacuum=softened_rates(d,1.,base,samples)
        for ratio in ratios:
            vals=softened_rates(d,ratio,base,samples)
            row={'nB_fm3':float(nb),'m_screen_over_m_pi':ratio,'pFn_MeV':d['pn'],
                 'Yp':d['yp'],'DU_gap_e_MeV':d['gap_e'],'DU_gap_mu_MeV':d['gap_mu']}
            for i,key in enumerate(('tri_e','tri_mu','nontri_e','nontri_mu')):
                row['Q_'+key]=vals[i]
                row['R_'+key]=vals[i]/vacuum[i] if vacuum[i]>0 else np.nan
            for name,start in [('tri',0),('nontri',2)]:
                row['Q_'+name+'_total']=sum(vals[start:start+2]);row['R_'+name+'_total']=sum(vals[start:start+2])/sum(vacuum[start:start+2])
            rows.append(row)
    path=output/'DD2_Bottaro_pion_screening.csv'
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    fig,ax=plt.subplots(figsize=(9,5.4))
    for ratio in ratios:
        sel=[r for r in rows if r['m_screen_over_m_pi']==ratio]
        if ratio==1.:
            ax.axhline(1,color='.4',ls=':',label='Unmodified OPE');continue
        line=ax.plot([r['nB_fm3'] for r in sel],[r['R_nontri_total'] for r in sel],lw=2,label=f'm_scr/m_pi={ratio:g}')[0]
        ax.plot([r['nB_fm3'] for r in sel],[r['R_tri_total'] for r in sel],color=line.get_color(),lw=1.3,ls='--')
    ax.set(xlabel=r'Baryon density $n_B$ (fm$^{-3}$)',ylabel='Ratio to same-method vacuum-pion result',title='PT(GRDF2-DD2): static pion screening sensitivity at $10^9$ K')
    ax.legend(title='Solid: nontriangular; dashed: triangular',fontsize=9);ax.grid(alpha=.22)
    fig.text(.10,.015,'DD2 composition and effective masses, NN vertex, rho, Landau couplings and beta=0.68 held fixed. NOT collective MMU.',fontsize=8,color='.3')
    fig.subplots_adjust(left=.12,right=.98,bottom=.18,top=.9)
    fig.savefig(output/'DD2_Bottaro_pion_screening.png',dpi=180)
    fig.savefig(output/'DD2_Bottaro_pion_screening.pdf')
    plt.close(fig)
    return rows

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--npoints',type=int,default=36)
    p.add_argument('--sobol-power',type=int,default=14)
    p.add_argument('--ratios',type=float,nargs='+',default=[1.,.75,.5,.25])
    p.add_argument('--output',default=None)
    a=p.parse_args()
    rows=run(a.npoints,a.sobol_power,a.ratios,a.output)
    print('Generated',len(rows),'density/mass samples')
    for density in (.16,.5,.9):
        rr=min([r for r in rows if r['m_screen_over_m_pi']==.5],key=lambda r:abs(r['nB_fm3']-density)) if .5 in a.ratios else None
        if rr:print(f"nB={rr['nB_fm3']:.3f}: R_tri={rr['R_tri_total']:.4f}, R_nontri={rr['R_nontri_total']:.4f}")
