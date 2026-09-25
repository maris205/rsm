#!/usr/bin/env python3
"""Render archived matching diagnostics, without refitting or solving backgrounds."""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
COLORS={'global':'#6c7079','canonical':'#c16526','shift':'#007f8b'}
LABELS={'global':'Global reference','canonical':'Canonical K','shift':'Shift K'}


def main():
    out=ROOT/'figures';out.mkdir(exist_ok=True)
    with (ROOT/'results/matched_curves.csv').open() as f:
        rows=list(csv.DictReader(f))
    data={m:{k:np.array([float(r[k]) for r in rows if r['model']==m])
             for k in rows[0] if k!='model'} for m in COLORS}
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,
                         'axes.spines.right':False,'legend.frameon':False})
    def save(fig,name):
        for ext in ('pdf','png'):
            fig.savefig(out/(name+'.'+ext),dpi=200,bbox_inches='tight')
        plt.close(fig)
    fig,ax=plt.subplots(1,2,figsize=(11.4,4.3),layout='constrained')
    for m,r in data.items():
        ax[0].plot(r['z'],1e6*r['holomorphic_uv_delta_alpha'],color=COLORS[m],
                   label=LABELS[m]+' (fixed holomorphic UV)',
                   ls='--' if m=='shift' else '-',lw=4 if m=='global' else 2)
    r=data['canonical']
    ax[0].plot(r['z'],1e6*r['physical_uv_delta_alpha'],color=COLORS['canonical'],
               ls=':',lw=2.5,label='Canonical K (fixed physical UV)')
    ax[0].set(ylabel=r'$10^6[\alpha_U(z)/\alpha_U(0)-1]$',
              title='(a) UV boundary determines the comparison')
    handles,labels=ax[0].get_legend_handles_labels()
    fig.legend(handles,labels,fontsize=9,loc='outside lower center',ncol=2)
    for key,label,color,ls in (
        ('mass_delta_g_inverse2','Physical mass contribution','#c16526',':'),
        ('anomaly_delta_g_inverse2','KL anomaly contribution','#6650a4','--'),
        ('matched_delta_g_inverse2','Sum: fixed holomorphic UV','#007f8b','-')):
        ax[1].plot(r['z'],1e3*4*np.pi*r[key],label=label,color=color,ls=ls,lw=2)
    ax[1].set(ylabel=r'$10^3\,\Delta(\alpha_U^{-1})$',
              title='(b) Canonical K: terms in the inverse coupling')
    ax[1].legend(fontsize=9,loc='lower right')
    for a in ax:
        a.set(xlabel='Redshift z (evolution toward the right)',xlim=(4.2,0))
        a.axhline(0,color='#aaaaaa',lw=.8);a.grid(alpha=.15)
    fig.suptitle('Specified minimal U(1) | one-loop matching | no observational fit',fontsize=13)
    save(fig,'matched_gauge_response')

    fig,ax=plt.subplots(1,2,figsize=(11.4,4.3),layout='constrained')
    r=data['shift']
    for key,label,color,ls,lw in (
        ('shape_matched','Matched coupling','#007f8b','-',3.5),
        ('shape_chi',r'$\chi^{-2}$ response','#555555','--',1.5),
        ('shape_time',r'$\ln^{-2}(t/t_*)$ response','#c16526',':',2.5)):
        ax[0].plot(r['z'],r[key],label=label,color=color,ls=ls,lw=lw)
    ax[0].set(ylabel='Endpoint-normalized response',title='(a) Shift K: similar, not identical shapes')
    ax[0].legend(fontsize=9,loc='upper right')
    for m,r in data.items():
        ax[1].plot(r['z'],100*(r['shape_chi']-r['shape_time']),
                   label=LABELS[m],color=COLORS[m],
                   lw=3.5 if m=='global' else 2,ls='--' if m=='shift' else '-')
    ax[1].set(ylabel='Shape difference (percentage points)',
              title=r'(b) $\chi^{-2}$ minus $\ln^{-2}(t/t_*)$')
    ax[1].legend(fontsize=9,loc='upper right')
    for a in ax:
        a.set(xlabel='Redshift z (evolution toward the right)',xlim=(4.2,0))
        a.grid(alpha=.15)
    fig.suptitle('Archived dynamics | fixed time anchor | endpoint normalization, no fitting',fontsize=12)
    save(fig,'response_shape_limits')
    print('Saved two figure pairs from the archived matching curves.')


if __name__=='__main__':main()
