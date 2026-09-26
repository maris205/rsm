#!/usr/bin/env python3
"""Plot frozen matching and selected static-response tables."""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'figures'
OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,
                     'axes.titleweight':'bold','savefig.dpi':180})

def rows(name):
    with (ROOT/'results'/name).open() as f: return list(csv.DictReader(f))

fig,axs=plt.subplots(1,2,figsize=(14,5.3),layout='constrained')
data=rows('chain_matching.csv')
colors=['#227c9d','#e28a39','#5c956d']
for q,color in zip((2,3,5),colors):
    rr=[r for r in data if float(r['q'])==q]
    axs[0].semilogy([int(r['n']) for r in rr],[float(r['tau']) for r in rr],
                     lw=2,color=color,label=f'q = {q}')
axs[0].axhline(2.9251519294813043e-61,color='black',ls=':',lw=1.1,
               label='Previous target at 1 TeV')
axs[0].scatter([128],[2.2617683795227874e-61],s=90,marker='*',color='#e28a39',zorder=5)
axs[0].annotate(r'$n=128$: $\tau=2.26\times10^{-61}$',
                 (128,2.2617683795227874e-61),xytext=(28,1e-42),
                 arrowprops={'arrowstyle':'->','color':'#6b6b6b'},fontsize=10)
axs[0].set(xlim=(2,260),ylim=(1e-120,2),xlabel='Number of mediator sites, n',
            ylabel=r'Endpoint cross coefficient $\tau$',
            title='(a) A chain generates a tiny cross coupling')
axs[0].legend(loc='lower left',fontsize=9)
data=rows('response_scan.csv')
for lam,color in [(1e12,'#8493a5'),(2e12,'#227c9d')]:
    rr=[r for r in data if float(r['q'])==3 and float(r['Lambda_eV'])==lam]
    x=np.array([float(r['mG_eV']) for r in rr])
    y=np.array([float(r['max_static_force_ratio']) for r in rr])
    ym=np.array([float(r['max_once_rematched_force_ratio'] or 'nan') for r in rr])
    axs[1].loglog(x,y,color=color,lw=2,label=fr'$\Lambda={lam/1e12:g}$ TeV: frozen matching')
    axs[1].loglog(x,ym,color=color,lw=1.4,ls='--',label=fr'$\Lambda={lam/1e12:g}$ TeV: once rematched')
axs[1].axhline(.1,color='black',ls=':',lw=1.1,label='Declared 10% response threshold')
axs[1].axvline(1e-6,color='#c1534c',ls=':',lw=1)
axs[1].text(1.8e-6,1e-7,'Previous\n$m_G$',color='#c1534c',fontsize=10)
axs[1].scatter([1e-15],[8.954500419269619e-6],color='#227c9d',s=45,zorder=5)
axs[1].set(xlim=(1e-20,1e-4),ylim=(1e-20,1e27),xlabel=r'$m_G$ [eV]',
            ylabel=r'Maximum selected static force / $2U$',
            title='(b) The hidden self response is not suppressed')
axs[1].legend(loc='upper left',fontsize=8.5)
for ax in axs: ax.grid(alpha=.15)
fig.supxlabel('Local endpoint couplings are order one; the complete matching retains '
              r'$-\,|X_c|^4/\Lambda^2$.'
              '\nRight panel: aligned phase, old clock interval, selected tree response only. '
              'The 1 TeV cases also fail the stated charged-mass gap condition.',fontsize=10)
for ext in ('png','pdf'):
    fig.savefig(OUT/f'chain_matching_and_response.{ext}',bbox_inches='tight')
plt.close(fig)
print('Rendered matching and response figure (PNG/PDF).')
