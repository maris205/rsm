#!/usr/bin/env python3
"""Figures from frozen threshold/RG output tables only."""
from pathlib import Path
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,
                     'axes.titleweight':'bold','savefig.dpi':180})
def rows(name):
    with (ROOT/'results'/name).open() as f:return list(csv.DictReader(f))

fig,axs=plt.subplots(1,2,figsize=(13.5,5.0),layout='constrained')
r=rows('scale_curve.csv');mu=np.array([float(x['mu_over_Lambda']) for x in r])
for key,color,label,style in [('cL_hat','#217e99','Running local coefficient','-'),
                             ('explicit_loop_force_hat','#d7853c','Explicit light-loop term','--'),
                             ('total_force_hat','#333333','Their sum: unchanged force','-')]:
    axs[0].semilogx(mu,[float(x[key]) for x in r],lw=2,color=color,ls=style,label=label)
summary=json.loads((ROOT/'results/summary.json').read_text())
p=summary['primary'];mx_ratio=p['mX_eV']/p['Lambda_eV']
axs[0].axvline(mx_ratio,color='#888888',ls=':',lw=1)
axs[0].text(mx_ratio*1.6,-.3,r'$\mu=m_X$',color='#666666',fontsize=10)
axs[0].set(xlabel=r'Renormalization scale $\mu/\Lambda$',
           ylabel=r'Source force coefficient: $\Lambda^4\partial_S\Delta V/(2S)$',
           title='(a) Changing scale moves the logarithm',xlim=(1e-14,1))
axs[0].legend(loc='lower left',fontsize=9)
r=rows('mass_scan.csv');mg=np.array([float(x['mG_eV']) for x in r])
axs[1].loglog(mg,[float(x['max_tree_force_ratio']) for x in r],color='#d7853c',ls='--',lw=2,
              label='New finite-auxiliary tree term')
axs[1].loglog(mg,[float(x['max_selected_total_force_ratio']) for x in r],color='#217e99',lw=2,
              label='Tree + selected one-loop sector')
axs[1].axhline(.1,color='#333333',ls=':',lw=1,label='Declared 10% condition')
axs[1].scatter([1e-15],[p['max_selected_total_force_ratio']],color='#217e99',s=50,zorder=5)
axs[1].annotate(r'$m_G=10^{-15}$ eV: $R=2.00\times10^{-5}$',
               (1e-15,p['max_selected_total_force_ratio']),xytext=(2e-19,2e-2),fontsize=9,
               arrowprops={'arrowstyle':'->','color':'#666666'})
bound=summary['selected_total_10percent_mG_eV']
axs[1].plot(bound,.1,'o',color='#333333',ms=5)
axs[1].set(xlabel=r'$m_G$ [eV]',ylabel=r'Maximum selected force / $2U$',
           title='(b) The low-scale selected force remains small',
           xlim=(1e-20,1e-12),ylim=(1e-20,1e5))
axs[1].legend(loc='lower right',fontsize=9)
for ax in axs:ax.grid(alpha=.15)
fig.supxlabel(r'$q=3$, $n=127$, $\Lambda=2$ TeV; fixed finite boundary, leading $S^2$, selected Gaussian sector.'
              '\nScale consistency is not complete UV matching. Chain locality and the full supergravity theory remain open.',fontsize=10)
for ext in ('png','pdf'):fig.savefig(OUT/f'threshold_running_and_force.{ext}',bbox_inches='tight')
plt.close(fig)
print('Rendered threshold/RG figure (PNG/PDF).')
