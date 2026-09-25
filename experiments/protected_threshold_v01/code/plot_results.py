#!/usr/bin/env python3
"""Scientific figures from archived conditional calculations."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
FIG=ROOT/'figures'; FIG.mkdir(exist_ok=True)
data=np.genfromtxt(ROOT/'results/trajectories.csv',names=True,delimiter=',',dtype=None,encoding='utf8')
inputs=json.loads((ROOT/'results/inputs.json').read_text())
budget=json.loads((ROOT/'results/budgets.json').read_text())
summary=json.loads((ROOT/'results/summary.json').read_text())
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                     'axes.grid':True,'grid.alpha':.15,'savefig.dpi':180})
colors=['#222222','#0072B2','#009E73','#E69F00','#C44745']
fig,ax=plt.subplots(1,2,figsize=(11.4,4.9))
for rec,color in zip(summary['summaries'],colors):
    r=data[data['case']==rec['case']]; eta=rec['eta_soft']
    lab='Clock F term only' if eta==0 else r'$\eta_{\rm soft}=$'+f'{eta:g}'
    ls='--' if eta==0 else '-'
    ax[0].plot(r['z'],r['u'],ls,color=color,lw=1.8,label=lab,zorder=6 if eta==0 else 3)
    ax[1].plot(r['z'],r['delta_alpha']*1e6,ls,color=color,lw=1.8,zorder=6 if eta==0 else 3)
for a in ax:
    a.set_xlim(4.2,0); a.set_xlabel('Redshift z (evolution proceeds to the right)')
ax[0].set_ylabel(r'Clock displacement $\chi-\chi_i$')
ax[1].set_ylabel(r'$[\alpha(z)/\alpha(0)-1]$ (ppm)')
ax[0].legend(frameon=False,fontsize=9,loc='upper left')
ax[0].set_title('The clock with a paired charged spectrum')
ax[1].set_title('The electromagnetic response survives')
fig.text(.5,.013,'100 GeV illustrative spectrum; preselected mass law; global-SUSY matter diagnostic, not an observationally validated cosmology.',
         ha='center',fontsize=8.5,color='#555555')
fig.tight_layout(rect=(0,.055,1,1))
for ext in ('png','pdf'): fig.savefig(FIG/f'protected_clock_response.{ext}')
plt.close(fig)

old=ROOT.parent/'charged_threshold_v01/results/radiative_budget.json'
unprotected=json.loads(old.read_text())['records'][-1]['max_loop_force_over_tree']
pure=next(r for r in budget['records'] if r['mass_eV']==1e11 and r['mu_over_mi']==1)
values=[unprotected,pure['max_F_force_over_tree'],budget['large_soft_example']['max_force_over_tree']]
fig,ax=plt.subplots(1,2,figsize=(11.4,5.3),gridspec_kw={'width_ratios':[1,1.2]})
names=['One complex\nscalar','Paired spectrum\n+ clock F term','Paired spectrum\n+ 1 GeV soft mass']
for i,(v,c) in enumerate(zip(values,['#888888','#0072B2','#C44745'])):
    ax[0].plot([i,i],[1,v],color=c,lw=5,alpha=.6)
    ax[0].scatter(i,v,color=c,s=58,zorder=3)
    ax[0].annotate(f'{v:.2e}',(i,v),xytext=(0,10),textcoords='offset points',ha='center',fontsize=9)
ax[0].set_yscale('log'); ax[0].set_ylim(1e-43,1e61)
ax[0].set_xticks(range(3),names); ax[0].set_xlim(-.45,2.45)
ax[0].axhline(1,color='#555555',lw=.8)
ax[0].set_ylabel(r'Maximum $|U_{1,\chi}|/|U_{{\rm tree},\chi}|$')
ax[0].set_title('100 GeV: conditional force budgets')
soft=np.logspace(-20,-12,301)
for mu,color in zip((.5,1.,2.),['#E69F00','#0072B2','#009E73']):
    r=next(r for r in budget['records'] if r['mass_eV']==1e11 and r['mu_over_mi']==mu)
    eta=soft**2*1e22/(8*np.pi**2*inputs['potential_unit_eV4'])
    ax[1].loglog(soft,r['max_soft_force_per_eta']*eta,color=color,lw=1.8,label=r'$\mu/m_i=$'+str(mu))
ax[1].axhline(.1,color='#333333',lw=.9,ls='--',label='0.1 force-budget criterion')
ax[1].set(xlabel='Common soft parameter '+r'$m_{\rm soft}$ (eV)',
          ylabel=r'Maximum $|U_{{\rm soft},\chi}|/|U_{{\rm tree},\chi}|$',
          title='Sensitivity to additional symmetry breaking')
ax[1].legend(frameon=False,fontsize=8,loc='upper left')
fig.text(.5,.015,'Budgets evaluated on the same tree trajectory. Fixed finite matching inputs at each scale; these curves are not experimental bounds.',
         ha='center',fontsize=8.5,color='#555555')
fig.tight_layout(rect=(0,.065,1,1))
for ext in ('png','pdf'): fig.savefig(FIG/f'protection_and_soft_budget.{ext}')
plt.close(fig)
print('Saved two PNG/PDF figure pairs.')
