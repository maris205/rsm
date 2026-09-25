#!/usr/bin/env python3
"""Plot archived outputs, without rerunning models or fitting observations."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
FIG=ROOT/'figures'
FIG.mkdir(exist_ok=True)
data=np.genfromtxt(ROOT/'results/trajectories.csv',delimiter=',',names=True,
                   dtype=None,encoding='utf8')
summary=json.loads((ROOT/'results/summary.json').read_text())
budget=json.loads((ROOT/'results/radiative_budget.json').read_text())
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                     'axes.grid':True,'grid.alpha':0.17,'figure.dpi':130,'savefig.dpi':180})
colors={'probe':'#333333','loop_1':'#0072B2','loop_100':'#E69F00','loop_10000':'#C44745'}
labels={'probe':'Tree / probe control','loop_1':r'$\ell=1$',
        'loop_100':r'$\ell=100$','loop_10000':r'$\ell=10^4$'}

fig,axes=plt.subplots(2,2,figsize=(11.4,7.8),sharex=True)
for j,kin in enumerate(('chi','R')):
    for suffix,color in colors.items():
        r=data[data['case']==f'{kin}_{suffix}']
        style='--' if suffix=='probe' else '-'
        for ax,key,scale in ((axes[0,j],'u',1),(axes[1,j],'delta_alpha',1e6)):
            ax.plot(r['z'],r[key]*scale,style,color=color,lw=2 if suffix=='probe' else 1.6,
                    label=labels[suffix],zorder=4 if suffix=='probe' else 3)
    axes[0,j].set_title(r'$K=1$ (canonical $\chi$)' if kin=='chi'
                        else r'$K=e^{2(\chi-\chi_i)}$ (constant $R$ kinetic term)')
    axes[1,j].set_xlabel('Redshift z (evolution proceeds to the right)')
    axes[1,j].set_xlim(4.2,0)
    axes[1,j].axhline(0,color='#555555',lw=.7)
    axes[0,j].axhline(0,color='#555555',lw=.7)
axes[0,0].set_ylabel(r'Clock displacement $\chi-\chi_i$')
axes[1,0].set_ylabel(r'$[\alpha(z)/\alpha(0)-1]$ (ppm)')
axes[0,0].legend(frameon=False,loc='upper left',fontsize=9)
fig.suptitle('A mass threshold transmits the response; its vacuum loop changes the clock',fontsize=14,y=.98)
fig.text(.5,.016,'Conditional vacuum EFT: 0.79–7.89 meV charged-mass examples; no observational fit or optical-spectrum prediction.',
         ha='center',fontsize=9,color='#555555')
fig.tight_layout(rect=(0,.043,1,.95))
for ext in ('png','pdf'): fig.savefig(FIG/f'clock_and_threshold.{ext}')
plt.close(fig)

fig,axes=plt.subplots(1,2,figsize=(11.4,4.7))
for kin,style in (('chi','-'),('R','--')):
    for suffix in ('loop_1','loop_100','loop_10000'):
        r=data[data['case']==f'{kin}_{suffix}']
        axes[0].semilogy(r['z'],abs(r['loopp'])/(2*r['tree']),style,color=colors[suffix],
                         lw=1.6,label=labels[suffix]+(', K=1' if kin=='chi' else ', R kinetic'))
axes[0].axhline(1,color='#333333',lw=.9)
axes[0].set(xlim=(4.2,0),xlabel='Redshift z',ylabel=r'$|U_{1,\chi}|/|U_{{\rm tree},\chi}|$',
            title='Feedback force on the computed trajectories')
axes[0].legend(fontsize=8,frameon=True,framealpha=.95,ncol=2,loc='lower left',edgecolor='none')
mass=np.logspace(-4,11,301)
alpha=1/137.035999084
ratio=3*mass**4*1e-6/(4*np.pi*alpha*budget['rho_crit_ref_eV4'])
axes[1].loglog(mass,ratio,color='#704C90',lw=2)
axes[1].axhline(1,color='#333333',lw=.9)
axes[1].axvspan(.0007887683448493857,.007887683448493858,color='#0072B2',alpha=.15,
                label='Masses used in ODE examples')
for rec,label in zip(budget['records'],('1 eV','1 MeV','100 GeV')):
    axes[1].scatter(rec['mass_initial_eV'],rec['ppm_budget_over_rho_crit_ref'],s=28,color='#704C90')
    axes[1].annotate(label,(rec['mass_initial_eV'],rec['ppm_budget_over_rho_crit_ref']),
                      xytext=(-4,8),textcoords='offset points',ha='right',fontsize=9)
axes[1].set(xlabel='Charged threshold mass (eV)',ylabel=r'$|\Delta U_1|/\rho_{\rm crit,ref}$',
            title='Unprotected vacuum-energy cost of a 1 ppm response')
axes[1].legend(fontsize=8,frameon=False,loc='lower right')
fig.text(.5,.013,r'Fixed $\overline{\rm MS}$ scale $\mu=m_i$; constant vacuum term adjusted once. No slope cancellation or protecting mechanism.',
         ha='center',fontsize=9,color='#555555')
fig.tight_layout(rect=(0,.047,1,1))
for ext in ('png','pdf'): fig.savefig(FIG/f'feedback_and_budget.{ext}')
plt.close(fig)
print('Saved two figures in PNG and PDF.')
