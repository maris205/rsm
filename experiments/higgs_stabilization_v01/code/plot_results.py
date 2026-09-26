#!/usr/bin/env python3
"""Export a standalone diagnostic figure from the archived tables."""
from pathlib import Path
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
def rows(name):
    with (ROOT/'results'/name).open() as f:return list(csv.DictReader(f))

def main():
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'figure.dpi':140,'savefig.dpi':190})
    fig,axs=plt.subplots(2,2,figsize=(11.5,8.1))
    colors=['#31688e','#35b779','#b65a30']
    phase=rows('phase_lifting.csv');gauge=rows('gauge_scan.csv');moment=rows('momentum_profile.csv')
    b=json.loads((ROOT/'results/summary.json').read_text())['benchmark']
    ax=axs[0,0]
    for c,q in zip(colors,(2,3,5)):
        r=[x for x in phase if int(x['q'])==q and x['v_over_cutoff']=='0.9']
        ax.plot([int(x['n']) for x in r],[float(x['log10_degree']) for x in r],'.-',c=c,label=f'q = {q}')
    ax.scatter([127],[np.log10(float(b['degree_scientific']))],c='black',s=24,zorder=5)
    ax.annotate(r'$D=5.90\times10^{60}$',xy=(127,60.77),xytext=(47,71),
                arrowprops={'arrowstyle':'->','lw':.8})
    ax.set(xlabel='Original gauge factors n',ylabel=r'$\log_{10}D$',title='A  Original-field phase lifting')
    ax.legend(frameon=False,loc='upper left')
    ax.text(.04,.05,r'$\mu/v\ \propto\ (v/M_*)^{D-2}$',transform=ax.transAxes,
            bbox={'facecolor':'white','edgecolor':'#ddd','pad':5})
    ax=axs[0,1]
    eta=np.array([float(x['eta']) for x in gauge])
    vec=np.array([float(x['lightest_vector_eV'])/1e12 for x in gauge])
    ax.plot(eta,vec,c=colors[0],label='Lightest vector')
    ax.axhline(float(b['radial_mass_lambda_over_g_3_eV'])/1e12,c=colors[1],ls='--',label=r'Radial/singlet: $\lambda/g=3$')
    ax.axhline(1,c='#555',ls=':',label='Inherited 1 TeV gap criterion')
    ax.scatter([2],[next(float(x['lightest_vector_eV'])/1e12 for x in gauge if float(x['eta'])==2)],c='black',s=24)
    ax.set(xlabel=r'Extra gauge coupling ratio $\eta$',ylabel='Mass (TeV)',ylim=(0,2.35),
           title=r'B  Free spectrum at $F=D=0$')
    ax.legend(frameon=False,loc='upper left',fontsize=9)
    ax=axs[1,0]
    old=float(b['original_tau']);new=float(b['abstract_tau_after_gauging'])
    labels=['Original chain\n(flat chiral retained)','Extra gauge\nlocal Higgs currents','Extra gauge\ndistributed source']
    ax.bar(range(3),[old/1e-61,0,new/1e-61],color=[colors[0],'#a0a0a0',colors[2]],width=.58)
    ax.scatter([1],[0],marker='x',s=70,c='#555',zorder=4)
    ax.text(1,.35,'exact zero',ha='center')
    ax.text(2,new/1e-61+.35,r'includes $q^{-n}$ overlap',ha='center',fontsize=9)
    ax.set(xticks=range(3),xticklabels=labels,ylabel=r'Static cross coefficient $\tau\,/10^{-61}$',
           ylim=(-.4,9.5),title='C  Sources must be matched as well as masses')
    ax.tick_params(axis='x',labelsize=8.5)
    ax=axs[1,1]
    for c,eta in zip(colors,(.5,1.,2.)):
        r=[x for x in moment if float(x['eta'])==eta]
        ax.loglog([float(x['p_over_M']) for x in r],[float(x['ratio_to_abstract_static']) for x in r],c=c,label=fr'$\eta={eta:g}$')
    ax.set(xlabel=r'Euclidean momentum $p_E/M$',ylabel=r'$\tau_{\rm local}(p_E)/\tau_{\rm site}(0)$',
           title='D  Derivative exchange reappears at finite momentum')
    ax.legend(frameon=False,loc='lower right')
    ax.grid(which='major',alpha=.15)
    fig.suptitle('Higgs stabilization: a heavy spectrum does not preserve the local endpoint contact',fontsize=13,y=.99)
    fig.text(.5,.008,'Conditional tree-level models; no observational fit. Panel D does not bound all auxiliary-field effects.',ha='center',fontsize=9,color='#444')
    fig.tight_layout(rect=(0,.035,1,.965),h_pad=2,w_pad=2)
    (ROOT/'figures').mkdir(exist_ok=True)
    for ext in ('png','pdf'):
        metadata={'CreationDate':None,'ModDate':None} if ext=='pdf' else None
        fig.savefig(ROOT/'figures'/('stabilization_tradeoff.'+ext),metadata=metadata)
    plt.close(fig)
    print('Exported stabilization_tradeoff.png and .pdf')

if __name__=='__main__':main()
