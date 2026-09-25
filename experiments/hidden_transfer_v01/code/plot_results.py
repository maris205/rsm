#!/usr/bin/env python3
"""Archived hidden-transfer diagnostics; no fits or new background evolution."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]


def read(name):
    with (ROOT/'results'/name).open() as f:return list(csv.DictReader(f))


def main():
    out=ROOT/'figures';out.mkdir(exist_ok=True)
    s=json.loads((ROOT/'results/summary.json').read_text());lim=s['limits']
    scan=read('scan.csv');valley=read('valleys.csv')
    def arr(rows,key):return np.array([float(r[key]) for r in rows])
    p={k:[r for r in scan if r['phase']==k and float(r['mG_eV'])>0] for k in ('zero','quadrature','pi')}
    m=arr(p['zero'],'mG_eV');eft=lim['mG_for_VA_MI_eV'];bound=lim['phases'][0]['charged_CW_limit_eV']
    plt.rcParams.update({'font.size':10.5,'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False})
    def save(fig,name):
        for ext in ('png','pdf'):fig.savefig(out/(name+'.'+ext),dpi=200,bbox_inches='tight')
        plt.close(fig)
    fig,ax=plt.subplots(1,2,figsize=(11.5,4.6),layout='constrained')
    for phase,color,ls,label in [('zero','#c16526','-',r'Axial extra force: $\theta=0,\pi$'),
                                 ('quadrature','#6650a4','--',r'Axial extra force: $\theta=\pi/2$')]:
        ax[0].loglog(m,arr(p[phase],'max_initial_axis_tree_force_ratio'),color=color,ls=ls,lw=1.6,label=label)
    ax[0].loglog(m,arr(p['quadrature'],'max_charged_CW_force_ratio'),color='#007f8b',lw=2.3,label=r'Charged CW force: $\theta=\pi/2$')
    ax[0].axhline(.1,color='black',ls=':',lw=1.1,label='Chosen 10% force budget')
    ax[0].set(xlabel=r'$m_G$ (eV)',ylabel='Force / old clock force',
              title='(a) Frozen-axis budgets; not new trajectories',xlim=(1e-36,1e3))
    ax[0].legend(fontsize=8,loc='upper left')
    ax[0].grid(alpha=.15)
    ax[1].loglog(m,arr(p['zero'],'VA_over_MI'),color='#335780',lw=2)
    ax[1].axhline(1,color='black',ls=':',lw=1)
    ax[1].axhline(10,color='#555555',ls='--',lw=1)
    ax[1].axvspan(1e-36,bound,color='#007f8b',alpha=.10,label='Formal CW budget passes')
    ax[1].axvspan(eft,1e3,color='#335780',alpha=.12,label='100 GeV below VA marker')
    ax[1].axvline(bound,color='#007f8b',ls='--',lw=1)
    ax[1].axvline(eft,color='#335780',ls='--',lw=1)
    ax[1].set(xlabel=r'$m_G$ (eV)',ylabel=r'$\Lambda_{\rm VA}/(100\ {\rm GeV})$',
              title='(b) Nominal EFT marker and small-force region',xlim=(1e-36,1e3))
    ax[1].legend(fontsize=8.5,loc='upper left')
    ax[1].text(.38,.13,r'CW: $m_G\lesssim2.0\times10^{-17}$ eV'+'\n'+r'VA marker: $m_G\gtrsim2.4\times10^{-6}$ eV',
               transform=ax[1].transAxes,fontsize=9)
    ax[1].grid(alpha=.15)
    fig.suptitle('Specified additive hidden sector | low-marker heavy loops are formal extrapolations',fontsize=12)
    save(fig,'hidden_force_and_domain')

    fig,ax=plt.subplots(1,2,figsize=(11.5,4.4),layout='constrained')
    for n,color in [(0,'#c16526'),(512,'#6650a4'),(1024,'#007f8b')]:
        r=[x for x in valley if int(x['index'])==n]
        ax[0].semilogx(arr(r,'beta'),arr(r,'VminusLambda_over_U'),marker='o',color=color,lw=1.5,
                       label=f"z = {abs(float(r[0]['z'])):.2f}")
    ax[0].axhline(.74985,color='#555555',ls='--',lw=1.2,label=r'$b-1/4=0.74985$')
    ax[0].set(xlabel=r'$m_G/H_{\rm ref}$',ylabel=r'$(V(y_*)-\rho_\Lambda)/U$',
              title=r'(a) $\theta=\pi/2$: a classical heavy-field valley')
    ax[0].ticklabel_format(axis='y',useOffset=False)
    ax[0].legend(fontsize=9,loc='lower right');ax[0].grid(alpha=.15)
    ax[1].loglog(m,arr(p['zero'],'max_charged_CW_force_ratio'),color='#007f8b',lw=2,label='Gravity-induced charged CW budget')
    ax[1].loglog(m,arr(p['zero'],'gauge_common_old_force_ratio_cX1'),color='#b56a20',lw=2,ls='--',label=r'Extra gauge channel, $c_X=1$: common $d$ only')
    ax[1].axhline(.1,color='black',ls=':',lw=1)
    ax[1].set(xlabel=r'$m_G$ (eV)',ylabel='Component force / old clock force',
              title='(b) Direct X-gauge coupling adds a transmission path',xlim=(1e-36,1e3))
    ax[1].legend(fontsize=8,loc='upper left');ax[1].grid(alpha=.15)
    fig.suptitle('Tree-level exception retained | gauge curve is a partial leading-log diagnostic',fontsize=12)
    save(fig,'heavy_valley_and_gauge_channel')
    print('Saved two figure pairs from archived results.')


if __name__=='__main__':main()
