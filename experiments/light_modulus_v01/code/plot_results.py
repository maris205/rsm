#!/usr/bin/env python3
"""Plot archived deterministic outputs; no observations or cosmic-time inference."""
from pathlib import Path
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
def read(name):
    with (ROOT/'results'/name).open() as f:return list(csv.DictReader(f))
def col(rows,key):return np.array([float(r[key]) for r in rows])

def main():
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'pdf.fonttype':42,'savefig.facecolor':'white'})
    fig,axs=plt.subplots(2,2,figsize=(11,7.8),layout='constrained')
    blue,orange,green='#1764a0','#d07518','#25816d'
    profile=read('mode_profile.csv');valley=read('physical_valley.csv')
    stress=read('source_stress.csv');transient=read('local_transient.csv')
    bench=json.loads((ROOT/'results/summary.json').read_text())['benchmark']
    ax=axs[0,0];ax.semilogy(col(profile,'link'),np.abs(col(profile,'u')),color=blue,lw=2)
    ax.set(xlabel='Higgs-pair index a',ylabel=r'$|u_a|$',title='(a) Retained light mode: endpoint asymmetry')
    ax.text(.04,.88,r'$u_0=2.40\times10^{-61}$'+'\n'+r'$u_n=-0.943$',transform=ax.transAxes)
    ax.set_ylim(1e-63,10);ax.grid(alpha=.15)
    ax=axs[0,1];ax.plot(col(valley,'r'),col(valley,'V_over_S'),color=blue,lw=2,
                         label=r'Benchmark $\kappa=5.00\times10^{-24}$')
    for kap,color in [('0.1',orange),('1',green)]:
        rows=[r for r in stress if r['kappa']==kap and float(r['r'])<=1]
        ax.plot(col(rows,'r'),col(rows,'V_over_S'),color=color,ls='--',
                label=rf'Stress test $\kappa={kap}$')
    ax.set(xlabel=r'$r=-z/(2v^2)$',ylabel=r'$V_{\rm valley}/S$',
           title='(b) Relaxed heavy fields: no stationary point')
    ax.legend(fontsize=8,loc='upper right');ax.grid(alpha=.15)
    ax=axs[1,0];ax.plot(col(valley,'r'),col(valley,'radial_kinetic_B'),color=blue,lw=2,label=r'Radial $B(r)$')
    ax.plot(col(valley,'r'),col(valley,'kahler_metric_C'),color=green,lw=2,label=r'$K_{C\bar C}=1/B(r)$')
    ax.set(xlabel=r'$r=-z/(2v^2)$',ylabel='Dimensionless metric',
           title='(c) Zero-source quotient geometry')
    ax.legend(fontsize=9);ax.grid(alpha=.15)
    ax=axs[1,1];ax.plot(col(transient,'local_elapsed_seconds')*1e17,
                        col(transient,'hidden_metric_change')*100,color=orange,lw=2)
    ax.axhline(10,ls=':',color='0.4');ax.set(xlabel=r'Local elapsed time ($10^{-17}$ s)',
        ylabel='Hidden metric change (%)',title='(d) Fixed-source Minkowski transient')
    ax.text(.05,.72,'Initially at rest; no expansion\n10% diagnostic endpoint\n'+
        rf'$\Delta t={bench["transient_event_seconds"]*1e17:.3f}\times10^{{-17}}$ s',transform=ax.transAxes,fontsize=9)
    ax.grid(alpha=.15)
    fig.suptitle('Original Higgs chain with its uneaten chiral mode retained',fontsize=14)
    (ROOT/'figures').mkdir(exist_ok=True)
    for ext in ('png','pdf'):
        metadata={'CreationDate':None,'ModDate':None} if ext=='pdf' else {'Software':'Matplotlib; archived deterministic outputs'}
        fig.savefig(ROOT/'figures'/('light_modulus_diagnostics.'+ext),dpi=180,metadata=metadata)
    plt.close(fig)
    print('Saved light_modulus_diagnostics.png and .pdf')

if __name__=='__main__':main()
