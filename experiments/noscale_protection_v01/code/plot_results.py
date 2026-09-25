#!/usr/bin/env python3
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
def read(p):
 with p.open() as f:return list(csv.DictReader(f))
def main():
 rows=read(ROOT/'results/scheme_scan.csv');limits=read(ROOT/'results/budget_limits.csv');reps=read(ROOT/'results/representative_curves.csv')
 old=read(REPO/'experiments/hidden_transfer_v01/results/scan.csv')
 def select(a=1,phase='zero',scheme='tree_eta',mu=1):
  return [r for r in rows if float(r['a'])==a and r['phase']==phase and r['scheme']==scheme and float(r['mu_factor'])==mu]
 def arr(rs,key):return np.array([float(r[key]) for r in rs])
 def save(fig,name):
  for ext in ['png','pdf']:fig.savefig(ROOT/'figures'/(name+'.'+ext),dpi=200,bbox_inches='tight')
  plt.close(fig)
 plt.rcParams.update({'font.size':10.5,'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False})
 teal='#007f8b';orange='#b56a20';purple='#7253a2';navy='#284b73'
 fig,ax=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
 r=select();m=arr(r,'mG_eV')
 o=[r for r in old if r['phase']=='zero' and 1e-12<=float(r['mG_eV'])<=1e3]
 ax[0].loglog(arr(o,'mG_eV'),arr(o,'max_charged_CW_force_ratio'),color='#999999',lw=2,label='Prior additive hidden model: charged CW')
 ax[0].loglog(m,arr(r,'max_charged_force_ratio'),color=teal,lw=2,label='Quartic geometry: charged CW')
 ax[0].loglog(m,arr(r,'max_modulus_force_ratio'),color=orange,lw=2,label='Quartic geometry: heavy-scalar CW')
 ax[0].axhline(.1,color='black',ls=':',lw=1.2,label='Chosen 10% local force budget')
 ax[0].set(xlim=(1e-12,1e3),xlabel=r'$m_G$ (eV)',ylabel='Component force / reference clock force',title='(a) The large charged feedback can be suppressed')
 ax[0].legend(fontsize=8.5,loc='upper left');ax[0].grid(alpha=.15)
 rr=select(scheme='scalar_retuned_eta');mm=arr(rr,'mG_eV')
 for key,color,style,label in [
 ('max_modulus_force_ratio',orange,'--','Heavy-scalar CW'),
 ('max_tree_force_ratio',purple,':',r'Tree force from retuned $\eta$'),
 ('max_charged_force_ratio',teal,'-.',r'Charged CW with retuned tree spectrum'),
 ('max_component_sum_ratio',navy,'-','Sum of component norms (diagnostic)')]:
  ax[1].loglog(mm,arr(rr,key),color=color,ls=style,lw=1.7 if key!='max_component_sum_ratio' else 2.5,label=label)
 ax[1].axhline(.1,color='black',ls=':',lw=1.1)
 bound=next(float(x['component_sum_10percent_mG_eV']) for x in limits if float(x['a'])==1 and x['phase']=='zero' and float(x['mu_factor'])==1 and x['scheme']=='scalar_retuned_eta')
 ax[1].axvline(bound,color=navy,ls=':',lw=1)
 ax[1].text(.55,.06,r'$m_G\simeq0.316$ MeV'+'\n'+'selected-term budget only',transform=ax[1].transAxes,fontsize=9)
 ax[1].set(xlim=(1e-12,1e12),xlabel=r'$m_G$ (eV)',ylabel='Component force / reference clock force',title='(b) Vacuum adjustment also changes charged masses')
 ax[1].legend(fontsize=8.1,loc='upper left');ax[1].grid(alpha=.15)
 fig.suptitle(r'Specified local finite contributions | $a=b=1$, $\theta=0$, fixed matching scales',fontsize=12)
 save(fig,'protection_and_radiative_vacuum')
 fig,ax=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
 grid=np.array([np.log10(arr(select(a=a,scheme='scalar_retuned_eta'),'max_component_sum_ratio')/.1) for a in (.1,1,10)])
 im=ax[0].imshow(grid,aspect='auto',origin='lower',extent=(-12.125,12.125,-.5,2.5),cmap='RdBu_r',norm=TwoSlopeNorm(vmin=-30,vcenter=0,vmax=25),interpolation='nearest')
 ax[0].set(yticks=[0,1,2],yticklabels=['0.1','1','10'],xlabel=r'$\log_{10}(m_G/{\rm eV})$',ylabel=r'Quartic coefficient $a=b$',title='(a) Retuned diagnostic, three sampled geometries')
 for j,a in enumerate((.1,1,10)):
  lim=next(float(x['component_sum_10percent_mG_eV']) for x in limits if float(x['a'])==a and x['phase']=='zero' and float(x['mu_factor'])==1 and x['scheme']=='scalar_retuned_eta')
  ax[0].plot(np.log10(lim),j,marker='|',color='black',ms=20,mew=2)
 cb=fig.colorbar(im,ax=ax[0],pad=.02,extend='both');cb.set_label(r'$\log_{10}(R_{\rm sum}/0.1)$',fontsize=9)
 points=[]
 for phase in ['zero','quadrature','pi']:
  points.append(next(r for r in reps if r['phase']==phase and float(r['mG_eV'])==1 and int(r['index'])==1024))
 scale=abs(float(points[0]['tree_eta_modulus_CW_x']))
 x=np.arange(3)
 ax[1].bar(x-.17,[float(p['tree_eta_modulus_CW_x'])/scale for p in points],.32,color=orange,label=r'$\chi$ force component')
 ax[1].bar(x+.17,[float(p['tree_eta_modulus_CW_y'])/scale for p in points],.32,color=purple,label=r'$y$ force component')
 ax[1].axhline(0,color='black',lw=.8)
 ax[1].set(xticks=x,xticklabels=[r'$0$',r'$\pi/2$',r'$\pi$'],xlabel=r'Phase $\theta$',ylabel='Heavy-scalar force / aligned magnitude',title='(b) Quadrature moves the force into the transverse field')
 ax[1].legend(fontsize=9,loc='upper center')
 ax[1].text(.39,.11,'Full geometric derivative:'+'\n'+r'transverse coefficient = $1/3$.',transform=ax[1].transAxes,fontsize=9)
 fig.suptitle('Conditional EFT diagnostics; no observational exclusion or complete-loop stability claim',fontsize=12)
 save(fig,'geometry_budget_and_phase')
 print('Rendered two figure pairs from archived CSV files.')
if __name__=='__main__':main()
