#!/usr/bin/env python3
"""Fixed conditional surrogate comparisons; ranks are not calibrated p-values."""
from pathlib import Path
import csv
import hashlib
import json
import platform

import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results'
FIG = ROOT / 'figures'
INPUT_HASHES = {}
CHECKS = []


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name):
    path = ROOT / name
    INPUT_HASHES[name] = sha(path)
    with np.load(path, allow_pickle=False) as data:
        return {k: data[k] for k in data.files}


def read_json(name):
    path=ROOT/name
    INPUT_HASHES[name]=sha(path)
    return json.loads(path.read_text())


def check(name, passed, **evidence):
    CHECKS.append(dict(name=name,passed=bool(passed),**evidence))


def empirical_rank(real, null):
    n=len(null)
    return float(min(1., 2*min((1+np.count_nonzero(null<=real))/(n+1),
                              (1+np.count_nonzero(null>=real))/(n+1))))


def describe(real, null, baseline):
    q=np.quantile(null,[.025,.25,.5,.75,.975])
    return dict(real_T=float(real),real_D=float(real-baseline),surrogates=len(null),
                null_q025=float(q[0]),null_q25=float(q[1]),null_median=float(q[2]),
                null_q75=float(q[3]),null_q975=float(q[4]),null_IQR=float(q[3]-q[1]),
                empirical_two_sided_rank=empirical_rank(real,null),
                inside_central_95=bool(q[0]<=real<=q[4]),
                real_minus_null_median=float(real-q[2]))


def main():
    FIG.mkdir(exist_ok=True)
    inp=load('data/inputs.npz')
    main=load('results/lattice_main.npz')
    fine=load('results/lattice_refined_selected.npz')
    r64=load('results/reference_N64_selected.npz')
    r128=load('results/reference_N128_selected.npz')
    ins=read_json('results/input_summary.json')
    refsummary=read_json('results/reference_summary.json')
    freeze=read_json('results/protocol_freeze.json')
    quality_path=OUT/'input_quality.csv'
    INPUT_HASHES['results/input_quality.csv']=sha(quality_path)
    with quality_path.open() as h:
        quality=list(csv.DictReader(h))
    check('frozen_protocol_unchanged',sha(ROOT/'protocol.md')==freeze['sha256'])
    check('case_alignment',np.array_equal(main['case_indices'],np.arange(1593))
          and np.array_equal(inp['case_family'],main['case_family'])
          and np.array_equal(inp['case_block'],main['case_block']))
    check('common_initial_energy',np.ptp(main['Htotal'][:,0])<1e-13)
    fam,blk,num=inp['case_family'],inp['case_block'],inp['case_surrogate']
    T=(main['Hclock'][:,-1]-main['Hclock'][:,0])/main['Htotal'][:,0]
    check('primary_reconstructed',np.max(np.abs(T-main['primary_T']))<1e-15)
    baseline=float(T[0]); rows=[]; agg={}; cubes={}
    for family in ('iaaft','exact_spectrum'):
        real=[]; null=[]
        for b in range(8):
            ids=np.flatnonzero((fam==family)&(blk==b))
            rid=np.flatnonzero((fam=='real')&(blk==b))
            check(f'{family}_block{b}_complete',len(ids)==99 and len(rid)==1
                  and np.array_equal(num[ids],np.arange(1,100)))
            x=float(T[rid[0]]); y=T[ids]
            real.append(x); null.append(y)
            quality_pass=int(np.sum(inp['primary_node_quality_pass'][ids])) if family=='iaaft' else None
            rows.append(dict(family=family,block=b+1,gap_start=4097+128*b,
                             gap_end=4224+128*b,node_quality_pass_count=quality_pass,
                             **describe(x,y,baseline)))
        cube=np.asarray(null)
        cubes[family]=cube
        agg[family]=describe(np.mean(real),np.mean(cube,axis=0),baseline)
    # Additional implementation diagnostic only: do not replace frozen primary.
    valid=np.ones(99,dtype=bool)
    for b in range(8):
        valid &= inp['primary_node_quality_pass'][(fam=='iaaft')&(blk==b)]
    quality_sensitivity={f:describe(agg[f]['real_T'],np.mean(cubes[f][:,valid],axis=0),baseline)
                         for f in cubes}
    primary_rows=[r for r in rows if r['family']=='iaaft']
    iqr=np.asarray([r['null_IQR'] for r in primary_rows])
    ids=fine['case_indices']
    fineT=(fine['Hclock'][:,-1]-fine['Hclock'][:,0])/fine['Htotal'][:,0]
    fineD=fineT-fineT[np.flatnonzero(ids==0)[0]]
    contrast_error=fineD-(T[ids]-baseline)
    scale=np.array([iqr[int(blk[i])] if i else np.nan for i in ids])
    refinement_ratio=float(np.max(np.abs(contrast_error[ids!=0])/scale[ids!=0]))
    check('time_error_below_five_percent_IQR',refinement_ratio<.05,value=refinement_ratio)
    ids64=r64['caseids']; ids128=r128['caseids']
    check('reference_same_cases_times',np.array_equal(ids64,ids128)
          and np.array_equal(main['t'],r64['t']) and np.array_equal(main['t'],r128['t']))
    tref=[]
    for data in (r64,r128):
        tc=(data['clock_energy'][:,-1]-data['clock_energy'][:,0])/data['total_energy'][:,0]
        tref.append(tc)
    d64=tref[0]-tref[0][np.flatnonzero(ids64==0)[0]]
    d128=tref[1]-tref[1][np.flatnonzero(ids128==0)[0]]
    spatial_rows=[]
    for j,cid in enumerate(ids64):
        if cid==0:
            continue
        b=int(blk[cid]); err=float(d128[j]-d64[j])
        spatial_rows.append(dict(case_index=int(cid),family=str(fam[cid]),block=b+1,
                                 contrast_N64=float(d64[j]),contrast_N128=float(d128[j]),
                                 spatial_contrast_change=err,change_over_IAAFT_IQR=abs(err)/iqr[b]))
    spatial_ratio=max(r['change_over_IAAFT_IQR'] for r in spatial_rows)
    # This is a diagnostic scale comparison, not a global convergence proof.
    check('selected_spatial_change_below_five_percent_IQR',spatial_ratio<.05,value=spatial_ratio)
    ode_error=float(np.max(np.abs(T[ids64]-tref[0])))
    check('main_vs_independent_primary',ode_error<1e-6,value=ode_error)
    quality_stats={}
    for family in ('iaaft','exact_spectrum'):
        qq=[q for q in quality if q['family']==family]
        def values(k): return np.array([float(q[k]) for q in qq])
        quality_stats[family]={k:dict(median=float(np.median(values(k))),min=float(np.min(values(k))),
                                    max=float(np.max(values(k)))) for k in
             ('fft_amplitude_relative_l2_no_dc','fft_power_relative_l2_no_dc',
              'sorted_relative_l2_centered_denominator','dense_windowed_fft_power_relative_l2_no_dc',
              'dense_derivative_fft_power_relative_l2_no_dc','dense_windowed_rms_relative_difference',
              'dense_derivative_rms_relative_difference')}
    realids=np.flatnonzero(fam=='real')
    c1=main['C1']
    c1err=np.sqrt(np.mean((c1[realids]-c1[0])**2,axis=1)/np.mean(c1[0]**2))
    secondary=dict(real_C1_relative_rms_vs_baseline=c1err.tolist(),
                   real_clock_gain_fractional_change_vs_baseline=((T[realids]-baseline)/baseline).tolist(),
                   meaning='Descriptive only; no endpoint replacement or second significance claim.')
    result=dict(status='COMPLETED_WITH_CONTROL_LIMITATIONS',scope='Conditional response discrimination, no observed physical targets or calibrated significance.',
                primary='T=(Hclock(20)-Hclock(0))/Htotal(0); D=T-T_baseline',baseline_T=baseline,
                rank_definition='min(1,2*min((1+count(null<=real))/(M+1),(1+count(null>=real))/(M+1))); conditional empirical rank, not calibrated p-value.',
                blocks=rows,aggregate=agg,secondary=secondary,
                quality_sensitivity=dict(note='Additional input-quality diagnostic; excludes across-block surrogate columns containing any node-quality failure, does not replace frozen primary.',
                                         retained_columns=int(valid.sum()),excluded_surrogate_numbers=(np.flatnonzero(~valid)+1).tolist(),results=quality_sensitivity),
                input_quality_summary=quality_stats,input_quality_failures=ins['surrogates']['failed_cases_retained'],
                numerical_precision=dict(max_time_contrast_change=float(np.max(np.abs(contrast_error))),
                                         max_time_contrast_change_over_IQR=refinement_ratio,
                                         max_selected_spatial_contrast_change=max(abs(r['spatial_contrast_change']) for r in spatial_rows),
                                         max_selected_spatial_change_over_IQR=spatial_ratio,
                                         max_main_vs_independent_T_error=ode_error,
                                         note='Selected resolution diagnostics; two spatial resolutions do not establish an order. Original reference accuracy failures are retained in reference_initial_*.'),
                spatial_selected=spatial_rows,checks=CHECKS,numerical_comparison_checks_passed=all(c['passed'] for c in CHECKS),
                all_input_matching_checks_passed=False,
                limitations=['Four IAAFT node-quality failures retained in primary ensemble.',
                             'Windowed continuous driver and derivative spectra are not strictly matched.',
                             'Conditional ranks are not calibrated p-values; no independent physical target.',
                             'Absence of an extreme response is not equivalence or an exclusion of other bridges.',
                             'Exact-spectrum controls are paired with IAAFT, not an independent replication.'],
                input_sha256=INPUT_HASHES,code_sha256=sha(Path(__file__)),
                software=dict(python=platform.python_version(),numpy=np.__version__,matplotlib=matplotlib.__version__))
    (OUT/'comparison_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    for filename,data in (('comparison_blocks.csv',rows),('comparison_spatial.csv',spatial_rows)):
        with (OUT/filename).open('w',newline='') as h:
            wr=csv.DictWriter(h,fieldnames=list(data[0]));wr.writeheader();wr.writerows(data)
    np.savez_compressed(OUT/'comparison_distributions.npz',
                        real=T[realids],iaaft=cubes['iaaft'],exact_spectrum=cubes['exact_spectrum'],
                        baseline=baseline,quality_common_columns=valid)
    plots(inp,main,rows,agg,cubes,quality,baseline)
    print(json.dumps({k:result[k] for k in ('aggregate','quality_sensitivity','numerical_precision','secondary',
                                         'numerical_comparison_checks_passed')},indent=2))
    print(f'{len(CHECKS)} numerical comparisons; input matching limitations retained.')
    if not result['numerical_comparison_checks_passed']:
        raise SystemExit('A numerical comparison gate failed; inspect saved diagnostics.')


def plots(inp,main,rows,agg,cubes,quality,baseline):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'axes.titleweight':'semibold','legend.frameon':False,'savefig.dpi':180})
    colors={'iaaft':'#2478a1','exact_spectrum':'#68985b','real':'#c45445'}
    labels_by_family={'iaaft':'IAAFT','exact_spectrum':'Exact node spectrum'}
    fig,axs=plt.subplots(1,2,figsize=(11.3,4.7),layout='constrained')
    fig.suptitle('Arithmetic input: a fixed energy-exchange comparison',fontsize=14,weight='semibold')
    for f,offset in (('iaaft',-.13),('exact_spectrum',.13)):
        sub=[r for r in rows if r['family']==f]
        lo=np.array([r['null_q025'] for r in sub]);mid=np.array([r['null_median'] for r in sub]);hi=np.array([r['null_q975'] for r in sub])
        axs[0].errorbar(np.arange(1,9)+offset,1e4*(mid-baseline),yerr=1e4*np.array([mid-lo,hi-mid]),
                        fmt='o',ms=4,capsize=3,color=colors[f],label=labels_by_family[f]+' central 95%')
    real=np.array([r['real_T'] for r in rows if r['family']=='iaaft'])
    axs[0].scatter(np.arange(1,9),1e4*(real-baseline),c=colors['real'],marker='D',s=32,zorder=4,label='Riemann residual')
    axs[0].axhline(0,color='0.4',ls=':',lw=1,label='Smooth envelope')
    axs[0].set(title='(a) All eight fixed evaluation blocks',xlabel='Block',ylabel=r'$10^4(T-T_0)$',xticks=np.arange(1,9))
    axs[0].legend(fontsize=8.5,loc='lower right')
    for f in cubes:
        mean=np.mean(cubes[f],axis=0)
        axs[1].hist(1e4*(mean-baseline),bins=np.linspace(-3.5,4,22),histtype='step',lw=1.8,
                    color=colors[f],label=labels_by_family[f])
    axs[1].axvline(1e4*(np.mean(real)-baseline),color=colors['real'],lw=2,label='Riemann mean')
    axs[1].axvline(0,color='0.4',ls=':',lw=1)
    ranks=', '.join(f"{labels_by_family[f]}: {agg[f]['empirical_two_sided_rank']:.2f}" for f in cubes)
    axs[1].set(title='(b) Prespecified mean across eight blocks',xlabel=r'$10^4(\overline{T}-T_0)$',ylabel='Count among 99 conditional controls')
    axs[1].legend(fontsize=9)
    axs[1].text(.02,.04,'Empirical tail ranks\n'+ranks,transform=axs[1].transAxes,fontsize=8.5,
                bbox=dict(facecolor='white',alpha=.9,edgecolor='none'))
    for ax in axs: ax.grid(alpha=.15)
    fig.supxlabel('Ranks describe the constructed controls; they are not calibrated p-values or physical evidence.',fontsize=9)
    for ext in ('png','pdf'):fig.savefig(FIG/f'arithmetic_primary_comparison.{ext}')
    plt.close(fig)
    fig,axs=plt.subplots(2,2,figsize=(11.3,7.6),layout='constrained')
    fig.suptitle('Input structure, field response and control limitations',fontsize=14,weight='semibold')
    selected=(1,2,101); labels=('Riemann block 1','IAAFT #1','Exact node spectrum #1')
    R=np.linspace(1,3,2049); w=np.sin(np.pi*(R-1)/2)**4
    for cid,label,color in zip(selected,labels,(colors['real'],colors['iaaft'],colors['exact_spectrum'])):
        axs[0,0].plot(np.arange(128),inp['eta'][cid],lw=.9,color=color,alpha=.8,label=label)
        cs=CubicSpline(inp['nodes'],inp['eta'][cid],bc_type='natural')
        axs[0,1].plot(R,.02*w*cs(R),lw=.9,color=color,label=label)
    axs[0,0].set(title='(a) Fixed first block and first controls',xlabel='Local gap index',ylabel=r'Bounded residual $\eta$')
    axs[0,0].legend(fontsize=8,ncol=1)
    axs[0,1].set(title='(b) Residual in the restoring coefficient',xlabel=r'Internal coordinate $R$',ylabel=r'$\delta M^2(R)$')
    realids=np.flatnonzero(inp['case_family']=='real')
    delta=main['C1'][realids]-main['C1'][0]
    bound=float(np.max(np.abs(delta)))
    im=axs[1,0].imshow(delta,aspect='auto',origin='lower',extent=(0,20,.5,8.5),cmap='RdBu_r',vmin=-bound,vmax=bound)
    axs[1,0].set(title='(c) First-mode response to real residuals',xlabel=r'Elapsed model time $\tau$',ylabel='Block',yticks=np.arange(1,9))
    fig.colorbar(im,ax=axs[1,0],label=r'$C_1-C_{1,0}$',shrink=.85)
    for f,offset in (('iaaft',-.12),('exact_spectrum',.12)):
        qq=[q for q in quality if q['family']==f]
        keys=['fft_power_relative_l2_no_dc','dense_windowed_fft_power_relative_l2_no_dc','dense_derivative_fft_power_relative_l2_no_dc']
        for j,key in enumerate(keys):
            vals=100*np.array([float(q[key]) for q in qq]);lo,mid,hi=np.quantile(vals,[.025,.5,.975])
            axs[1,1].errorbar(j+offset,mid,yerr=[[mid-lo],[hi-mid]],fmt='o',color=colors[f],capsize=4,
                            label=labels_by_family[f] if j==0 else None)
    axs[1,1].set(title='(d) Control matching by stage',ylabel='Power-spectrum relative L2 difference [%]',
                 xticks=[0,1,2],xticklabels=['Input nodes','Windowed drive','Drive derivative'])
    axs[1,1].legend(fontsize=9)
    axs[1,1].grid(alpha=.15)
    fig.supxlabel('Dimensionless synthetic dynamics. Actual continuous driver spectra remain unmatched.',fontsize=9)
    for ext in ('png','pdf'):fig.savefig(FIG/f'arithmetic_response_and_controls.{ext}')
    plt.close(fig)


if __name__=='__main__':
    main()
