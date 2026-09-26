#!/usr/bin/env python3
"""Render archived results only; no recomputation of scientific quantities."""
from pathlib import Path
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'figures'
OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.size': 11, 'axes.spines.top': False,
                     'axes.spines.right': False, 'savefig.dpi': 180,
                     'axes.titleweight': 'bold', 'figure.facecolor': 'white'})

def rows(name):
    with (ROOT / 'results' / name).open() as f:
        return list(csv.DictReader(f))

def save(fig, name):
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / (name + '.' + ext), bbox_inches='tight')
    plt.close(fig)

data = [r for r in rows('scan.csv') if r['phase'] == 'zero']
x = np.unique([float(r['mG_eV']) for r in data])
y = np.unique([float(r['rH']) for r in data])
budget = np.empty((len(y), len(x)))
passed = np.empty_like(budget)
for r in data:
    i = np.searchsorted(y, float(r['rH']))
    j = np.searchsorted(x, float(r['mG_eV']))
    # A failed local-mass branch has no quadratic-valley budget. Preserve its
    # missing value; do not replace it with zero or a successful classification.
    budget[i, j] = float(r['selected_feedback_sum_ratio'] or 'nan')
    passed[i, j] = r['stated_conditions_pass'] == 'True'
fig, ax = plt.subplots(figsize=(9.2, 6.0), layout='constrained')
cmap = plt.get_cmap('magma').copy()
cmap.set_bad('#d5d9df')
mesh = ax.pcolormesh(x, y, np.log10(budget), shading='nearest',
                     cmap=cmap, norm=Normalize(-4, 8), rasterized=True)
ax.set_xscale('log'); ax.set_yscale('log')
ax.contour(x, y, passed, levels=[.5], colors=['#51e5b5'], linewidths=2.2)
ax.axhline(2, color='#84ccff', linestyle='--', linewidth=1.5)
ax.scatter([1e-6], [4], marker='*', s=210, color='white',
           edgecolor='black', linewidth=1, zorder=5)
ax.annotate('Worked example', (1e-6, 4), xytext=(1e-10, 30), color='white',
            arrowprops={'arrowstyle': '->', 'color': 'white'})
ax.text(2e-12, 2e9, 'Green outline: all stated local conditions\n'
        'Blue line: leading persistent-mass boundary\n'
        'Grey: budget undefined after local mass failure', color='white', fontsize=10)
ax.set(xlabel=r'$m_G$ [eV]', ylabel=r'$r_H=\tau P/(A\Lambda^2)$',
       title='A conditional window with persistent transverse restoration')
cb = fig.colorbar(mesh, ax=ax, pad=.025)
cb.set_label(r'$\log_{10}\,\max_\chi R_{\rm selected}$ (conservative budget)')
fig.supxlabel(r'$M_V=\Lambda=1$ TeV, $m_{\rm KK}=100$ TeV, $\sigma=1$, aligned phase'
              '\nThe example needs an unexplained coupling: '
              r'$\tau=2.93\times10^{-61}$. No observational likelihood.', fontsize=10)
save(fig, 'persistent_window')

future = rows('future_curves.csv')
fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.1), layout='constrained')
for rh, color, style, label in [(0, '#8493a5', '--', 'No persistent current'),
                               (4, '#167d9a', '-', r'Persistent current ($r_H=4$)')]:
    rr = [r for r in future if float(r['rH']) == rh]
    chi = np.array([float(r['chi']) for r in rr])
    mass = np.array([float(r['local_my2_eV2']) for r in rr])
    axes[0].plot(chi, mass / 1e-41, color=color, ls=style, lw=2, label=label)
axes[0].axhline(0, color='black', lw=.8)
axes[0].set_yscale('symlog', linthresh=2)
axes[0].set(xlim=(139.794, 155), xlabel=r'Clock coordinate $\chi$',
            ylabel=r'Local $m_y^2$ / $10^{-41}$ eV$^2$',
            title='(a) The transverse mass can remain positive')
axes[0].legend(loc='upper right', fontsize=9)
rr = [r for r in future if float(r['rH']) == 4]
chi = np.array([float(r['chi']) for r in rr])
axes[1].semilogy(chi, [float(r['moduli_budget']) for r in rr],
                 color='#167d9a', lw=2, label='Conservative selected budget')
axes[1].semilogy(chi, [float(r['signed_chi_budget']) for r in rr],
                 color='#db7542', lw=2, ls='--', label='Signed selected real-axis force')
axes[1].axhline(.1, color='black', lw=1, ls=':', label='Declared 10% threshold')
summary = json.loads((ROOT / 'results/summary.json').read_text())
limit = next(r for r in summary['future_limits'] if r['rH'] == 4 and r['phase'] == 'zero')
for key, color, text_y in [('chi_moduli_budget_10percent', '#167d9a', .006),
                          ('chi_signed_chi_budget_10percent', '#db7542', .0004)]:
    xc = limit[key]
    axes[1].plot(xc, .1, 'o', color=color, ms=6)
    axes[1].annotate(f'{xc:.3f}', (xc, .1), xytext=(xc+.8, text_y), color=color,
                      arrowprops={'arrowstyle': '-', 'color': color}, fontsize=10)
axes[1].set(xlim=(139.794, 155), ylim=(1e-5, 1e5), xlabel=r'Clock coordinate $\chi$',
            ylabel=r'Selected loop-force ratio to $2U$',
            title='(b) The quantum slope eventually dominates')
axes[1].legend(loc='upper left', fontsize=9)
for ax in axes:
    ax.grid(alpha=.15, which='major')
fig.supxlabel(r'Fixed inputs: $m_G=10^{-6}$ eV, $r_H=4$, $\sigma=1$, aligned phase. '
              r'$U\propto e^{-2\chi}$ and the mass law are prescribed.'
              '\nThe horizontal coordinate is not a cosmic date. Curves show declared sectors and matching only.', fontsize=10)
save(fig, 'persistent_future_limit')
print('Saved four figure files from frozen CSV/JSON results.')
