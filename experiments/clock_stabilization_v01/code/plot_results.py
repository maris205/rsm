#!/usr/bin/env python3
"""Render archived conditional results; the right panel is a coordinate extrapolation."""
from pathlib import Path
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
res = ROOT / 'results'
inp = json.loads((res / 'inputs.json').read_text())
summary = json.loads((res / 'summary.json').read_text())
with (res / 'scan.csv').open() as stream:
    rows = [r for r in csv.DictReader(stream) if r['phase'] == 'zero']
with (res / 'conditional_frontiers.csv').open() as stream:
    front = [r for r in csv.DictReader(stream) if float(r['mKK_eV']) == 1e14]
x = np.log10(inp['masses_eV'])
y = np.log10(np.array(inp['vector_masses_eV']) / 1e12)
color = np.log10(np.array([float(r['known_axis_feedback_without_valley']) for r in rows]).reshape(len(x), len(y)).T)
passed = np.array([r['necessary_conditions_pass'] == 'True' for r in rows]).reshape(len(x), len(y)).T
plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False})
fig, axs = plt.subplots(1, 2, figsize=(14, 5.3), gridspec_kw={'width_ratios': [1.2, 1]})
ax = axs[0]
im = ax.pcolormesh(x, y, color, shading='nearest', cmap='viridis', vmin=-8, vmax=-1)
cb = fig.colorbar(im, ax=ax, pad=.025)
cb.set_label(r'$\log_{10} R_{\rm selected,axis}$')
ax.contour(x, y, passed.astype(float), levels=[.5], colors=['#00e7bb'], linewidths=1.8)
ax.plot(np.log10([float(r['local_curvature_mG_max_eV']) for r in front]), y,
        color='white', lw=2, ls='--', label='Local curvature boundary')
ax.plot(np.log10([float(r['auxiliary_mG_max_eV']) for r in front]), y,
        color='#ffb94e', lw=2, label='Adopted auxiliary hierarchy')
ax.axhspan(-1, 0, color='gray', alpha=.5, hatch='//')
ax.axhspan(1, 2, color='gray', alpha=.5, hatch='//')
ax.scatter([-6], [0], s=65, color='#ff7474', edgecolor='black', zorder=5)
ax.text(-11.6, .45, 'Outlined region passes\nthe stated local conditions', color='white')
ax.text(-11.6, -0.65, r'$M_Q/M_V>0.1$', color='white')
ax.text(-11.6, 1.55, r'$M_V/m_{\rm KK}>0.1$', color='white')
ax.set(xlim=(x[0], x[-1]), ylim=(-1, 2), xlabel=r'$\log_{10}(m_G/{\rm eV})$',
       ylabel=r'$\log_{10}(M_V/{\rm TeV})$', title=r'Finite-interval conditions: $m_{\rm KK}=100$ TeV')
ax.legend(loc='lower right', fontsize=8, framealpha=.9, facecolor='#303d46', labelcolor='white')

ax = axs[1]
p = inp['parent_parameters']
sample = summary['sample_aligned']
chi = np.linspace(inp['chi_first'], sample['chi_critical_leading'] + 2, 700)
U = inp['Umax_eV4'] * np.exp(-2 * (chi - inp['chi_first']))
negative = 24 * sample['A_5d'] * sample['mG_eV'] ** 2
mass_ratio = (4 * (U + p['rho']) / (3 * p['P']) + 12 * U / sample['Lambda_eV'] ** 2) / negative - 1
old = chi <= inp['chi_last']
ax.plot(chi[old], mass_ratio[old], color='#147d64', lw=2.5, label='Tested coordinate interval')
ax.plot(chi[~old], mass_ratio[~old], color='#147d64', lw=2.5, ls='--', label='Fixed-parameter extrapolation')
ax.axvspan(inp['chi_first'], inp['chi_last'], color='#147d64', alpha=.08)
ax.axhline(0, color='black', lw=.8)
ax.axhline(-1, color='#b24e61', lw=1.2, ls=':', label=r'$U\to0$ limit (approximately)')
ax.axvline(sample['chi_critical_leading'], color='#b24e61', lw=1, ls=':')
ax.set_yscale('symlog', linthresh=.2)
ax.set_yticks([-1, 0, 1, 1e2, 1e4])
ax.set_yticklabels([r'$-1$', r'$0$', r'$1$', r'$10^2$', r'$10^4$'])
ax.set(xlim=(chi[0], chi[-1]), ylim=(-1.4, 2e5), xlabel=r'Clock coordinate $\chi$ (not a new cosmic time solution)',
       ylabel=r'$m_y^2/(24 A m_G^2)$', title=r'Benchmark: $m_G=10^{-6}$ eV, $M_V=1$ TeV')
ax.text(.07, .30, 'The restoring term falls with U.\nStability is lost in this extrapolation.', transform=ax.transAxes, fontsize=10)
ax.legend(loc='upper right', fontsize=8)
ax.grid(alpha=.13)
fig.suptitle('A shared-current quartic can stabilize a finite interval, with companion matter feedback', fontsize=14)
fig.text(.5, .012, r'Declared matching: $g_Z=\sqrt{2}$, $\sigma=1$. Selected EFT terms only; no observational allowed region or new cosmic trajectory.', ha='center', fontsize=9, color='#555555')
fig.tight_layout(rect=(0, .04, 1, .93))
for suffix in ('png', 'pdf'):
    fig.savefig(ROOT / 'figures' / ('finite_interval_stabilization.' + suffix), dpi=180, bbox_inches='tight')
