#!/usr/bin/env python3
"""Scientific exports of the conditional transfer and local-curvature tests."""
from pathlib import Path
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]


def main():
    with (ROOT / 'results/scan.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    summary = json.loads((ROOT / 'results/summary.json').read_text())
    inputs = json.loads((ROOT / 'results/inputs.json').read_text())
    masses = np.array(inputs['masses_eV'])
    kk = np.array(inputs['KK_eV'])
    selected = [r for r in rows if int(r['kappa']) == 1]
    transfer = np.array([float(r['selected_c_2loop_force_ratio']) for r in selected]).reshape(len(masses), len(kk)).T
    plt.rcParams.update({'font.size': 10, 'font.family': 'DejaVu Sans',
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.8), gridspec_kw={'width_ratios': [1.22, 1]})
    mesh = axes[0].pcolormesh(np.log10(masses), np.log10(kk / 1e12), np.log10(transfer),
                             shading='auto', cmap='viridis', vmin=-50, vmax=6, rasterized=True)
    colorbar = fig.colorbar(mesh, ax=axes[0], pad=.02)
    colorbar.set_label(r'$\log_{10} R_{Q,\mathrm{selected}}$')
    force_line = np.log10(summary['force_product_bound_eV2'] / kk)
    stability_line = np.log10(summary['kappa1_stability_product_bound_eV2'] / kk)
    axes[0].plot(force_line, np.log10(kk / 1e12), color='#ff9b40', lw=2.2)
    axes[0].plot(stability_line, np.log10(kk / 1e12), color='white', lw=2.2, linestyle='--')
    axes[0].axhline(0, color='#cce8f3', lw=1.5, linestyle=':')
    axes[0].axhspan(-1, 0, color='grey', alpha=.33)
    axes[0].scatter([-18, 0], [0, 0], marker='o', s=46, facecolors=['#69d4ea', '#f58b85'], edgecolors='black', zorder=6)
    axes[0].text(-23.3, 1.62, 'Positive local curvature\nto the left of white line', color='white', fontsize=8.5)
    axes[0].text(-23.3, -.75, r'$M_Q/m_{\rm KK}>0.1$ (excluded by convention)', color='white', fontsize=8)
    axes[0].set(xlim=(-24, 2), ylim=(-1, 2), xlabel=r'$\log_{10}(m_G/\mathrm{eV})$',
                ylabel=r'$\log_{10}(m_{\rm KK}/\mathrm{TeV})$', title=r'Conditional clock extension: $\kappa=1$')
    axes[0].legend(handles=[Line2D([], [], color='#ff9b40', lw=2, label='Selected force = 10%'),
                            Line2D([], [], color='white', lw=2, linestyle='--', label='Local curvature changes sign')],
                   facecolor='#303d5b', labelcolor='white', framealpha=.85, loc='lower right', fontsize=7.9)
    colors = {-1: '#cc7a29', 0: '#267e54', 1: '#9b2743'}
    for kappa in (-1, 0, 1):
        r = [r for r in rows if int(r['kappa']) == kappa and np.isclose(float(r['mKK_eV']), 1e12, rtol=1e-12)]
        y = np.array([float(item['min_my2_over_Href2']) for item in r])
        label = {1: r'$\kappa=1$: same-sign extension', 0: r'$\kappa=0$: no clock term', -1: r'$\kappa=-1$: sign control'}[kappa]
        axes[1].plot(np.log10(masses), y, color=colors[kappa], lw=2, label=label)
    axes[1].set_yscale('symlog', linthresh=10, linscale=.7)
    tick_powers = (35, 25, 15, 5)
    axes[1].set_yticks([-10. ** x for x in tick_powers] + [0.] + [10. ** x for x in reversed(tick_powers)])
    axes[1].set_yticklabels([r'$-10^{%d}$' % x for x in tick_powers] + ['0'] + [r'$10^{%d}$' % x for x in reversed(tick_powers)])
    axes[1].axhline(0, color='black', lw=.8)
    axes[1].axvline(np.log10(summary['kappa1_stability_product_bound_eV2'] / 1e12), color='#555555', lw=1, linestyle=':')
    axes[1].set(xlim=(-24, 2), xlabel=r'$\log_{10}(m_G/\mathrm{eV})$',
                ylabel=r'$\min(V_{yy}/G_{yy})/H_{\rm ref}^2$', title=r'Local clock curvature at $m_{\rm KK}=1$ TeV')
    axes[1].legend(frameon=False, loc='upper left', fontsize=8)
    axes[1].grid(alpha=.15)
    axes[1].annotate('1 eV example:\nnegative curvature,\nlocal e-fold time 20.4 s',
                     xy=(0, summary['sample_1eV_1TeV']['min_my2_over_Href2']), xytext=(-13, -1e20),
                     arrowprops={'arrowstyle': '->', 'color': '#555555'}, fontsize=9)
    fig.suptitle('Small charged-sector feedback does not guarantee a stable clock direction', fontsize=13)
    fig.text(.5, -.018, 'Finite 5D sector + conditional 4D embedding; local slices, no new cosmological solution or observational fit.',
             ha='center', fontsize=9, color='#555555')
    fig.tight_layout()
    for extension in ('png', 'pdf'):
        fig.savefig(ROOT / 'figures' / ('transfer_and_clock_curvature.' + extension), dpi=200, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    main()
