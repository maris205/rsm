#!/usr/bin/env python3
"""Export publication-style plots of the explicitly conditional EFT diagnostics."""
from pathlib import Path
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'figures'


def read(name):
    with (ROOT / 'results' / name).open() as stream:
        return list(csv.DictReader(stream))


def save(fig, name):
    fig.savefig(OUT / (name + '.png'), dpi=200, bbox_inches='tight')
    fig.savefig(OUT / (name + '.pdf'), bbox_inches='tight')
    plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False})
    rows = [r for r in read('coefficient_budgets.csv') if r['phase'] == 'zero']
    masses = np.array([float(r['mG_eV']) for r in rows])
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.4))
    for key, label, color, marker in [('c_10percent', r'$|c|$', '#176c9c', 'o'),
                                      ('abs_beta_10percent', r'$|\beta|$', '#cc6d23', 's'),
                                      ('abs_nu_10percent', r'$|\nu_0|$', '#6e5195', '^')]:
        axes[0].loglog(masses, [float(r[key]) for r in rows], label=label, color=color, marker=marker)
    axes[0].set(xlabel=r'$m_G$ [eV]', ylabel='One-coefficient 10% tolerance',
                title=r'Conditional boundary: $\nu_0=0$ for $c,\beta$')
    axes[0].legend(frameon=False)
    axes[0].grid(alpha=.18, which='major')
    rg = [r for r in read('rg_closure.csv') if float(r['mG_eV']) == 1 and r['phase'] == 'zero']
    for case, label, color, marker in [('c', r'$c$ alone', '#176c9c', 'o'),
                                      ('beta2', r'$\beta^2$ alone', '#cc6d23', 's')]:
        selected = [r for r in rg if r['case'] == case]
        mu = [float(r['mu_over_Mi']) for r in selected]
        axes[1].loglog(mu, [float(r['incorrectly_frozen_nu_max_force_ratio']) for r in selected],
                       '--', color=color, marker=marker, label=label + ', fixed counterterm')
    selected = [r for r in rg if r['case'] == 'c']
    axes[1].loglog([float(r['mu_over_Mi']) for r in selected],
                   [float(r['matched_max_force_ratio']) for r in selected],
                   '-', color='#267e54', linewidth=2.4, label='Matched running (both cases)')
    axes[1].set(xlabel=r'$\mu/M_i$', ylabel=r'$\max\|\Delta\nabla V\|/(2U)$',
                title=r'Same finite boundary, $m_G=1$ eV')
    axes[1].legend(frameon=False, fontsize=8.4, loc='best')
    axes[1].grid(alpha=.18, which='major')
    fig.suptitle('Local EFT sensitivities and selected one-loop scale closure', fontsize=13)
    fig.text(.5, -.015, 'Specified matching conditions; these are not observational limits or full SUGRA loop results.',
             ha='center', fontsize=9, color='#505050')
    fig.tight_layout()
    save(fig, 'wilson_budgets_and_rg')

    kernel = [r for r in read('response_kernels.csv') if r['phase'] == 'zero']
    chi = np.array([float(r['chi']) for r in kernel])
    U = np.array([float(r['U_eV4']) for r in kernel])
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.3))
    for key, label, color, style in [('c', r'$K_c$', '#176c9c', '-'),
                                    ('beta2', r'$K_{\beta^2}$', '#cc6d23', '--'),
                                    ('nu', r'$K_{\nu}$', '#6e5195', ':')]:
        k = np.array([float(r[key + '_force_chi_eV4']) for r in kernel]) / (2 * U)
        axes[0].plot(chi - chi[0], k / np.max(np.abs(k)), style, color=color, label=label, linewidth=2)
    axes[0].set(xlabel=r'$\chi-\chi_i$ (old coordinate grid)', ylabel='Signed, individually normalized response',
                title='Different operators, limited shape information')
    axes[0].legend(frameon=False)
    axes[0].grid(alpha=.18)
    result = next(r for r in json.loads((ROOT / 'results/identifiability.json').read_text())
                  if r['mG_eV'] == 1 and r['phase'] == 'zero')
    sv = result['normalized_singular_values']
    axes[1].bar([1, 2, 3], sv, color=['#176c9c', '#cc6d23', '#aaaaaa'], width=.55)
    axes[1].set(yscale='log', ylim=(1e-17, 10), xticks=[1, 2, 3], xlabel='Singular-value index',
                ylabel='Singular value (unit-norm columns)', title='Three coefficients, only two combinations')
    axes[1].text(.05, .9, r'$K_c+2K_{\beta^2}-K_\nu=0$', transform=axes[1].transAxes, fontsize=13)
    axes[1].annotate('Roundoff; analytic zero', xy=(3, sv[2]), xytext=(1.2, 1e-11),
                     arrowprops={'arrowstyle': '->', 'color': '#555555'}, fontsize=9)
    axes[1].grid(axis='y', alpha=.18)
    fig.suptitle('Response degeneracy in the selected local operator basis', fontsize=13)
    fig.text(.5, -.015, 'A fitted response cannot independently determine c, beta, and the finite local potential coefficient.',
             ha='center', fontsize=9, color='#505050')
    fig.tight_layout()
    save(fig, 'kernel_shapes_and_degeneracy')


if __name__ == '__main__':
    main()
