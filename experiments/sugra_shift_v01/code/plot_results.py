#!/usr/bin/env python3
"""Plot archived SUGRA test results; does not solve or fit a model."""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'figures'
COLORS = {'global': '#666b73', 'canonical': '#c46527', 'shift': '#007f8b'}
LABELS = {'global': 'Global reference', 'canonical': 'Canonical K', 'shift': 'Shift K'}


def read(name):
    with (ROOT / 'results' / name).open() as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([r[k] if r[k] != '' else np.nan for r in rows], dtype=str if k in
            ('case', 'model') else float) for k in rows[0]}


def select(table, mask):
    return {k: v[mask] for k, v in table.items()}


def save(fig, name):
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'{name}.{ext}', dpi=200, bbox_inches='tight')
    plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size': 11, 'axes.titlesize': 12,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'legend.frameon': False, 'savefig.facecolor': 'white'})
    b = read('trajectories.csv')
    m = read('linear_modes.csv')
    fig, ax = plt.subplots(1, 2, figsize=(11.2, 4.1), layout='constrained')
    for model in ('global', 'canonical', 'shift'):
        r = select(b, b['case'] == model + '_axis')
        kw = dict(color=COLORS[model], label=LABELS[model],
                  lw=3.5 if model == 'global' else 2,
                  ls='--' if model == 'shift' else '-',
                  zorder=4 if model == 'shift' else 2)
        ax[0].plot(r['z'], r['u'], **kw)
        ax[1].plot(r['z'], 1e6 * r['delta_alpha'], **kw)
    ax[0].set(ylabel=r'$\chi-\chi_i$', title='(a) Matched physical initial velocity')
    ax[1].set(ylabel=r'$10^6\Delta_{\rm th}$ (local threshold proxy)',
              title='(b) Charged mass threshold contribution')
    ax[1].axhline(0, color='#bbbbbb', lw=.8)
    ax[0].legend(loc='upper left')
    ax[1].text(.04, .40, 'Shift: +17.878 ppm\nCanonical: -22.024 ppm at z = 4.2',
               transform=ax[1].transAxes, fontsize=10)
    for a in ax:
        a.set(xlabel='Redshift z (evolution toward the right)', xlim=(4.2, 0))
        a.grid(alpha=.15)
    fig.suptitle('Specified K, W models | local charged threshold only | no observational fit', fontsize=12)
    save(fig, 'background_and_threshold')

    fig, ax = plt.subplots(2, 2, figsize=(11.2, 8), layout='constrained')
    for model in ('global', 'canonical', 'shift'):
        r = select(b, b['case'] == model + '_axis')
        q = select(m, (m['model'] == model) & (m['kappa'] == 0))
        ax[0, 0].plot(r['z'], 1e4*r['transverse_mass2_over_H2'],
                      color=COLORS[model], label=LABELS[model], lw=2)
        ax[0, 1].plot(q['z'], 1e4*(1-q['D0']), color=COLORS[model],
                      label=LABELS[model], lw=2)
    ax[0, 0].set(ylabel=r'$10^4 m_y^2/H^2$', title='(a) Positive curvature, very light field')
    ax[0, 0].legend()
    ax[0, 1].set(ylabel=r'$10^4(1-D_0)$', title='(b) Homogeneous displacement barely decays')
    ax[0, 1].text(.04, .78, 'Shift retains 99.9871%\nof its initial displacement',
                  transform=ax[0, 1].transAxes, fontsize=10)
    colors = ['#007f8b', '#6650a4', '#c46527', '#54722b']
    for k, color in zip((0., 1., 10., 100.), colors):
        r = select(m, (m['model'] == 'shift') & (m['kappa'] == k))
        label = rf'$k_{{\rm com}}/H_{{\rm ref}}={k:g}$'
        ax[1, 0].plot(r['z'], r['D0'], color=color, label=label, lw=1.3)
        ax[1, 1].semilogy(r['z'], r['energy0']/r['energy0'][0],
                          color=color, label=label, lw=1.6)
    ax[1, 0].set(ylabel=r'$D_0$', title='(c) Shift K: finite-wavelength modes')
    ax[1, 1].set(ylabel=r'$\mathcal{E}_0(N)/\mathcal{E}_0(N_i)$',
                  title='(d) Shift K: mode energy decreases')
    ax[1, 0].legend(fontsize=9, ncol=2, loc='lower left')
    for a in ax.flat:
        a.set(xlabel='Redshift z (evolution toward the right)', xlim=(4.2, 0))
        a.grid(alpha=.15)
    fig.suptitle('Classical transverse sector | finite interval | not a rapid attractor', fontsize=13)
    save(fig, 'transverse_modes')
    print('Saved 2 figure pairs from archived CSVs.')


if __name__ == '__main__':
    main()
