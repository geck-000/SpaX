#!/usr/bin/env python3
"""Figure: the bending size effect against a matched no-microstructure baseline.

Panel (a) puts the channelled sweep next to homogeneous cells solved at the SAME
six cell sizes. Those cells contain no inclusions and no channels, so any size
dependence they show is the cube-versus-plate kinematics of the bending
extraction plus discretisation -- it cannot be a material length scale. The
baseline rises 4.8% across the sweep, larger than the 2.9% trend in the
channelled cells it was supposed to be contaminating.

Panel (b) regresses the apparent modulus on 1/L^2, the form both nonclassical
families take,

    E_app/E_inf = 1 + 12 l^2 / L^2      (III.A gradient / couple stress)
    E_inf/E_app = 1 + (e0a)^2 / L^2     (III.B Eringen integral nonlocal)

before and after subtracting that baseline point-by-point. Uncorrected the slope
is already not significant (p=0.14); corrected it is flat (p=0.82, R^2=0.002).

Run from results/:  python3 ../analysis/plot_sizeeffect.py
"""
import csv
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as fs
from scipy import stats

fs.apply()
import matplotlib.pyplot as plt

D = 0.08                     # mean inclusion diameter


def load(path):
    g = defaultdict(list)
    for r in csv.DictReader(open(path)):
        try:
            v = float(r['E_bending'])
        except (ValueError, KeyError, TypeError):
            continue
        if v > 0:
            g[float(r['L'])].append(v / 1e9)
    return g


def main():
    ch = load('results_eringen.csv')
    ho = load('results_eringen_homog.csv')
    Ls = sorted(ch)
    href = {L: float(np.mean(v)) for L, v in ho.items()}
    ref = href[max(href)]                      # widest cell = the asymptote
    rel = {L: href[L] / ref for L in href}

    ld = np.array([L / D for L in Ls])
    mean = np.array([np.mean(ch[L]) for L in Ls])
    sd = np.array([np.std(ch[L], ddof=0) for L in Ls])
    corr = np.array([np.mean(ch[L]) / rel[L] for L in Ls])
    csd = np.array([np.std(ch[L], ddof=0) / rel[L] for L in Ls])

    fig, ax = plt.subplots(1, 2, figsize=(12.4, 4.9))

    # ---- (a) the three curves ------------------------------------------
    a = ax[0]
    a.errorbar(ld, mean, yerr=sd, marker='o', color=fs.BLUE, capsize=3,
               label='channelled RVE (raw)')
    a.errorbar(ld, corr, yerr=csd, marker='s', color=fs.GREEN, capsize=3,
               label='channelled, baseline-corrected')
    a.set_xlabel(r'cell size $L/d$')
    a.set_ylabel(r'$E_{\mathrm{bending}}$  [GPa]')
    a.set_title('(a) apparent bending modulus')
    a.legend(loc='lower right', frameon=True, fontsize=11.5)

    # baseline on a twin axis: it is a different material (no inclusions), so
    # plotting it on the same GPa scale would invite a false comparison.
    b = a.twinx()
    rise = 100 * (max(rel.values()) / min(rel.values()) - 1)
    b.plot(ld, [rel[L] for L in Ls], marker='^', ls='--', color=fs.VERM,
           label=('homogeneous baseline:\n'
                  'no microstructure, yet rises %.1f%%' % rise))
    b.set_ylabel('baseline, relative to widest cell', color=fs.VERM)
    b.tick_params(axis='y', colors=fs.VERM)
    b.grid(False)
    b.legend(loc='upper left', frameon=True, fontsize=11)


    # ---- (b) regression on 1/L^2, before and after ----------------------
    c = ax[1]
    for tag, colour, marker, use_corr in (('raw', fs.BLUE, 'o', False),
                                          ('baseline-corrected', fs.GREEN, 's', True)):
        x, y = [], []
        for L in Ls:
            f = rel[L] if use_corr else 1.0
            for v in ch[L]:
                x.append(1.0 / L ** 2)
                y.append(v / f)
        x, y = np.array(x), np.array(y)
        r = stats.linregress(x, y)
        c.plot(x, y, marker, color=colour, alpha=0.45, ms=5, ls='none')
        xx = np.linspace(0, max(x) * 1.05, 50)
        c.plot(xx, r.intercept + r.slope * xx, '-', color=colour, lw=2.2,
               label=r'%s: slope $%+.3f\pm%.3f$, $p=%.2f$'
                     % (tag, r.slope, r.stderr, r.pvalue))
    c.set_xlabel(r'$1/L^2$')
    c.set_ylabel(r'$E_{\mathrm{bending}}$  [GPa]')
    c.set_title(r'(b) regression on $1/L^{2}$')
    c.legend(loc='lower right', frameon=True, fontsize=11)
    c.annotate('a couple-stress length scale\nrequires a positive slope here',
               xy=(0.55, 0.90), xycoords='axes fraction', fontsize=11.5,
               color=fs.BLACK, ha='left', va='top')

    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig('rve_sizeeffect_baseline.%s' % ext, dpi=200)
    print('wrote rve_sizeeffect_baseline.{png,pdf}')

    print('\n%5s %9s %9s %9s' % ('L/d', 'raw', 'baseline', 'corrected'))
    for L in Ls:
        print('%5.0f %9.3f %9.4f %9.3f'
              % (L / D, np.mean(ch[L]), rel[L], np.mean(ch[L]) / rel[L]))


if __name__ == '__main__':
    main()
