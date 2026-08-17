r"""The size-effect measurement, both morphologies, one form.

r(L) is the apparent bending stiffness of a heterogeneous cell divided by that
of a matched phi=0 cell at the same edge and the same element size. Referring
to the control is the whole method: a cubic cell forced into a plate-like
curvature carries an extraction artefact of several percent, and it is present
in both, so it divides out. What survives is a property of the microstructure.

Fitting r(L) = r_inf [1 + (l/L)^2]:
  r_inf  the classical, size-independent value -- the intercept at 1/L^2 = 0
  l      the length scale, the edge at which the apparent stiffness is 2 r_inf

    python3 analysis/plot_size_effect.py [outdir]
"""
import collections
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as fs

fs.apply()
import matplotlib.pyplot as plt

RES = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   'results')


def load(fname, col='E_bending'):
    g = collections.defaultdict(list)
    with open(os.path.join(RES, fname), encoding='utf8') as fh:
        for r in csv.DictReader(fh):
            v = (r.get(col) or '').strip()
            if v in ('', 'MISSING', 'ERROR'):
                continue
            try:
                val, L = float(v), float(r['L'])
            except ValueError:
                continue
            if val > 0:
                g[L].append(val)
    return g


def ratio(hetf, homf):
    het, hom = load(hetf), load(homf)
    Ls = np.array(sorted(L for L in het if L in hom))
    y = np.array([np.mean(het[L]) / np.mean(hom[L]) for L in Ls])
    sd = np.array([np.std(np.array(het[L]) / np.mean(hom[L]), ddof=0)
                   for L in Ls])
    return Ls, y, sd


def fit(Ls, y):
    x = 1.0 / Ls ** 2
    s, c = np.polyfit(x, y, 1)
    r2 = 1 - ((y - (s * x + c)) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    return s, c, r2


def r2_of(Ls, y):
    x = 1.0 / Ls ** 2
    s, c = np.polyfit(x, y, 1)
    return 1 - ((y - (s * x + c)) ** 2).sum() / ((y - y.mean()) ** 2).sum()


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else '.'
    P = ratio('results_eringen.csv', 'results_eringen_homog.csv')
    Q = ratio('results_eringen_layer.csv', 'results_eringen_layer_homog.csv')
    sP, cP, r2P = fit(P[0], P[1])
    sQ, cQ, r2Q = fit(Q[0], Q[1])

    fig, ax = plt.subplots(1, 2, figsize=(12.2, 5.0))

    # ---- (a) the measurement, on its own scale for each morphology ---------
    a = ax[0]
    xg = np.linspace(0, 1.0 / P[0].min() ** 2 * 1.06, 100)
    for (Ls, y, sd), s, c, col, mk, lab in (
            (P, sP, cP, fs.ORANGE, 'o', 'pocket and channel'),
            (Q, sQ, cQ, fs.BLUE, 's', 'layered')):
        a.errorbar(1.0 / Ls ** 2, y, yerr=sd, fmt=mk, color=col, ms=7,
                   capsize=3, lw=0, elinewidth=1.2, label=lab)
        a.plot(xg, s * xg + c, color=col, lw=1.8, alpha=0.85)
        a.plot([0], [c], marker='<', color=col, ms=9, clip_on=False)
        a.annotate(r'$r_\infty=%.3f$' % c, xy=(0, c), xytext=(2.4, c - 0.055),
                   fontsize=10.5, color=col,
                   arrowprops=dict(arrowstyle='->', color=col, lw=0.9))
    a.set_xlim(0, xg[-1]); a.set_ylim(0.15, 0.80)
    a.set_xlabel(r'$1/L^{2}$   [model units$^{-2}$]')
    a.set_ylabel(r'$r(L)=E_{\rm bend}\,/\,E_{\rm bend}^{\,\phi=0}$')
    a.text(0.015, 0.965, '(a)', transform=a.transAxes, fontsize=13,
           fontweight='bold', va='top')
    a.legend(loc='center right', fontsize=9.5)

    # ---- (b) each referred to its own classical limit ----------------------
    b = ax[1]
    for (Ls, y, sd), s, c, col, mk, lab, micro, mname in (
            (P, sP, cP, fs.ORANGE, 'o', 'pocket and channel', 0.08, 'd'),
            (Q, sQ, cQ, fs.BLUE, 's', 'layered', 0.12, r'a_0')):
        b.errorbar(1.0 / Ls ** 2, y / c, yerr=sd / c, fmt=mk, color=col, ms=7,
                   capsize=3, lw=0, elinewidth=1.2, label=lab)
        b.plot(xg, (s * xg + c) / c, color=col, lw=1.8, alpha=0.85)
        l = np.sqrt(max(s, 0) / c)
        txt = (r'$l=%.2f\,%s$' % (l / micro, mname) if r2_of(Ls, y) > 0.5
               else r'no trend: $R^2=%.2f$, $l<%.2f\,%s$'
               % (r2_of(Ls, y), 0.58, mname))
        b.text(0.97, 0.42 if col == fs.BLUE else 0.30, txt,
               transform=b.transAxes, ha='right', fontsize=10.5, color=col)
    b.axhline(1.0, color='0.45', lw=1.0, ls=':')
    b.text(xg[-1] * 0.02, 1.02, 'classical limit', fontsize=9, color='0.4',
           ha='left', va='bottom')
    b.set_xlim(0, xg[-1]); b.set_ylim(0.92, 1.95)
    b.set_xlabel(r'$1/L^{2}$   [model units$^{-2}$]')
    b.set_ylabel(r'$r(L)\,/\,r_\infty$')
    b.text(0.015, 0.965, '(b)', transform=b.transAxes, fontsize=13,
           fontweight='bold', va='top')
    b.legend(loc='lower right', fontsize=9.5)

    fig.tight_layout()
    p = os.path.join(out, 'fig_size_effect.png')
    fig.savefig(p, dpi=175)
    print('wrote %s' % p)
    print('  pocket : r_inf %.4f  slope %+.5f  R2 %.4f' % (cP, sP, r2P))
    print('  layered: r_inf %.4f  slope %+.5f  R2 %.4f  l = %.2f a_0'
          % (cQ, sQ, r2Q, np.sqrt(sQ / cQ) / 0.12))


if __name__ == '__main__':
    main()
