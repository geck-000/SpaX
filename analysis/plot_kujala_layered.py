r"""Figure 18 redrawn with the layered basal zone and the drainage bracket.

Keeps the structure of plot_kujala.py -- profile on the left, the shape metrics
that follow from it on the right -- and adds what the layered morphology does to
both. The point of the original figure is that the disagreement is one thing
seen three ways: the base modulus, the grading parameter alpha and the neutral
plane all fail together because the computed profile is too flat. So the test of
the layered zone has to be made on all three, not on the base modulus alone.

The drained edge of the bracket is drawn HATCHED and labelled provisional. At
fixed plate spacing the drained cell's transverse modulus still falls with cell
size as n^-1.14, so it is not yet a converged material property; the campaign
that aligns bridges between layers is running and decides it. The undrained
edge passed the same test at CV 0.9% and is drawn solid.

Run from results/:  python3 ../analysis/plot_kujala_layered.py
"""
import os
import statistics as st
import sys
from collections import defaultdict
import csv

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as fs

fs.apply()
import matplotlib.pyplot as plt

from gogolaze_layered import pocket, layered, blend

K_TOP = np.array([7.18, 8.16, 8.25, 8.60])
K_BOT = np.array([0.86, 1.25, 1.56, 1.42])
K_Z0 = np.array([0.37, 0.38, 0.39, 0.38])
K_EFF_MEAN, K_EFF_SD = 4.1, 0.5

LAYER_TOP = 0.75          # where percolation switches the morphology over


def load(path, col='E_eff'):
    g = defaultdict(list)
    for r in csv.DictReader(open(path, encoding='utf8', errors='replace')):
        try:
            v = float(r[col])
        except (ValueError, KeyError, TypeError):
            continue
        if v > 0:
            g[r['run_id'].split('_s')[0]].append(v / 1e9)
    ks = sorted(g, key=lambda k: int(''.join(c for c in k if c.isdigit()) or 0))
    return (np.array([st.mean(g[k]) for k in ks]),
            np.array([st.pstdev(g[k]) for k in ks]))


def neutral(E, z):
    return float(np.trapz(E * z, z) / np.trapz(E, z))


def flex(E, z):
    z0 = neutral(E, z)
    D = np.trapz(E * (z - z0) ** 2, z)
    return float(12.0 * D / (z[-1] - z[0]) ** 3)


def main():
    m, s = load('results_column_ensemble.csv')
    zc = np.linspace(0.05, 0.95, len(m))

    # the same column on a fine grid, so the layered blend can be applied
    zf = np.linspace(0, 1, 400)
    Ep = np.interp(zf, zc, m)
    phi = np.interp(zf, [0, .29, .63, .79, .96, 1.0],
                    [0.104, 0.086, 0.128, 0.168, 0.227, 0.227])
    w = np.clip((zf - LAYER_TOP) / (1.0 - LAYER_TOP), 0.0, 1.0)
    # scale the blended layered profile so its cold end matches this column
    Edr = blend(w, phi, True) * (Ep[0] / pocket(phi)[0])
    Eun = blend(w, phi, False) * (Ep[0] / pocket(phi)[0])

    fig, ax = plt.subplots(1, 2, figsize=(12.6, 5.6))

    # ---- (a) profiles ----------------------------------------------------
    a = ax[0]
    a.fill_betweenx(zf, Edr, Eun, color=fs.SKY, alpha=0.35, zorder=1,
                    label='layered base: drainage bracket')
    a.plot(Eun, zf, color=fs.BLUE, lw=2.2, zorder=3,
           label='layered base, undrained (converged)')
    a.plot(Edr, zf, color=fs.BLUE, lw=2.0, ls=(0, (4, 3)), zorder=3,
           label='layered base, drained (provisional)')
    a.errorbar(m, zc, xerr=s, marker='o', ms=5, capsize=2, color=fs.GREEN,
               zorder=4, label='pocket column (current Fig. 18)')

    zz = np.array([0.0, 1.0])
    a.fill_betweenx(zz, [K_TOP.min(), K_BOT.min()], [K_TOP.max(), K_BOT.max()],
                    color=fs.VERM, alpha=0.16, zorder=0,
                    label='Kujala 1990, inferred linear $E(z)$')
    a.plot([K_TOP.mean(), K_BOT.mean()], zz, color=fs.VERM, lw=2.0, ls='--',
           zorder=2)
    a.errorbar([K_TOP.mean()], [0.0], xerr=[K_TOP.std(ddof=0)], marker='^',
               ms=9, color=fs.VERM, capsize=3, zorder=5)
    a.errorbar([K_BOT.mean()], [1.0], xerr=[K_BOT.std(ddof=0)], marker='v',
               ms=9, color=fs.VERM, capsize=3, zorder=5)
    fs.depth_axis(a)
    a.set_xscale('log'); a.set_xlim(0.2, 14)
    a.set_xlabel("Effective Young's modulus [GPa]")
    a.set_title('(a) computed profile vs the inferred one')
    a.legend(loc='lower left', fontsize=9.5, framealpha=0.95)

    # ---- (b) the shape metrics ------------------------------------------
    b = ax[1]
    def metrics(E):
        at = E[-1] / E[0]
        return [at, neutral(E, zf), E[0] / flex(E, zf)]
    at_k = float((K_BOT / K_TOP).mean())
    theirs = [at_k, float(K_Z0.mean()),
              3 * (1 + at_k) / (at_k ** 2 + 4 * at_k + 1)]
    pk = metrics(np.interp(zf, zc, m))
    dr, un = metrics(Edr), metrics(Eun)

    labels = [r'$\alpha=E_b/E_t$', r'$z_0/H$', r'$E_t/E_{\rm flex}$']
    x = np.arange(3)
    bw = 0.26
    b.bar(x - bw, pk, bw, color=fs.GREEN, label='pocket column')
    for i in range(3):
        lo, hi = min(dr[i], un[i]), max(dr[i], un[i])
        b.bar(x[i], hi - lo, bw, bottom=lo, color=fs.SKY, alpha=0.55,
              hatch='//', edgecolor=fs.BLUE,
              label='layered bracket' if i == 0 else None)
    b.bar(x + bw, theirs, bw, color=fs.VERM, label='Kujala 1990')
    b.set_xticks(x); b.set_xticklabels(labels, fontsize=13)
    b.set_ylabel('value')
    b.set_title('(b) the three ways the disagreement shows')
    b.legend(fontsize=10.5)
    b.grid(axis='x', alpha=0)

    for i, (p, t) in enumerate(zip(pk, theirs)):
        b.text(x[i] - bw, p + 0.03, '%.2f' % p, ha='center', fontsize=10)
        b.text(x[i] + bw, t + 0.03, '%.2f' % t, ha='center', fontsize=10)

    fig.tight_layout()
    out = sys.argv[1] if len(sys.argv) > 1 else 'kujala_layered.png'
    fig.savefig(out, dpi=170)
    print('wrote %s' % out)
    print('  pocket   alpha %.3f  z0/H %.3f  Et/Eflex %.3f' % tuple(pk))
    print('  drained  alpha %.3f  z0/H %.3f  Et/Eflex %.3f' % tuple(dr))
    print('  undrained alpha %.3f  z0/H %.3f  Et/Eflex %.3f' % tuple(un))
    print('  Kujala   alpha %.3f  z0/H %.3f  Et/Eflex %.3f' % tuple(theirs))


if __name__ == '__main__':
    main()
