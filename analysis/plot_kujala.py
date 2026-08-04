#!/usr/bin/env python3
"""Figure: the computed depth profile against the four-point bending decomposition.

Kujala et al. (1990) are the only in-situ dataset that resolves the modulus
through the thickness rather than returning one effective value per beam, which
makes them the natural comparison for a depth-resolved homogenisation. Their
Table 2 gives, for the four strain-gauged beams, a top-surface modulus, a
bottom-surface modulus and a neutral-axis position, obtained by assuming E(z)
linear and fitting two parameters to deflection and surface strain.

The comparison has to be made carefully, because the two determinations are not
the same object. Ours is a computed profile; theirs is a two-parameter fit whose
functional form is assumed. The panels separate what agrees from what does not.

(a) The profiles. Our cold-end modulus reaches the top of the range of their
    measured top-surface values with no rescaling -- it exceeds their stiffest
    of four beams by 2% and their softest by 22% -- which is the closest
    agreement in the comparison, though at the stiff edge of their range rather
    than central to it. The base does not agree at all: their inferred
    bottom-surface modulus is far softer than anything a periodic cell can be
    made to produce.

(b) What follows from the shape. Their linear profile places the neutral plane
    at z0/H = 0.37-0.39; ours, being convex, places it near 0.44-0.47 whatever
    the endpoints. The same convexity makes our profile less gradient-affected
    than a straight line, so the flexural correction E_top/E_flex is 1.17 for
    our column against ~2.1 for theirs.

Run from results/:  python3 ../analysis/plot_kujala.py
"""
import csv
import os
import statistics as st
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as fs

fs.apply()
import matplotlib.pyplot as plt

# Kujala et al. (1990), Table 2: the four beams carrying surface strain gauges.
K_TOP = np.array([7.18, 8.16, 8.25, 8.60])
K_BOT = np.array([0.86, 1.25, 1.56, 1.42])
K_Z0 = np.array([0.37, 0.38, 0.39, 0.38])
K_EFF_MEAN, K_EFF_SD = 4.1, 0.5          # all 34 beams, thickness-averaged


def load(path):
    g = defaultdict(list)
    for r in csv.DictReader(open(path)):
        try:
            v = float(r['E_eff'])
        except (ValueError, KeyError, TypeError):
            continue
        if v > 0:
            g[r['run_id'].split('_s')[0]].append(v / 1e9)
    key = lambda k: int(''.join(c for c in k if c.isdigit()) or 0)
    ks = sorted(g, key=key)
    return np.array([st.mean(g[k]) for k in ks]), np.array([st.pstdev(g[k]) for k in ks])


def neutral(E, H=1.0):
    n = len(E)
    t = H / n
    z = np.array([(i + 0.5) * t for i in range(n)])
    return float((E * t * z).sum() / (E * t).sum()) / H


def flex(E, H=1.0):
    n = len(E)
    t = H / n
    z = np.array([(i + 0.5) * t for i in range(n)])
    zb = float((E * t * z).sum() / (E * t).sum())
    D = float((E * (t ** 3 / 12.0 + t * (z - zb) ** 2)).sum())
    return 12.0 * D / H ** 3


def main():
    prof = {}
    # The ensemble column is the profile the manuscript tabulates and assembles
    # by CLT, so panel (b) has to be computed on it or the figure prints a
    # different E_top/E_flex than Section 4.3.2.
    for path, lab in (('results_column_ensemble.csv', 'C-shape column'),
                      ('results_steep_column.csv', 'steep monotonic column')):
        if os.path.isfile(path):
            m, s = load(path)
            if len(m) >= 10:
                prof[lab] = (m[:10], s[:10])
    if not prof:
        raise SystemExit('no column results found')

    z = np.linspace(0.05, 0.95, 10)
    fig, ax = plt.subplots(1, 2, figsize=(11.6, 5.0))

    # ---- (a) profiles ---------------------------------------------------
    a = ax[0]
    cols = [fs.BLUE, fs.GREEN]
    for (lab, (m, s)), c in zip(prof.items(), cols):
        a.errorbar(m, z, xerr=s, marker='o', ms=5, capsize=2, color=c, label=lab)

    # their linear profile, drawn as a band over the four beams
    zz = np.array([0.0, 1.0])
    for lo, hi, c, lab in ((K_TOP.min(), K_TOP.max(), fs.VERM, None),):
        a.fill_betweenx([0, 0], 0, 0, color=c)  # placeholder to keep colour order
    lin_lo = np.array([K_TOP.min(), K_BOT.min()])
    lin_hi = np.array([K_TOP.max(), K_BOT.max()])
    a.fill_betweenx(zz, lin_lo, lin_hi, color=fs.VERM, alpha=0.18,
                    label='Kujala 1990, inferred linear $E(z)$')
    a.plot([K_TOP.mean(), K_BOT.mean()], zz, color=fs.VERM, lw=2.0, ls='--')
    a.errorbar([K_TOP.mean()], [0.0], xerr=[K_TOP.std(ddof=0)], marker='^',
               ms=9, color=fs.VERM, capsize=3)
    a.errorbar([K_BOT.mean()], [1.0], xerr=[K_BOT.std(ddof=0)], marker='v',
               ms=9, color=fs.VERM, capsize=3)
    a.axvspan(K_EFF_MEAN - K_EFF_SD, K_EFF_MEAN + K_EFF_SD, color=fs.ORANGE,
              alpha=0.16, zorder=0)
    a.annotate('their thickness-averaged\n$E_f=4.1\\pm0.5$ GPa', xy=(K_EFF_MEAN, 0.30),
               xytext=(0.6, 0.10), fontsize=11, color=fs.ORANGE,
               arrowprops=dict(arrowstyle='->', color=fs.ORANGE, lw=1.2))
    fs.depth_axis(a)
    a.set_xlabel('Effective Young\'s modulus [GPa]')
    a.set_title('(a) computed profile vs the inferred one')
    a.legend(loc='center left', fontsize=10.5, framealpha=0.95)

    # ---- (b) what follows from the shape --------------------------------
    b = ax[1]
    labels, ours, theirs = [], [], []
    m0 = prof['C-shape column'][0] if 'C-shape column' in prof else list(prof.values())[0][0]
    labels.append(r'$\alpha=E_{b}/E_{t}$')
    ours.append(m0[-1] / m0[0]); theirs.append(float((K_BOT / K_TOP).mean()))
    labels.append(r'$z_0/H$')
    ours.append(neutral(m0)); theirs.append(float(K_Z0.mean()))
    labels.append(r'$E_{t}/E_{\mathrm{flex}}$')
    ours.append(m0[0] / flex(m0))
    at = float((K_BOT / K_TOP).mean())
    theirs.append(3 * (1 + at) / (at ** 2 + 4 * at + 1))

    x = np.arange(len(labels))
    w = 0.36
    b.bar(x - w / 2, ours, w, color=fs.BLUE, label='this work (computed, convex)')
    b.bar(x + w / 2, theirs, w, color=fs.VERM,
          label='Kujala 1990 (fitted, linear)')
    for i, (o, t) in enumerate(zip(ours, theirs)):
        b.text(i - w / 2, o + 0.03, '%.2f' % o, ha='center', fontsize=11)
        b.text(i + w / 2, t + 0.03, '%.2f' % t, ha='center', fontsize=11)
    b.set_xticks(x); b.set_xticklabels(labels)
    b.set_ylabel('value')
    b.set_title('(b) consequences of the profile shape')
    b.legend(loc='upper left', fontsize=10.5)
    b.set_ylim(0, max(max(ours), max(theirs)) * 1.28)

    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig('kujala_comparison.%s' % ext, dpi=200)
    print('wrote kujala_comparison.{png,pdf}')
    for lab, (m, s) in prof.items():
        print('  %-24s E_top=%.2f  E_base=%.2f  alpha=%.3f  z0/H=%.3f  Et/Eflex=%.2f'
              % (lab, m[0], m[-1], m[-1] / m[0], neutral(m), m[0] / flex(m)))
    print('  %-24s E_top=%.2f  E_base=%.2f  alpha=%.3f  z0/H=%.3f'
          % ('Kujala 1990', K_TOP.mean(), K_BOT.mean(),
             float((K_BOT / K_TOP).mean()), float(K_Z0.mean())))


if __name__ == '__main__':
    main()
