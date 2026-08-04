#!/usr/bin/env python3
"""Figure: separate the level from the shape in the field comparison.

Section 4.3.2 corrects the model onto the vibrating-beam data with one scalar,
the vibrating-beam-effective matrix factor, and then improves the depth-RMS
misfit by changing the imposed salinity profile. Those two moves are easy to
confuse, and the confusion matters: if the salinity profile were only
compensating for a badly chosen scalar, the improvement would say nothing about
salinity.

The separation is straightforward. A scalar can only slide a profile up or
down, so whatever misfit survives after each column has been given its own
best-fit scalar is shape, and the rest is level. Panel (b) reports both parts.

The point of panel (a) is that no scalar fixes the shape: the adopted factor
anchors the cold surface, matching the column mean would take a different one,
and neither can bend the interior plateau onto a monotonic fit.

Run from results/:  python3 ../analysis/plot_field_decomp.py
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

# Marchenko (2024) first-year vibrating-beam fit, as used throughout.
M_EBOT, M_M, M_N = 1.67, 2.63, 0.5

CAMPAIGNS = [('MSEED', 'C-shape', fs.BLUE),
             ('NSEED', 'monotonic', fs.GREEN),
             ('N2SEED', 'steep monotonic', fs.VERM)]


def prof(pref, path='results_fieldseeds.csv'):
    g = defaultdict(list)
    for r in csv.DictReader(open(path)):
        v = r.get('E_eff')
        if not v or v == 'MISSING' or float(v) <= 0:
            continue
        rid = r['run_id'].rsplit('_s', 1)[0]
        if not rid.startswith(pref):
            continue
        g[rid].append(float(v) / 1e9)
    ks = sorted(g, key=lambda k: int(''.join(c for c in k.rsplit('_z', 1)[-1]
                                            if c.isdigit())))
    return (np.array([st.mean(g[k]) for k in ks]),
            np.array([st.pstdev(g[k]) for k in ks]))


def main():
    z = np.arange(0.05, 1.0, 0.1)
    Em = M_EBOT * ((M_M - 1) * (1 - z) ** M_N + 1.0)

    fig, ax = plt.subplots(1, 2, figsize=(12.2, 5.0))

    # ---- (a) the three profiles against the field fit -------------------
    a = ax[0]
    a.plot(Em, z, 'k--', lw=2.0, label='Marchenko 2024 fit')
    rows = []
    for pref, lab, col in CAMPAIGNS:
        E, s = prof(pref)
        a.errorbar(E, z, xerr=s, marker='o', ms=5, capsize=2, color=col, label=lab)
        rel = E / Em
        rms_as = 100 * np.sqrt(np.mean((rel - 1) ** 2))
        f = np.sum(rel) / np.sum(rel ** 2)          # least-squares rescale
        rms_shape = 100 * np.sqrt(np.mean((f * rel - 1) ** 2))
        rows.append((lab, col, rms_as, rms_shape, f))
    fs.depth_axis(a)
    a.set_xlabel("Effective Young's modulus [GPa]")
    a.set_title('(a) recalibrated columns vs the field fit')
    a.legend(loc='lower left', fontsize=10.5, framealpha=0.95)

    # ---- (b) how much of the misfit is level, how much is shape ---------
    b = ax[1]
    x = np.arange(len(rows))
    w = 0.6
    shape = [r[3] for r in rows]
    level = [r[2] - r[3] for r in rows]
    b.bar(x, shape, w, color=[r[1] for r in rows], label='shape (survives rescaling)')
    b.bar(x, level, w, bottom=shape, color='0.80', edgecolor='0.45',
          label='level (a scalar can remove)')
    for i, r in enumerate(rows):
        b.text(i, r[2] + 0.4, '%.1f%%' % r[2], ha='center', fontsize=11)
        b.text(i, r[3] / 2, '%.1f%%' % r[3], ha='center', va='center',
               fontsize=11, color='white')
    b.set_xticks(x)
    b.set_xticklabels([r[0] for r in rows])
    b.set_ylabel('depth-RMS misfit [%]')
    b.set_title('(b) level and shape parts of the misfit')
    b.legend(loc='upper right', fontsize=10.5)
    b.set_ylim(0, max(r[2] for r in rows) * 1.30)

    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig('field_decomposition.%s' % ext, dpi=200)
    print('wrote field_decomposition.{png,pdf}')
    print('%-18s %10s %10s %10s' % ('column', 'RMS', 'shape', 'best f'))
    for lab, col, r_as, r_sh, f in rows:
        print('%-18s %9.1f%% %9.1f%% %10.3f' % (lab, r_as, r_sh, f))


if __name__ == '__main__':
    main()
