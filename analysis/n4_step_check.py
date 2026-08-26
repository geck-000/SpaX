# -*- coding: utf-8 -*-
"""Does the step at b = 0.312 survive four bridges? No.

The closure's ramp between phi_c and phi_sat was fitted to two-bridge cells. In
those cells n = ln(E_x/E_pocket)/ln(b) jumps by 0.256 across a 1.2% change in b,
which reads as a sharp material transition and is what phi_sat was derived from.

It is not a material transition. At b ~ 0.314 each of two bridges spans about
44% of the cell edge, so the pair and their periodic images make and break a
connected ice path across the plane. Dividing the same ice area into four
removes the coincidence, and with it the step: six cells spanning b = 0.388 to
0.311 -- sixteen times the range over which the two-bridge step occurs, and
crossing the split -- hold n to a range of 0.051.

Conventions follow ramp_exponent.exponents() exactly so the numbers are
comparable: phi is phi_inclusion (the packer under-reports; see the meshed-vs-
packed note), and b is the bridge_fraction the cell was BUILT with.

    python3 analysis/n4_step_check.py
"""
import csv
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

E_ICE, DRAIN = 9.37, 1.04
SPLIT = 0.312


def e_pocket(phi):
    """The drained pocket law, GPa. Eq. (2) divided by the drainage factor."""
    return E_ICE * (1.0 - 1.65 * phi) / DRAIN


def n_of(E, phi, b):
    return math.log(E / e_pocket(phi)) / math.log(b)


def load(name):
    path = os.path.join(ROOT, 'results', name)
    if not os.path.exists(path):
        return []
    rows = []
    for r in csv.DictReader(open(path)):
        if not r.get('E_x_GPa'):
            continue
        phi = float(r['phi_inclusion'])
        b = float(r['bridge_fraction'])
        E = float(r['E_x_GPa'])
        rows.append((r['run_id'], phi, b, int(r['n_bridges']), E,
                     n_of(E, phi, b)))
    return rows


def mean(v):
    return sum(v) / len(v)


def sd(v):
    m = mean(v)
    return (sum((x - m) ** 2 for x in v) / len(v)) ** 0.5   # population, ddof=0


def main():
    rows = load('results_gapcells.csv') + load('results_n4redo.csv')
    if not rows:
        print('no gap-cell or n4redo results present')
        return 1

    hdr = '%-24s %7s %7s %3s %8s %7s' % ('run_id', 'phi', 'b', 'N', 'E_x', 'n')
    for N in (2, 4):
        sel = sorted([r for r in rows if r[3] == N], key=lambda x: -x[2])
        if not sel:
            continue
        print('=' * 66)
        print('N = %d' % N)
        print(hdr)
        print('-' * 66)
        for rid, phi, b, n_br, E, n in sel:
            mark = '  <-- below the split' if N == 2 and b < SPLIT else ''
            print('%-24s %7.5f %7.4f %3d %8.4f %7.4f%s'
                  % (rid, phi, b, n_br, E, n, mark))
        print()

    hi = [r[5] for r in rows if r[3] == 2 and r[2] > SPLIT]
    lo = [r[5] for r in rows if r[3] == 2 and r[2] < SPLIT]
    n4 = [r[5] for r in rows if r[3] == 4]
    print('=' * 66)
    if hi and lo:
        print('N=2 above the split : n = %.3f  (%d cells)' % (mean(hi), len(hi)))
        print('N=2 below the split : n = %.3f  (%d cells)   STEP = %.3f'
              % (mean(lo), len(lo), mean(lo) - mean(hi)))
    if n4:
        print('N=4 across ALL b    : n = %.3f  (%d cells)   range %.3f, s.d. %.3f'
              % (mean(n4), len(n4), max(n4) - min(n4), sd(n4)))
    if hi and lo and n4:
        print()
        print('The step is %.1fx the entire N=4 spread, over 1/16 of the range in b.'
              % ((mean(lo) - mean(hi)) / (max(n4) - min(n4))))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
