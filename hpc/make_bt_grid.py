#!/usr/bin/env python3
r"""Is the bridge factor a function of b alone? A grid in (b, t) says.

The closure writes the bridge factor as g = b^n with n constant. LAYERB tested
that by sweeping b at fixed SLAB brine fraction and the exponent came out 0.93,
1.03, 0.85 at b = 0.10, 0.20, 0.30 -- non-monotonic, thirty times the seed
scatter, and on no power law at all.

That reading is incomplete, because the sweep was not a pure sweep in b. The
layer thickness is tied to the other two by

    phi_slab = t (1 - b),

so holding phi_slab while raising b THICKENS the layer: t ran 0.167, 0.188,
0.214 across those three cells, a 29% change. The two push opposite ways -- more
bridge stiffens the cell, a thicker brine layer softens it -- so the apparent
b-dependence of n is two geometric variables moving at once.

Only two of the three are independent, so no deck can vary b at fixed phi AND
fixed t. What a deck CAN do is sample the (b, t) plane properly instead of
walking one line across it. The closure walks the Assur line, b = 1 -
sqrt(phi/phi_0); LAYERB walks the fixed-phi_slab line; the two disagree about
the exponent, which is what one expects if E is a surface and b^n is a slice.

This deck is a 3x3 grid, b in (0.10, 0.20, 0.30) crossed with t in (0.14, 0.18,
0.22), two seeds each. slab_vof is set to t(1-b) so each cell realises its
intended pair. If a single exponent describes the surface, n extracted at fixed
t will be constant along b -- and the earlier result was an artefact of the path.
If it is not, the bridge factor needs two arguments and the closure's g(phi)
cannot be recovered from b alone.

SIZING. These cells are far cheaper than the RAMP ones and deliberately so.
Every layer here is thick -- t = 0.14 to 0.22 against 0.11 to 0.15 in RAMP -- so
mesh_for's own rule returns a coarse element, and the drained convergence point
measured by rve_layermesh is 1.6 elements across rather than the 2.5 a
mixed-drainage deck needs. The result is 72k-96k cells to the edge against
140k-350k, small enough to generate and solve on a workstation.

The usual objection to shrinking a cell does not apply here. RAMP and SUBC had
to keep L = 0.5 because their numbers are compared against LCOL cells at that
size, so a size bias that would have cancelled becomes differential. This deck
is compared only against itself: every point on the grid carries the same size
bias, and the question is the SHAPE of the surface, not its level.

Drained only, and only the x load case is solved -- n is a drained transverse
quantity and nothing here reads E_z.

    python3 make_bt_grid.py <outdir>
"""
import csv
import os
import sys

from make_layer_decks import BASE, COLS, K, LM_MAX, LM_MIN

# Drained convergence, measured: rve_layermesh puts the drained response within
# 1% by 1.47 elements across the brine layer and 0.3% by 2.21. 1.6 sits past the
# first and short of the second.
ELEM_ACROSS = 1.6

B_VALUES = (0.10, 0.20, 0.30)
T_VALUES = (0.14, 0.18, 0.22)
SEEDS = (1, 2)

L = 0.50
N_SLABS = 4


def slab_vof(t, b):
    """The deck's brine fraction that realises this (t, b) pair."""
    return t * (1.0 - b)


def layer_thickness(t):
    """Layer thickness in model units. t is the fraction of the cell it spans,
    shared over n_slabs layers, so the individual layer is t*L/n_slabs."""
    return t * L / N_SLABS


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else '.'
    os.makedirs(out, exist_ok=True)
    p = os.path.join(out, 'rve_btgrid.csv')

    rows = []
    with open(p, 'w', newline='', encoding='utf8') as fh:
        w = csv.writer(fh)
        w.writerow(COLS)
        for t in T_VALUES:
            for b in B_VALUES:
                sv = slab_vof(t, b)
                th = layer_thickness(t)
                lm = min(max(th / ELEM_ACROSS, LM_MIN), LM_MAX)
                for s in SEEDS:
                    r = dict(BASE)
                    r['run_id'] = 'BT_t%03d_b%03d_s%d' % (round(t * 1000),
                                                          round(b * 100), s)
                    r['slab_vof'] = '%.4f' % sv
                    r['bridge_fraction'] = '%.4f' % b
                    r['K_inclusion'] = K['drn']
                    r['L_mesh'] = '%.4f' % lm
                    r['min_distance'] = '0.005'
                    w.writerow([r[c] for c in COLS])
                rows.append((t, b, sv, th, lm))

    n = len(rows) * len(SEEDS)
    print('wrote %s  (%d cells: %d conditions x %d seeds)'
          % (p, n, len(rows), len(SEEDS)))
    print('  %-7s %-7s %-9s %-9s %-9s %-8s %s'
          % ('t', 'b', 'slab_vof', 'thickness', 'L_mesh', 'across', 'elements'))
    for t, b, sv, th, lm in rows:
        print('  %-7.2f %-7.2f %-9.4f %-9.5f %-9.5f %-8.1f %.0fk'
              % (t, b, sv, th, lm, th / lm, (L / lm) ** 3 / 1e3))

    print()
    print('What the grid can return:')
    print('  n constant along b at fixed t  -> b^n survives; the LAYERB result')
    print('                                    was the path, not the form')
    print('  n still varying along b        -> the bridge factor needs both')
    print('                                    arguments and g(phi) cannot be')
    print('                                    recovered from b alone')
    return 0


if __name__ == '__main__':
    sys.exit(main())
