# -*- coding: utf-8 -*-
"""The drained/undrained bracket, and which end of it the column sits at.

A periodic cell filled with brine at K = 2.2 GPa returns the UNDRAINED modulus:
the pore fluid is sealed in, so it resists compression and the cell is stiff.
Release K and the fluid carries nothing, which is the DRAINED modulus. Real ice
lies between, and where it lies is set by whether brine can leave the loaded
region in the time the load is applied.

Two things make this a prediction rather than a knob.

First, the WIDTH of the bracket is morphology, not a parameter. A spherical
pocket softens the matrix mainly in shear, and shear is untouched by the fill's
bulk modulus, so drained and undrained differ by ~3% (Mori-Tanaka, confirmed
against the cells). A cell-spanning brine layer is fully confined, so its normal
compliance is bulk-dominated and the two limits differ by up to 6.9x (measured
this session, at fixed geometry). Cold ice holds isolated pockets and warm ice
holds connected layers, so the bracket is narrow at the top of the column and
wide at the base without anyone choosing that.

Second, WHICH end applies follows from permeability, and sea ice permeability
switches on at the percolation threshold -- Golden et al. (1998), and Pringle
et al. (2009) for the single-crystal value. Above it brine drains far faster
than a flexural test lasts; below it the ice is sealed. The transition is
therefore placed by the threshold, not fitted to the moduli.
"""
import csv
import os
import re
import sys

import numpy as np

E_ICE = 9.37
PHI_C = 0.05            # Golden rule of fives; Pringle 4.6 +/- 0.7% vertical

# Measured at fixed geometry this session (b = 0.03, phi = 0.179, two layers).
LAYER_UNDRAINED, LAYER_DRAINED = 4.605, 0.671

# brine transport properties, for the Deborah number
MU_BRINE = 2.0e-3       # Pa s
K_FILL = 2.2e9          # Pa


def drainage_time(k_perm, L_path, phi):
    """Poroelastic diffusion time over L_path.

    c = k*M/(mu*phi) with M the fluid storage modulus, taken as K_fill. Order
    of magnitude only -- permeability spans decades in the field -- which is
    all that is needed, because the answer is not marginal.
    """
    c = k_perm * K_FILL / (MU_BRINE * phi)
    return L_path ** 2 / c


def main():
    print('=' * 68)
    print('IS A FLEXURAL TEST DRAINED OR UNDRAINED?')
    print('=' * 68)
    print('poroelastic drainage time, for a 0.05 m path to the nearest channel')
    print('%14s %16s %14s' % ('k (m^2)', 'drain time (s)', 'vs 30 s test'))
    for k in (1e-13, 1e-12, 1e-11, 1e-10, 1e-9):
        t = drainage_time(k, 0.05, 0.15)
        print('%14.0e %16.4g %14s' % (
            k, t, 'DRAINED' if t < 3.0 else 'undrained'))
    print('\nSea ice permeability above the percolation threshold is 1e-12 to')
    print('1e-9 m^2, so brine leaves the loaded region in well under a second')
    print('against tests lasting tens of seconds. Permeable ice is DRAINED,')
    print('and the conclusion does not depend on pinning k closely.')

    print('\n' + '=' * 68)
    print('BRACKET WIDTH IS SET BY MORPHOLOGY')
    print('=' * 68)
    print('%-26s %11s %11s %8s' % (
        'microstructure', 'undrained', 'drained', 'ratio'))
    print('%-26s %11.3f %11.3f %8.2f' % (
        'isolated pockets (MT)', 6.110, 5.903, 6.110 / 5.903))
    print('%-26s %11.3f %11.3f %8.2f' % (
        'spanning layers (FE)', LAYER_UNDRAINED, LAYER_DRAINED,
        LAYER_UNDRAINED / LAYER_DRAINED))
    print('\nSo drainage is nearly irrelevant to a pocket and decisive for a')
    print('layer. The column changes morphology with depth, which is what')
    print('turns a materials question into a depth-dependent one.')

    columns = [
        ('Marchenko column', 'results_column.csv', 'VoF_incl_sphere', 'E_eff'),
        ('Gogolaze column', 'results_gogo_column.csv', 'phi_soft_total', 'E_x'),
    ]
    root = os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), 'results')
    for name, fn, pcol, ecol in columns:
        path = os.path.join(root, fn)
        if not os.path.exists(path):
            continue
        rows = {}
        for r in csv.DictReader(open(path)):
            m = re.search(r'z(\d+)', r['run_id'])
            if not m:
                continue
            try:
                p, e = float(r[pcol]), float(r[ecol])
            except (ValueError, TypeError, KeyError):
                continue
            rows.setdefault(int(m.group(1)), []).append((p, e / 1e9))
        if not rows:
            continue
        print('\n' + '=' * 68)
        print('%s: where the threshold falls' % name.upper())
        print('=' * 68)
        print('%8s %10s %11s %12s' % ('depth', 'phi', 'E cell GPa', 'regime'))
        for z in sorted(rows):
            p = float(np.mean([a for a, _ in rows[z]]))
            e = float(np.mean([b for _, b in rows[z]]))
            reg = 'DRAINED' if p > PHI_C else 'undrained (sealed)'
            print('%8d %10.4f %11.3f %12s' % (z, p, e, reg))
        ps = [float(np.mean([a for a, _ in rows[z]])) for z in sorted(rows)]
        n_perm = sum(1 for p in ps if p > PHI_C)
        print('\n%d of %d slices lie above the percolation threshold %.2f'
              % (n_perm, len(ps), PHI_C))
        if n_perm == len(ps):
            print('This column is permeable throughout, so the threshold does')
            print('NOT place a transition inside it: every slice is drained and')
            print('the depth dependence must come from morphology alone.')
        elif n_perm == 0:
            print('This column is sealed throughout.')
        else:
            print('The transition falls INSIDE this column, so the profile')
            print('carries a drained base and a sealed top, which is a')
            print('monotonic change with depth over a C-shaped porosity.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
