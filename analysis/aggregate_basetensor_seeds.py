#!/usr/bin/env python3
"""Aggregate the base-slice full-tensor replicates into an ensemble statement.

Reads elasticity_tensor_<prefix>_z95_s*.csv (one per packing, six load cases each,
written by SpaX_PostProcess.extract_elasticity_tensor) and reports the mean and
population standard deviation of the engineering constants across packings.

The question this answers: the single base cell of the depth sweep splits
E_x = 4.85 against E_y = 5.02 GPa, a 3.3% in-plane difference. Is that a
resolution limit -- too few channels in one cell to average the in-plane
directions -- or a preferred in-plane direction? check_channel_isotropy.py
settles the second: the generator is unbiased. So a non-zero mean here is read
as the cell being too small, and the test is whether the scatter falls when the
cell grows (L=0.50 holds 3-5 channels, L=0.80 holds 10-11).

Scatter is the population standard deviation (ddof=0), the convention used
throughout the paper.

Usage: python3 aggregate_basetensor_seeds.py <dir> [out.csv] [run_id_prefix]
       (prefix selects the campaign, e.g. BTEN for L=0.50 or BT80 for L=0.80)
"""
import glob
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from aggregate_coltensor import read_C, eng_constants     # noqa: E402


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    d = sys.argv[1]
    out_csv = sys.argv[2] if len(sys.argv) > 2 else 'results_basetensor_seeds.csv'

    pat = sys.argv[3] if len(sys.argv) > 3 else '*'
    files = sorted(glob.glob(os.path.join(
        d, 'elasticity_tensor_%s_z95_s*.csv' % pat)))
    if not files:
        print('no elasticity_tensor_%s_z95_s*.csv in %s' % (pat, d))
        return 1

    rows = []
    for fp in files:
        rid = os.path.basename(fp).replace('elasticity_tensor_', '').replace('.csv', '')
        ec = eng_constants(read_C(fp))
        ec['run_id'] = rid
        ec['E_p'] = 0.5 * (ec['E_x'] + ec['E_y'])
        ec['G_ax'] = 0.5 * (ec['G_xz'] + ec['G_yz'])
        ec['inplane_ratio'] = ec['E_y'] / ec['E_x']       # the statistic of interest
        ec['E_ratio'] = ec['E_z'] / ec['E_p']
        ec['G_ratio'] = ec['G_ax'] / ec['G_xy']
        rows.append(ec)

    cols = ['run_id', 'E_x', 'E_y', 'E_z', 'E_p', 'inplane_ratio', 'E_ratio',
            'G_xy', 'G_xz', 'G_yz', 'G_ax', 'G_ratio', 'nu_xy', 'nu_xz', 'nu_yz']
    with open(out_csv, 'w') as f:
        f.write(','.join(cols) + '\n')
        for r in rows:
            f.write(','.join(r['run_id'] if c == 'run_id'
                             else '%.6g' % r[c] for c in cols) + '\n')
    print('wrote %s  (%d packings)' % (out_csv, len(rows)))

    def stat(key):
        v = np.array([r[key] for r in rows])
        return v.mean(), v.std(ddof=0)

    print('\n  per packing:')
    for r in rows:
        print('    %-14s E_x=%.3f E_y=%.3f E_z=%.3f GPa   E_y/E_x=%.4f'
              % (r['run_id'], r['E_x'] / 1e9, r['E_y'] / 1e9, r['E_z'] / 1e9,
                 r['inplane_ratio']))

    print('\n  ensemble (mean +/- population s.d., n=%d):' % len(rows))
    for key, lab, scale in (('E_x', 'E_x (GPa)', 1e9), ('E_y', 'E_y (GPa)', 1e9),
                            ('E_z', 'E_z (GPa)', 1e9),
                            ('inplane_ratio', 'E_y/E_x', 1.0),
                            ('E_ratio', 'E_z/E_xy', 1.0),
                            ('G_ratio', 'G_axial/G_xy', 1.0)):
        m, s = stat(key)
        print('    %-14s %.4f +/- %.4f' % (lab, m / scale, s / scale))

    m, sdev = stat('inplane_ratio')
    sem = sdev / math.sqrt(max(len(rows) - 1, 1))
    print('\n  in-plane isotropy: E_y/E_x = %.4f +/- %.4f (per-cell spread)'
          % (m, sdev))
    print('    mean is %.1f SEM from unity; %d/%d packings above it'
          % (abs(m - 1) / sem if sem else float('inf'),
             sum(1 for r in rows if r['inplane_ratio'] > 1), len(rows)))
    if abs(m - 1) <= 2 * sem:
        print('    -> in-plane isotropy resolved at this cell size')
    else:
        print('    -> not resolved: the per-cell spread is what a cell this size'
              ' can average.\n       Compare the spread against a larger cell'
              ' rather than adding packings;\n       the generator itself is'
              ' unbiased (see check_channel_isotropy.py).')
    return 0


if __name__ == '__main__':
    sys.exit(main())
