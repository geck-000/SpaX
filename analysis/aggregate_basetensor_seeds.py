#!/usr/bin/env python3
"""Aggregate the base-slice full-tensor replicates into an ensemble statement.

Reads elasticity_tensor_BTEN_z95_s*.csv (one per packing, six load cases each,
written by SpaX_PostProcess.extract_elasticity_tensor) and reports the mean and
population standard deviation of the engineering constants across packings.

The question this answers: the single base cell of the depth sweep splits
E_x = 4.85 against E_y = 5.02 GPa, a 3.3% in-plane difference. Is that a
realisation effect -- too few channels in one cell to average the in-plane
directions -- or is the base genuinely orthotropic in plane? With five packings
E_y/E_x becomes an ensemble statement: scatter about unity means the former,
a systematic offset means the latter and E_xy would need qualifying.

Scatter is the population standard deviation (ddof=0), the convention used
throughout the paper.

Usage: python3 aggregate_basetensor_seeds.py <dir_with_tensor_csvs> [out.csv]
"""
import glob
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

    files = sorted(glob.glob(os.path.join(d, 'elasticity_tensor_BTEN_z95_s*.csv')))
    if not files:
        print('no elasticity_tensor_BTEN_z95_s*.csv in %s' % d)
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

    m, s = stat('inplane_ratio')
    print('\n  in-plane isotropy: E_y/E_x = %.4f +/- %.4f -> %.1f s.d. from unity'
          % (m, s, abs(m - 1) / s if s else float('inf')))
    print('  %s' % ('CONSISTENT with in-plane isotropy (single-cell split was a '
                    'realisation effect)' if abs(m - 1) <= 2 * s else
                    'SYSTEMATIC in-plane split -- E_xy needs qualifying'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
