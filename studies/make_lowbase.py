"""C-shape column with the basal brine reduced to what the field data imply.

Section 4.4.2 compares our realised brine profile against the one Marchenko's
modulus curve stands for, recovered by inverting the correlation he builds it
from. The two agree in magnitude over the upper nine slices and diverge sharply
at the base: we realise a brine fraction of 0.326 where his curve implies
0.114. That is a statement about the imposed basal microstructure, not about
the homogenisation, and it is testable.

This deck is the test. It is the C-shape recalibrated column (MSEED) with one
change: the two channelled slices carry the basal brine his profile implies
instead of the one Frankenstein-Garner gives at the imposed salinity. Every
other slice, and every other parameter, is echoed verbatim, so any change in
the depth profile is attributable to the basal parameterisation alone.

The targets have to be set on the *realised* fractions rather than the imposed
ones, because the packer overshoots its target systematically -- by a factor of
about two in these two slices (see RATIO below, measured from the campaign this
deck is derived from). Setting the imposed value to the desired realised value
would leave the base roughly twice as brine-rich as intended, which is the
error this deck exists to correct.

    cd params && python3 ../studies/make_lowbase.py    # -> rve_lowbase.csv

Then, from a directory holding SpaX_Standalone.py:

    SPAX_SEED=20260806 python3 SpaX_Standalone.py params/rve_lowbase.csv out_lowbase/

First-order (full_tensor=No): two load cases per RVE, so
10 slices x 5 replicates x 2 solves = 100 solves.
"""
import csv
import os

from make_ice_studies import BASE, COLS, write

N_SEED = 5
RES = os.path.join(os.path.dirname(__file__), '..', 'results')
SRC = 'results_marchenko.csv'
PREFIX = 'LBASE'

# Brine fraction the inverted Marchenko correlation implies, by slice depth,
# and the realised/imposed ratio measured on the campaign being modified.
# Only the two channelled slices are altered; the rest are echoed unchanged.
TARGET_REALISED = {85: 0.0992, 95: 0.1137}
RATIO = {85: 2.08, 95: 2.17}

ECHOED = ['L', 'L_mesh', 'Is_Porous', 'E_matrix', 'nu_matrix', 'VoF_sphere',
          'r_avg', 'VoF_void_sphere', 'VoF_incl_sphere', 'E_sphere_inclusion',
          'sphericity_avg', 'PBC_Method', 'Bending_PBC_Type', 'Growth_Direction',
          'generate_channels', 'channel_vof_target']


def depth_of(run_id):
    tail = run_id.rsplit('_z', 1)[-1]
    return int(round(float(''.join(c for c in tail if c.isdigit()))))


def rescale_base(r, z):
    """Reduce the basal brine, keeping the pocket/channel split as it was."""
    pockets = float(r['VoF_incl_sphere'])
    chans = float(r.get('channel_vof_target', 0) or 0)
    old = pockets + chans
    new = TARGET_REALISED[z] / RATIO[z]
    f = new / old
    r['VoF_incl_sphere'] = '%.4f' % (pockets * f)
    r['channel_vof_target'] = '%.4f' % (chans * f)
    r['VoF_sphere'] = '%.4f' % (float(r['VoF_void_sphere']) + pockets * f)
    print('    z%02d: imposed brine %.3f -> %.4f  (targets realised %.4f)'
          % (z, old, new, TARGET_REALISED[z]))
    return r


def study_lowbase():
    path = os.path.join(RES, SRC)
    if not os.path.isfile(path):
        raise SystemExit('missing %s -- run from params/' % SRC)
    rows, n = [], 0
    for x in csv.DictReader(open(path)):
        if not x.get('E_eff') or float(x['E_eff']) <= 0:
            continue
        z = depth_of(x['run_id'])
        r = dict(BASE)
        for k in ECHOED:
            if k in x and x[k] != '':
                r[k] = x[k]
        r['Growth_Concentration'] = '%.2f' % (0.40 + 0.32 * (z / 100.0))
        if z in TARGET_REALISED:
            r = rescale_base(r, z)
        for s in range(1, N_SEED + 1):
            rr = dict(r)
            rr['run_id'] = '%s_z%02d_s%d' % (PREFIX, z, s)
            rows.append(rr)
        n += 1
    write('rve_lowbase.csv',
          [{c: row.get(c, BASE.get(c, '')) for c in COLS} for row in rows])
    print('wrote rve_lowbase.csv: %d slices x %d replicates = %d RVEs'
          % (n, N_SEED, len(rows)))


if __name__ == '__main__':
    study_lowbase()
