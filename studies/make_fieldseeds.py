"""Replicate decks for the three field-comparison columns of Section 4.3.2.

The production column is now the five-packing ensemble mean at every depth
(analysis/build_ensemble_column.py). The columns the field comparison rests on
were not: results_marchenko.csv, results_mono.csv and results_mono2.csv are
single-packing campaigns, so the two halves of that section were on different
conventions. The RMS misfits quoted there (19%, 15% and 8%) and the end
agreement are therefore single-realisation numbers compared against an ensemble
column.

This replicates each of those three campaigns five times per slice, exactly as
make_colseeds.py did for the production column, so the whole of Section 4.3.2
can be put on one convention:

    MSEED  <- results_marchenko.csv   C-shape column, vibrating-beam matrix
    NSEED  <- results_mono.csv        monotonic salinity
    N2SEED <- results_mono2.csv       steeply monotonic salinity

Each campaign is rebuilt from its own results CSV rather than re-derived from
the physics, which guarantees the replicate means are directly comparable with
the single-packing values already published -- the same argument make_colseeds
makes. Every deck parameter is echoed verbatim; only the packing seed differs
between replicates, which is what a distinct run_id delivers under a fixed
SPAX_SEED.

    cd params && python3 ../studies/make_fieldseeds.py   # -> rve_fieldseeds.csv

Then, from a directory holding SpaX_Standalone.py:

    SPAX_SEED=20260723 python3 SpaX_Standalone.py params/rve_fieldseeds.csv out_fieldseeds/

First-order (full_tensor=No): two load cases per RVE, so
3 campaigns x 10 slices x 5 replicates x 2 solves = 300 solves.
"""
import csv
import os

from make_ice_studies import BASE, COLS, write

N_SEED = 5
RES = os.path.join(os.path.dirname(__file__), '..', 'results')

# (results file, run_id prefix for the replicates)
CAMPAIGNS = [
    ('results_marchenko.csv', 'MSEED'),
    ('results_mono.csv',      'NSEED'),
    ('results_mono2.csv',     'N2SEED'),
]

# Deck fields echoed verbatim in the results CSVs (results-col -> deck-col).
# E_matrix is included deliberately: the marchenko and mono2 campaigns were run
# with the vibrating-beam-effective matrix baked into the deck rather than
# applied afterwards, and that has to be preserved or the replicates would not
# reproduce the published levels.
ECHOED = ['L', 'L_mesh', 'Is_Porous', 'E_matrix', 'nu_matrix', 'VoF_sphere',
          'r_avg', 'VoF_void_sphere', 'VoF_incl_sphere', 'E_sphere_inclusion',
          'sphericity_avg', 'PBC_Method', 'Bending_PBC_Type', 'Growth_Direction',
          'generate_channels', 'channel_vof_target']


def depth_of(run_id):
    """MARCH_z95 / MONO_z95 / MON2_z95 -> 95"""
    tail = run_id.rsplit('_z', 1)[-1]
    return int(round(float(''.join(c for c in tail if c.isdigit()))))


def study_fieldseeds():
    rows = []
    for fname, prefix in CAMPAIGNS:
        path = os.path.join(RES, fname)
        if not os.path.isfile(path):
            print('  skipping %s (not found)' % fname)
            continue
        n = 0
        for x in csv.DictReader(open(path)):
            if not x.get('E_eff') or float(x['E_eff']) <= 0:
                continue
            z = depth_of(x['run_id'])
            r = dict(BASE)
            for k in ECHOED:
                if k in x and x[k] != '':
                    r[k] = x[k]
            r['Growth_Concentration'] = '%.2f' % (0.40 + 0.32 * (z / 100.0))
            for s in range(1, N_SEED + 1):
                rr = dict(r)
                rr['run_id'] = '%s_z%02d_s%d' % (prefix, z, s)
                rows.append(rr)
            n += 1
        print('  %-26s %2d slices x %d replicates' % (fname, n, N_SEED))

    write('rve_fieldseeds.csv',
          [{c: row.get(c, BASE.get(c, '')) for c in COLS} for row in rows])
    print('wrote rve_fieldseeds.csv: %d RVEs' % len(rows))


if __name__ == '__main__':
    study_fieldseeds()
