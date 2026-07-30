"""Weibull / weakest-link sensitivity of the stress-localisation measure.

Why
---
The manuscript reports localisation as the 99th percentile of the matrix
stress-concentration field, P99, chosen because the absolute peak of a
discretised field is mesh-limited and drifts with refinement. That is a
numerical argument, not a physical one, and for a quasi-brittle material the
strength is governed by the extreme tail that a percentile truncates. The
revised Section 3.3 therefore states P99 for what it is and gives the
principled generalisation, a Weibull integral over the cell,

    P_f(V) = 1 - exp[ -(1/V0) * INT (sigma_1/sigma_0)^m dV ]

in which the effective concentration factor is an m-weighted norm of the same
SCF field, and P99 corresponds to one particular m.

Two things are needed to make that concrete, and only one of them is new
simulation.

1. The m-sweep itself is free. It needs the per-element (sigma_1, V_e) pairs
   that analysis/scf_extract.py already builds and currently discards after
   reducing them to percentiles. That script now takes an optional --dump
   argument writing those arrays to a .npz, so the expensive ODB pass happens
   once and the whole m-sweep runs offline with no Abaqus licence
   (analysis/weibull_sensitivity.py).

2. The ensemble is new. Section 3.3 also notes that a single periodic cell
   cannot generate a weakest-link size effect: periodicity replicates one cell,
   so the extreme of the tessellated field is the extreme of that cell, and the
   scatter a Weibull analysis requires has to come from the ensemble of
   packings. The published SCF table is one packing per case. This deck
   replicates each of its six microstructures N_SEED times so the extreme-value
   distribution -- and hence a fitted m -- is measurable rather than assumed.

The six cases are reconstructed from the decks that produced the published
table rather than re-derived, so the replicate means stay comparable to the
values already in the manuscript. Only utx is solved: the SCF field is taken
from the uniaxial-tension ODB.
"""
import csv, glob, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_ice_studies import BASE, COLS

HERE   = os.path.dirname(os.path.abspath(__file__))
PARAMS = os.path.join(HERE, os.pardir, 'params')

N_SEED = 5

# published SCF-table label  ->  run_id of the deck row that produced it
CASES = [
    ('CTRL',  'GAS_v00'),        # control, 2% brine
    ('POCK',  'BRINE_iso_G1'),   # brine pockets, 5%
    ('ELON',  'MORF_s60g60'),    # elongated pockets
    ('GAS',   'GAS_v10'),        # gas voids, 10%
    ('CHAN',  'BRINE_chan_G1'),  # brine + channels, 8%
    ('BASE',  'SEAS_w06_z95'),   # warm base, high VoF
]


def find_row(run_id):
    """Locate a deck row by run_id across params/, ignoring decks we generate."""
    for path in sorted(glob.glob(os.path.join(PARAMS, '*.csv'))):
        if os.path.basename(path).startswith('rve_weibull'):
            continue
        try:
            for r in csv.DictReader(open(path)):
                if r.get('run_id') == run_id:
                    return r, os.path.basename(path)
        except Exception:
            continue
    return None, None


def main():
    rows, missing = [], []
    for tag, src_id in CASES:
        src, deck = find_row(src_id)
        if src is None:
            missing.append(src_id)
            continue
        print('  %-5s <- %-16s (%s)' % (tag, src_id, deck))
        for s in range(1, N_SEED + 1):
            r = {c: src.get(c, BASE.get(c, '')) for c in COLS}
            r['run_id'] = 'WBL_%s_s%d' % (tag, s)
            # SCF is read from the uniaxial-tension ODB; no second load case,
            # no bending, no full tensor -- one solve per packing.
            r['Mode'], r['Disp'] = 'Uniaxial Tension X', '0.005'
            r['Mode2'], r['Disp2'] = '', ''
            r['Kappa'] = '0'
            r['full_tensor'] = 'No'
            rows.append(r)
    if missing:
        sys.exit('could not locate deck rows: %s' % ', '.join(missing))

    path = os.path.join(PARAMS, 'rve_weibull.csv')
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        w.writerows(rows)
    print('\nwrote params/rve_weibull.csv  (%d RVEs, %d solves)' % (len(rows), len(rows)))


if __name__ == '__main__':
    main()
