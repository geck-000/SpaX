#!/usr/bin/env python3
"""Build the production depth column as the five-packing ensemble mean.

Why this replaces recentre_column.py
------------------------------------
The column used to be the single reference packing at every depth, with only
the warm base swapped for its five-packing mean because that one packing came
out ~6 sigma low. Checking the rest of the column against the same replicate
campaign shows the base was not the only offender:

    z/H     single    ensemble   offset
    0.65     8.530      8.388    +2.4 sd
    0.85     7.348      7.002    +2.4 sd
    0.95     4.469      4.849    -5.8 sd

Two further slices sit at 2.4 sd, and z/H=0.85 is 4.9% high -- larger than the
anisotropy effects the paper resolves. Re-centring one slice and not the others
is not defensible: the same argument applies to all of them, and the replicate
campaign supplies the mean at every depth. So the column is now the ensemble
mean throughout, which is also what the SRVE framing of Section 2.2 demands.

What has to be reconstructed
----------------------------
The replicate campaign (results_colseeds.csv) solved only the x and z load
cases, so it carries E_x, E_z and their Poisson ratios but no E_y or G_xy. The
laminated-plate model needs those, and they come from the independent
full-tensor sweep (results_coltensor.csv), which solves all six load cases at
every depth, rescaled slice by slice onto the ensemble level:

    E_y_ens = E_y_tensor * (E_x_ens / E_x_tensor)

The rescaling is small -- the two campaigns agree on E_x to better than 0.5% at
every depth -- so this transfers a genuine six-load-case in-plane ratio rather
than inventing one. It matters most at the base, where the single reference
packing was the 6 sigma outlier and its in-plane ratio (E_y/E_x = 1.045) is not
the ensemble's; the full-tensor solve gives 1.033 there.

    python3 build_ensemble_column.py results_column.csv results_colseeds.csv \
            results_column_ensemble.csv [results_coltensor.csv]
"""
import csv
import os
import statistics as st
import sys
from collections import defaultdict

SCALED = ('E_y', 'G_xy', 'G_xz', 'G_yz', 'G_eff')   # re-levelled with E_x
AVERAGED = ('E_x', 'E_z', 'E_eff')                  # straight ensemble means
COPIED_MEAN = ('nu_x', 'nu_z', 'nu_eff', 'nu_y')    # ensemble mean where present


def depth_key(rid):
    """CSEED_z95_s3 -> 95 ; ICE_z95 -> 95"""
    core = rid.split('_s')[0]
    return int(''.join(c for c in core.rsplit('_z', 1)[-1] if c.isdigit()))


def main():
    if len(sys.argv) not in (4, 5):
        print(__doc__)
        return 1
    single_path, seeds_path, out_path = sys.argv[1:4]
    tensor_path = sys.argv[4] if len(sys.argv) == 5 else 'results_coltensor.csv'

    tensor = {}
    if os.path.isfile(tensor_path):
        with open(tensor_path) as f:
            for r in csv.DictReader(f):
                tensor[depth_key(r['run_id'])] = r
        print('in-plane E_y and G_xy taken from %s (%d slices), rescaled to the '
              'ensemble level\n' % (tensor_path, len(tensor)))
    else:
        print('WARNING: %s not found; E_y and G_xy will be rescaled from the '
              'single-packing column instead\n' % tensor_path)

    single = {}
    with open(single_path) as f:
        rdr = csv.DictReader(f)
        fields = list(rdr.fieldnames)
        for r in rdr:
            single[depth_key(r['run_id'])] = r

    seeds = defaultdict(list)
    with open(seeds_path) as f:
        for r in csv.DictReader(f):
            seeds[depth_key(r['run_id'])].append(r)

    missing = sorted(set(single) - set(seeds))
    if missing:
        print('no replicates for depth(s) %s -- refusing to build a column that '
              'mixes conventions' % missing)
        return 1

    out = []
    print('%6s %10s %10s %8s   %s'
          % ('z/H', 'single', 'ensemble', 'offset', 'n'))
    for k in sorted(single):
        base = dict(single[k])
        reps = seeds[k]
        ex_s = float(base['E_x'])
        ex_vals = [float(r['E_x']) for r in reps]
        ex_e = st.mean(ex_vals)
        sd = st.pstdev(ex_vals)
        ratio = ex_e / ex_s

        for col in AVERAGED:
            vals = [float(r[col]) for r in reps if r.get(col)]
            if vals:
                base[col] = '%.6e' % st.mean(vals)
        for col in COPIED_MEAN:
            vals = [float(r[col]) for r in reps if r.get(col)]
            if vals:
                base[col] = '%.6f' % st.mean(vals)
        # In-plane quantities: prefer the full-tensor sweep, rescaled onto the
        # ensemble level, over the single packing's own values. At the base the
        # single packing is the 6 sigma outlier, and its in-plane ratio is not
        # the one a representative cell gives.
        src, src_ratio = base, ratio
        if k in tensor:
            src = tensor[k]
            src_ratio = ex_e / float(src['E_x'])
        for col in SCALED:
            if src.get(col):
                base[col] = '%.6e' % (float(src[col]) * src_ratio)

        # keep the derived ratio columns consistent with the new levels
        if base.get('E_z') and base.get('E_x'):
            ez, ex = float(base['E_z']), float(base['E_x'])
            if 'E_anisotropy' in base:
                base['E_anisotropy'] = '%.6f' % (ez / ex)
            if 'E_z_over_xy' in base:
                ey = float(base.get('E_y') or ex)
                base['E_z_over_xy'] = '%.6f' % (ez / (0.5 * (ex + ey)))

        print('%6.2f %10.3f %10.3f %+7.1f sd   %d'
              % (k / 100.0, ex_s / 1e9, ex_e / 1e9,
                 (ex_s - ex_e) / sd if sd else 0.0, len(reps)))
        out.append(base)

    with open(out_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in out:
            w.writerow({k: r.get(k, '') for k in fields})
    print('\nwrote %s (%d slices, ensemble mean at every depth)'
          % (out_path, len(out)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
