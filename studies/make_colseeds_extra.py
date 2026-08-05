"""Additional packings per column slice, extending the replicate ensemble.

Section 4.1.1 shows that five packings locate a mean well but estimate its
spread only to about a third of its own value, which is why the base slice is
measured twice and pooled. That left the column non-uniform: ten packings at
the base and five everywhere else. This script builds the extra packings needed
to bring every depth to the same ensemble size.

The microstructure of each slice is reconstructed from results_column.csv, so
the added packings are drawn from exactly the same target statistics as the
original five and the two sets pool without reservation. Only the run_id
differs between replicates; independence comes from the packing seed.

    cd results && python3 ../studies/make_colseeds_extra.py [first] [last]

defaults to seeds 6-10, writing rve_colseeds_extra.csv.

NOTE on seeding. Under hpc/generate_array.sh each task is handed a one-row
slice of the deck, and SpaX_Standalone derives its seed from the row index
within the file it is given. Before that script folded the global row index
into SPAX_SEED, every task packed from the same seed and replicates that
differed only in run_id came out near-identical -- their spread was some sixty
times narrower than a serially generated ensemble. Check the realised scatter
of any replicate campaign before trusting it: for this column it should be of
order 0.3-2% in E_x, not 0.01%.
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_ice_studies import BASE, COLS, write

COLUMN_CSV = 'results_column.csv'
ECHOED = ['L', 'L_mesh', 'Is_Porous', 'E_matrix', 'nu_matrix', 'VoF_sphere',
          'r_avg', 'VoF_void_sphere', 'VoF_incl_sphere', 'E_sphere_inclusion',
          'sphericity_avg', 'PBC_Method', 'Bending_PBC_Type', 'Growth_Direction',
          'generate_channels', 'channel_vof_target']


def main():
    first = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    last = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    rows = []
    for x in csv.DictReader(open(COLUMN_CSV)):
        z = int(round(float(x['run_id'].split('z')[-1])))
        r = dict(BASE)
        r.update({k: x[k] for k in ECHOED})
        r['Growth_Concentration'] = '%.2f' % (0.40 + 0.32 * (z / 100.0))
        for s in range(first, last + 1):
            rr = dict(r)
            rr['run_id'] = 'CSEED_z%02d_s%d' % (z, s)
            rows.append(rr)
    write('rve_colseeds_extra.csv',
          [{c: row.get(c, BASE.get(c, '')) for c in COLS} for row in rows])
    print('wrote rve_colseeds_extra.csv: %d RVEs (seeds %d-%d)'
          % (len(rows), first, last))


if __name__ == '__main__':
    main()
