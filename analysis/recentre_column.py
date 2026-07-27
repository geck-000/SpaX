#!/usr/bin/env python3
"""Build the re-centred depth column used by the laminated-plate macro model.

The base slice of results_column.csv (z/H=0.95) comes from a single packing that
the five-packing replicate campaign later showed to be a ~6 sigma low outlier
(4.47 GPa against a 4.85 +/- 0.07 GPa ensemble mean, see results_colseeds.csv).
Every profile in the paper therefore reports the base as the ensemble mean, and
the macro model has to use the same column or it is built on the outlier.

The replicate campaign solved only the x and z load cases, so it carries no E_y
or G_xy for the base layer. The independent full-tensor sweep does: its base
slice returns E_x = 4.853 GPa, within 0.1% of the five-packing mean, so its row
supplies a consistent, fully-solved base layer rather than a scaled guess.
(Scaling the original base row by the ensemble ratio instead changes B/sqrt(AD)
by less than 0.0001, so the choice is immaterial -- it is made this way because
every number then comes from a real six-load-case solve.)

Usage: python3 recentre_column.py results_column.csv results_coltensor.csv out.csv
"""
import csv
import sys


def main():
    if len(sys.argv) != 4:
        print(__doc__)
        return 1
    col_path, tensor_path, out_path = sys.argv[1:4]

    with open(col_path) as f:
        rows = list(csv.DictReader(f))
        fields = rows[0].keys()
    rows.sort(key=lambda r: r['run_id'])

    with open(tensor_path) as f:
        base = [r for r in csv.DictReader(f) if abs(float(r['z']) - 0.95) < 1e-6][0]

    # replace the base layer in place; only the fields the macro model reads
    tgt = rows[-1]
    for key in ('E_x', 'E_y', 'E_z', 'G_xy'):
        tgt[key] = base[key]
    if 'E_eff' in tgt:
        tgt['E_eff'] = base['E_x']

    with open(out_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(fields))
        w.writeheader()
        w.writerows(rows)
    print("wrote %s (base layer %s -> E_x=%.4g Pa from the full-tensor sweep)"
          % (out_path, tgt['run_id'], float(base['E_x'])))
    return 0


if __name__ == '__main__':
    sys.exit(main())
