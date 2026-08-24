"""Report sphere_element_ratio.sh: CalculiX C3D4 and U5+U6 against Abaqus C3D4H.

The sphere campaigns are UNDRAINED ONLY -- Abaqus stored no drained sphere
twin -- so the R = E_und/E_drn trick that carries the layered comparison is
not available here, and E_x has to be compared directly.

That is normally blocked by the packing seed not being recorded: a fresh run
cannot reproduce the stored cell's geometry.  What makes it readable anyway
is that the reference DOES record phi_inclusion per realization, and it
stores several seeds per stem.  So the stored seeds give both a mean and a
SPREAD, and that spread -- Abaqus disagreeing with itself over nothing but
the packing -- is the noise floor any CalculiX-vs-Abaqus difference has to
clear before it means anything about the element.

    python3 calculix/report_sphere_undrained.py <root> <ref.csv> <case>...
"""
import csv
import os
import sys


def val(path, col='E_x'):
    if not os.path.isfile(path):
        return float('nan')
    for r in csv.DictReader(open(path)):
        try:
            return float(r[col])
        except (KeyError, ValueError, TypeError):
            return float('nan')
    return float('nan')


def main(argv):
    root, ref, cases = argv[1], argv[2], argv[3:]
    rows = list(csv.DictReader(open(ref)))

    print('-- undrained E_x: CalculiX against the Abaqus C3D4H seed population --')
    print('%-12s %-8s %14s %14s %9s %9s %11s'
          % ('case', 'element', 'E_x', 'abq mean', 'floor', 'excess', 'verdict'))
    for stem in cases:
        seeds = [r for r in rows if r['run_id'].startswith(stem + '_')
                 and (r.get('E_x') or '').strip()]
        if not seeds:
            print('%-12s %-8s %14s' % (stem, '', 'no reference'))
            continue
        ex = [float(r['E_x']) for r in seeds]
        mean = sum(ex) / len(ex)
        floor = 100.0 * (max(ex) - min(ex)) / mean
        phi = [float(r['phi_inclusion']) for r in seeds if (r.get('phi_inclusion') or '').strip()]

        for arm, tag in (('', 'C3D4'), ('_u6', 'U5+U6')):
            f = os.path.join(root, '%s_und%s.out.csv' % (stem, arm))
            e = val(f)
            if e != e:
                continue
            exc = 100.0 * (e - mean) / mean
            verdict = 'inside' if abs(exc) <= floor else (
                'LOCKING' if exc > 0 else 'too soft')
            print('%-12s %-8s %14.6e %14.6e %8.2f%% %+8.2f%% %11s'
                  % (stem, tag, e, mean, floor, exc, verdict))
        # The geometry that licenses reading the numbers side by side.
        pc = val(os.path.join(root, stem + '_und.out.csv'), 'phi_inclusion')
        print('%-12s %-8s  phi_inclusion  ccx %.4f   abq %.4f..%.4f (n=%d)'
              % ('', '', pc, min(phi), max(phi), len(phi)) if phi else '')
        print('%-12s %-8s  gaps  C3D4 %.1e   U5+U6 %.1e'
              % ('', '',
                 val(os.path.join(root, stem + '_und.out.csv'), 'equilibrium_gap'),
                 val(os.path.join(root, stem + '_und_u6.out.csv'), 'equilibrium_gap')))
        print()

    print('floor  = full range of Abaqus E_x across its stored seeds.')
    print('excess = CalculiX minus the Abaqus mean, as a percent of it.')
    print('An excess inside the floor is indistinguishable from repacking.')


if __name__ == '__main__':
    main(sys.argv)
