"""Report layered_abaqus_ratio.sh: CalculiX C3D4 against Abaqus C3D4H.

Split out of the shell script so it can be re-run over finished output without
re-solving, and so editing it can never disturb a running run.

    python3 calculix/report_abaqus_ratio.py <root> <ref.csv> "<seeds>" <case>...
"""
import csv
import os
import sys


def val(path, col):
    if not os.path.isfile(path):
        return float('nan')
    for r in csv.DictReader(open(path)):
        try:
            return float(r[col])
        except (KeyError, ValueError, TypeError):
            return float('nan')
    return float('nan')


def abq_val(row, col):
    try:
        return float(row[col])
    except (KeyError, ValueError, TypeError):
        return float('nan')


def main(argv):
    root, ref, seeds, cases = argv[1], argv[2], argv[3].split(), argv[4:]
    abq = {r['run_id']: r for r in csv.DictReader(open(ref))}

    # The direct comparison first: it is only readable because the achieved
    # phase fractions match, which is what says the two packings are equivalent.
    print('-- E_x directly, and the geometry that licenses reading it --')
    print('%-18s %-11s %14s %14s %9s %9s'
          % ('case', 'run', 'E_x', 'E_z', 'porosity', 'phi_soft'))
    for stem in cases:
        for st in ('und', 'drn'):
            f = os.path.join(root, '%s_%s.out.csv' % (stem, st))
            print('%-18s %-11s %14.6e %14.6e %9.5f %9.5f'
                  % (stem, 'ccx ' + st, val(f, 'E_x'), val(f, 'E_z'),
                     val(f, 'porosity'), val(f, 'phi_soft_total')))
            for sd in seeds:
                r = abq.get('%s_%s_%s' % (stem, st, sd))
                if not r or not (r.get('E_x') or '').strip():
                    continue
                print('%-18s %-11s %14.6e %14.6e %9.5f %9.5f'
                      % ('', 'abq %s %s' % (st, sd), abq_val(r, 'E_x'),
                         abq_val(r, 'E_z'), abq_val(r, 'porosity'),
                         abq_val(r, 'phi_soft_total')))
        print()

    print('-- and the ratio, where the packing cancels --')
    print('%-18s %9s %9s %9s %9s %11s'
          % ('case', 'R_ccx', 'R_abq', 'spread', 'excess', 'verdict'))
    for stem in cases:
        eu = val(os.path.join(root, stem + '_und.out.csv'), 'E_x')
        ed = val(os.path.join(root, stem + '_drn.out.csv'), 'E_x')
        rc = eu / ed if ed else float('nan')
        ra = []
        for sd in seeds:
            u, d = abq.get('%s_und_%s' % (stem, sd)), abq.get('%s_drn_%s' % (stem, sd))
            if not u or not d:
                continue
            a, b = abq_val(u, 'E_x'), abq_val(d, 'E_x')
            if b:
                ra.append(a / b)
        if not ra:
            print('%-18s %9.4f %9s' % (stem, rc, 'no ref'))
            continue
        mean = sum(ra) / len(ra)
        # Population spread of the Abaqus R across seeds, as a percentage of
        # its mean: the noise floor the excess has to clear.
        var = sum((x - mean) ** 2 for x in ra) / len(ra)
        sd_pct = 100.0 * var ** 0.5 / mean
        full = 100.0 * (max(ra) - min(ra)) / mean
        exc = 100.0 * (rc - mean) / mean
        verdict = 'LOCKING' if exc > full else ('inside' if abs(exc) <= full
                                                else 'below ref')
        print('%-18s %9.4f %9.4f %8.2f%% %+8.2f%% %11s'
              % (stem, rc, mean, full, exc, verdict))
        print('%-18s %9s %9s %8s   (s.d. %.2f%%, n=%d)   gaps %.1e %.1e'
              % ('', '', '', '', sd_pct, len(ra),
                 val(os.path.join(root, stem + '_und.out.csv'), 'equilibrium_gap'),
                 val(os.path.join(root, stem + '_drn.out.csv'), 'equilibrium_gap')))

    print()
    print('spread = full range of the Abaqus R across seeds, the noise floor.')
    print('excess = how much higher CalculiX C3D4 reads than Abaqus C3D4H.')
    print('Above the floor is locking the hybrid element would have removed.')


if __name__ == '__main__':
    main(sys.argv)
