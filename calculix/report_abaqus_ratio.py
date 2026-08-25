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
    print('%-18s %-12s %9s %9s %9s %9s %11s %10s'
          % ('case', 'element', 'R_ccx', 'R_abq', 'spread', 'excess',
             'verdict', 'max gap'))
    for stem in cases:
        # Abaqus R across seeds sets the noise floor, and is the same
        # reference for both CalculiX element arms.
        ra = []
        for sd in seeds:
            u, d = abq.get('%s_und_%s' % (stem, sd)), abq.get('%s_drn_%s' % (stem, sd))
            if not u or not d:
                continue
            a, b = abq_val(u, 'E_x'), abq_val(d, 'E_x')
            if b:
                ra.append(a / b)
        mean = sum(ra) / len(ra) if ra else float('nan')
        if ra:
            var = sum((x - mean) ** 2 for x in ra) / len(ra)
            sd_pct = 100.0 * var ** 0.5 / mean
            full = 100.0 * (max(ra) - min(ra)) / mean
        else:
            sd_pct = full = float('nan')

        # Two arms on the SAME mesh and the same equations: the stock
        # displacement tet, and the deviatoric tet plus its nodal B-bar
        # patches.  Only the element differs, which is what R isolates.
        # _fbar<c>: F-barES-FEM-T4 at c cyclic smoothings (elements_ccx/
        # fbares.py, patches 0009/0010).  Its E_x comes from the reference-
        # point reaction alone -- an F-bar deck carries no element stress --
        # so its `max gap` column is empty by construction, not by omission.
        for arm, tag in (('', 'C3D4'), ('_u6', 'U5+U6 brine'),
                         ('_u6all', 'U5+U6 both'),
                         ('_fbar0', 'F-barES c=0'), ('_fbar1', 'F-barES c=1'),
                         ('_fbar2', 'F-barES c=2')):
            fu = os.path.join(root, '%s_und%s.out.csv' % (stem, arm))
            # THE DENOMINATOR IS ALWAYS PLAIN C3D4, for every arm.
            #
            # R is only the same quantity in both codes if it is built the
            # same way, and Abaqus's R is E_und(C3D4H) / E_drn(C3D4): it
            # substitutes the hybrid element in the UNDRAINED cell only.  So
            # the element under test is substituted only there too, and the
            # drained twin stays the plain C3D4 both codes share.
            #
            # This is not a formality.  Patching both phases softens the ice
            # by ~4%, and in the drained cell the load runs through thin ice
            # bridges, so E_drn moves 4.3% -- straight into R.  Dividing by
            # that arm's own drained value made the same element read +3.49%
            # at b020 and +0.62% at b040 purely from which denominator was
            # available, which is a reporting artefact, not an element
            # property.  The drained cell is also the CONTROL for mesh
            # equivalence against Abaqus, so it has to be left alone.
            fd = os.path.join(root, stem + '_drn.out.csv')
            eu, ed = val(fu, 'E_x'), val(fd, 'E_x')
            if eu != eu or not ed:
                continue
            rc = eu / ed
            gap = max(val(fu, 'equilibrium_gap'), val(fd, 'equilibrium_gap'))
            if not ra:
                print('%-18s %-12s %9.4f %9s' % (stem, tag, rc, 'no ref'))
                continue
            exc = 100.0 * (rc - mean) / mean
            verdict = ('LOCKING' if exc > full else
                       'inside' if abs(exc) <= full else 'below ref')
            print('%-18s %-12s %9.4f %9.4f %8.2f%% %+8.2f%% %11s %10s'
                  % (stem, tag, rc, mean, full, exc, verdict,
                     '-' if gap != gap else '%.1e' % gap))
        if ra:
            print('%-18s %-12s %9s %9s %8s   (s.d. %.2f%%, n=%d)'
                  % ('', '', '', '', '', sd_pct, len(ra)))

    print()
    print('spread = full range of the Abaqus R across seeds, the noise floor.')
    print('excess = how much higher CalculiX C3D4 reads than Abaqus C3D4H.')
    print('Above the floor is locking the hybrid element would have removed.')
    print('U5+U6 is the nodal-averaged B-bar arm: if it works, its excess')
    print('sits inside the floor where C3D4 sat above it.')


if __name__ == '__main__':
    main(sys.argv)
