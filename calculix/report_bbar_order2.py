"""Report layered_bbar_order2.sh: can C3D10 + B-bar reach Abaqus's C3D4H?

    python3 calculix/report_bbar_order2.py <root> <ref.csv> <stem> [conv_stem]

Compares three CalculiX configurations on ONE tetrahedralisation against two
Abaqus targets: the hybrid result at the same mesh, and the hybrid result at
the finest mesh the campaign ran, which is the converged answer both codes are
trying to reach.
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


def abq_mean(abq, stem, st, col='E_x'):
    v = []
    for sd in ('s1', 's2', 's3'):
        r = abq.get('%s_%s_%s' % (stem, st, sd))
        if not r or not (r.get(col) or '').strip():
            continue
        try:
            v.append(float(r[col]))
        except ValueError:
            pass
    return sum(v) / len(v) if v else float('nan')


def main(argv):
    root, ref, stem = argv[1], argv[2], argv[3]
    conv = argv[4] if len(argv) > 4 else 'LMESH_m0p0060'
    abq = {r['run_id']: r for r in csv.DictReader(open(ref))}

    def c(tag, col='E_x'):
        return val(os.path.join(root, '%s_%s.out.csv' % (stem, tag)), col)

    a_und, a_drn = abq_mean(abq, stem, 'und'), abq_mean(abq, stem, 'drn')
    k_und, k_drn = abq_mean(abq, conv, 'und'), abq_mean(abq, conv, 'drn')

    print('== %s: three CalculiX elements on one mesh ==' % stem)
    print('%-34s %14s %14s %9s %11s'
          % ('', 'E_x', 'E_z', 'R', 'eq gap'))
    rows = [
        ('CalculiX C3D4        (measured)', None, None, None),
        ('CalculiX C3D10', 'und_o2', 'drn_o2', None),
        ('CalculiX C3D10 + B-bar', 'und_bbar', 'drn_o2', 'gap invalid'),
    ]
    # C3D4 numbers come from the order-1 run of layered_abaqus_ratio.sh
    o1 = os.path.join(os.path.dirname(root.rstrip('/')) or '.', 'out_layerabq')
    e1u = val(os.path.join(o1, '%s_und.out.csv' % stem), 'E_x')
    e1d = val(os.path.join(o1, '%s_drn.out.csv' % stem), 'E_x')
    print('%-34s %14.5e %14.5e %9.4f %11s'
          % ('CalculiX C3D4', e1u, e1d, e1u / e1d if e1d else float('nan'), ''))
    for label, ut, dt, note in rows[1:]:
        eu, ed = c(ut), c(dt)
        print('%-34s %14.5e %14.5e %9.4f %11s'
              % (label, eu, c(ut, 'E_z'), eu / ed if ed else float('nan'),
                 note or ('%.1e' % c(ut, 'equilibrium_gap'))))
    print()
    print('%-34s %14.5e %14.5e %9.4f'
          % ('Abaqus C3D4H, same mesh', a_und, a_drn,
             a_und / a_drn if a_drn else float('nan')))
    print('%-34s %14.5e %14.5e %9.4f'
          % ('Abaqus C3D4H, converged (%s)' % conv, k_und, k_drn,
             k_und / k_drn if k_drn else float('nan')))

    print()
    print('-- distance from each Abaqus target, on undrained E_x and on R --')
    print('%-34s %12s %12s %12s %12s'
          % ('', 'vs same E_x', 'vs same R', 'vs conv E_x', 'vs conv R'))
    r_a = a_und / a_drn if a_drn else float('nan')
    r_k = k_und / k_drn if k_drn else float('nan')
    cand = [('CalculiX C3D4', e1u, e1u / e1d if e1d else float('nan'))]
    for label, ut, dt, _ in rows[1:]:
        eu, ed = c(ut), c(dt)
        cand.append((label.replace('CalculiX ', 'CalculiX '), eu,
                     eu / ed if ed else float('nan')))
    for label, eu, rr in cand:
        print('%-34s %11.2f%% %11.2f%% %11.2f%% %11.2f%%'
              % (label, 100 * (eu - a_und) / a_und, 100 * (rr - r_a) / r_a,
                 100 * (eu - k_und) / k_und, 100 * (rr - r_k) / r_k))
    print()
    print('The same-mesh columns ask whether CalculiX matches Abaqus element')
    print('for element. The converged columns ask the question that matters:')
    print('whether CalculiX can reach the right answer at a mesh it can solve.')
    print('B-bar runs carry no valid equilibrium_gap -- only the stiffness is')
    print('patched, so the reaction cross-check is inconsistent by construction.')


if __name__ == '__main__':
    main(sys.argv)
