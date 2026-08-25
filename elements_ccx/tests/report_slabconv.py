"""Report slabconv.sh over finished output, without re-solving.

    report_slabconv.py <root> <K/G> <load axis> <n>...

Split out of the shell driver for the same reason report_abaqus_ratio.py was:
a sweep that takes an hour must not have to be repeated to re-read it, and
editing the reporting must not be able to disturb a running solve.
"""
import os
import re
import sys

EPS = 1.0e-3


def faces(dat, ax):
    """(C1111, equilibrium residual) from the driven face reactions."""
    if not os.path.isfile(dat) or os.path.getsize(dat) == 0:
        return float('nan'), float('nan')
    blocks, cur, kind = {}, None, None
    for ln in open(dat):
        m = re.match(r'\s*(forces|displacements).*for set (\S+)', ln)
        if m:
            kind, cur = m.group(1), m.group(2).upper()
            blocks.setdefault((kind, cur), [])
            continue
        if cur and ln.strip() and not ln.lstrip().startswith(('*', 'S T')):
            f = ln.split()
            if len(f) >= 4:
                try:
                    blocks[(kind, cur)].append([float(x) for x in f[1:4]])
                except ValueError:
                    cur = None
    f1 = blocks.get(('forces', 'X1'), [])
    f0 = blocks.get(('forces', 'X0'), [])
    if not f1:
        return float('nan'), float('nan')
    r1 = sum(r[ax] for r in f1)
    bal = abs(r1 + sum(r[ax] for r in f0)) / abs(r1) if r1 else float('nan')
    return r1 / EPS, bal


def main():
    root, kg, load = sys.argv[1], sys.argv[2], sys.argv[3]
    ns = sys.argv[4:]
    ax = {'x': 0, 'y': 1, 'z': 2}[load]
    print('\nK/G = %s   load along %s   R = C1111(und)/C1111(drn), '
          'denominator always plain C3D4' % (kg, load))
    print('%-4s %-8s %13s %13s %13s %8s %8s %9s'
          % ('n', 'el/slab', 'C1111 drn', 'und C3D4', 'und fbar c=1',
             'R C3D4', 'R fbar', 'C3D4 exc'))
    for n in ns:
        b = os.path.join(root, 'n' + n)
        d, bd = faces(os.path.join(b, 'drn', 'm_ccx.dat'), ax)
        u, bu = faces(os.path.join(b, 'und', 'm_ccx.dat'), ax)
        fb, bf = faces(os.path.join(b, 'und_fbar1', 'm_ccx.dat'), ax)
        # How much stiffer plain C3D4 reads than the unlocked element on the
        # SAME mesh: the locking, isolated.  It is the quantity that has to
        # fall as h -> 0 if both elements are consistent.
        exc = 100.0 * (u - fb) / fb if fb == fb and fb else float('nan')
        print('%-4s %-8.0f %13.6e %13.6e %13.6e %8.4f %8.4f %8.2f%%'
              % (n, 0.2 * int(n), d, u, fb,
                 u / d if d else float('nan'),
                 fb / d if d else float('nan'), exc))
        for tag, v in (('drn', bd), ('und', bu), ('fbar', bf)):
            if v == v and v > 1e-6:
                print('     WARNING %s: the two faces carry %.2e of the driven '
                      'reaction out of balance' % (tag, v))


if __name__ == '__main__':
    main()
