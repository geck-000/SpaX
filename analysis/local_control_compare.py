# -*- coding: utf-8 -*-
"""One increment against ten: read the pairs of ODBs and report the difference.

Runs under `abaqus python` (Python 2.7 with the ODB API), so it is written to
that dialect and imports nothing the Abaqus interpreter does not ship.

What is compared is the reaction force at the driven reference point in the
last frame. That is the quantity the homogenised modulus is built from -- the
first-order extractor takes the volume-averaged stress and the RP displacement
and divides -- and it needs no mesh-dependent averaging to read, so any
difference between the pair is the increment size and nothing else. The two
ODBs of a pair come from the same .inp with one line changed, so the mesh, the
packing and the periodic equations are identical.

    abaqus python analysis/local_control_compare.py <workdir>
"""
import glob
import os
import sys

from odbAccess import openOdb

# The reference point carrying the imposed macro deformation, by load case.
RP = {'utx': 'RP-1', 'uty': 'RP-2', 'utz': 'RP-3'}


def rp_force(path):
    """(frames, frameValue, reaction force at the driven RP in the last frame).

    frameValue is returned and checked by the caller because an ODB being read
    while its job is still running looks entirely healthy: the file opens, the
    last frame has a full field, and the only thing wrong with it is that the
    step has not reached the end. Comparing a complete solve against one is a
    guaranteed disagreement that says nothing about the physics -- and it is the
    same trap as the truncated ODB that csc_solve_array.sh's skip-guard mistakes
    for success, so it is worth catching here rather than discovering it in the
    numbers.
    """
    odb = openOdb(path, readOnly=True)
    try:
        step = odb.steps[odb.steps.keys()[0]]
        n = len(step.frames)
        frame = step.frames[-1]
        t_end = frame.frameValue
        mode = os.path.basename(path).split('-')[-1].split('.')[0]
        name = RP.get(mode, 'RP-1')
        rf = frame.fieldOutputs['RF']
        total = None
        for v in rf.values:
            inst = v.instance
            label = None
            if inst is not None:
                for setname, nset in inst.nodeSets.items():
                    if setname.upper() == name.upper():
                        label = [nd.label for nd in nset.nodes]
                        break
            if label is not None and v.nodeLabel in label:
                total = v.data
                break
        if total is None:
            # Fall back to the largest reaction in the model: on these cells the
            # driven RP carries the whole imposed load, so it is that node.
            best, mag = None, -1.0
            for v in rf.values:
                m = sum(c * c for c in v.data)
                if m > mag:
                    best, mag = v.data, m
            total = best
        return n, t_end, total
    finally:
        odb.close()


def main():
    work = sys.argv[1] if len(sys.argv) > 1 else '.'
    ones = sorted(glob.glob(os.path.join(work, 'Job-LCTL_*.odb')))
    ones = [p for p in ones if '_TEN' not in p]
    if not ones:
        print('no LCTL ODBs in %s' % work)
        return 1

    print('%-26s %6s %6s %16s %16s %10s'
          % ('cell', 'fr(1)', 'fr(10)', 'RF 1 inc', 'RF 10 inc', 'rel diff'))
    worst = 0.0
    incomplete = 0
    compared = 0
    for one in ones:
        ten = one.replace('.odb', '_TEN.odb')
        base = os.path.basename(one)[4:-4]
        if not os.path.exists(ten):
            print('%-26s  -- ten-increment twin missing' % base)
            continue
        if (os.path.exists(one.replace('.odb', '.lck'))
                or os.path.exists(ten.replace('.odb', '.lck'))):
            print('%-26s  -- STILL SOLVING (.lck present); skipped' % base)
            incomplete += 1
            continue
        n1, t1, f1 = rp_force(one)
        n10, t10, f10 = rp_force(ten)
        if abs(t1 - 1.0) > 1e-9 or abs(t10 - 1.0) > 1e-9:
            print('%-26s  -- INCOMPLETE step (end time %.3f vs %.3f); skipped'
                  % (base, t1, t10))
            incomplete += 1
            continue
        # Compare the driven component: the largest by magnitude.
        i = max(range(len(f1)), key=lambda k: abs(f1[k]))
        a, b = float(f1[i]), float(f10[i])
        rel = abs(a - b) / abs(b) if b else float('nan')
        worst = max(worst, rel)
        compared += 1
        print('%-26s %6d %6d %16.8e %16.8e %10.2e'
              % (base, n1, n10, a, b, rel))

    print('')
    if incomplete:
        print('%d pair(s) skipped as incomplete -- rerun once they finish.'
              % incomplete)
    if not compared:
        print('nothing compared yet; no conclusion either way.')
        return 0
    print('compared %d pair(s); largest relative difference: %.3e'
          % (compared, worst))
    if worst < 1e-6:
        print('The two are the same solve. The nine discarded frames were')
        print('discardable, and the campaign can drop them without a caveat.')
    else:
        print('NOT the same solve. Something in these cells is not linear and')
        print('the single-increment results are not interchangeable with the')
        print('LCOL cells they would be pooled with. Do not use the flag.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
