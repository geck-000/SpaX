# -*- coding: utf-8 -*-
"""E_x and E_z for the lamellar base test, straight from the ODBs.

    abaqus python lam_extract.py <dir>

Reads the reaction at the driving reference point and the applied displacement,
which is what the production post-processor does, but without needing the deck.
"""
import glob
import os
import re
import sys

from odbAccess import openOdb

RP = {'utx': 'RP-1', 'utz': 'RP-3'}
FACE = {'utx': 1, 'utz': 3}


def modulus(odb_path, mode, L):
    odb = openOdb(path=odb_path, readOnly=True)
    step = odb.steps[odb.steps.keys()[-1]]
    fr = step.frames[-1]
    rpname = None
    for nm in odb.rootAssembly.nodeSets.keys():
        if RP[mode] in nm.upper():
            rpname = odb.rootAssembly.nodeSets[nm]
            break
    dof = FACE[mode] - 1
    R = fr.fieldOutputs['RF'].getSubset(region=rpname).values[0].data[dof]
    u = fr.fieldOutputs['U'].getSubset(region=rpname).values[0].data[dof]
    odb.close()
    if not u:
        return float('nan')
    # macroscopic stress = force / face area ; strain = u / L
    return abs((R / L ** 2) / (u / L))


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else '.'
    L = float(sys.argv[2]) if len(sys.argv) > 2 else 0.50
    # Third argument selects the campaign: LAM_ for the ellipsoid aspect-ratio
    # sweep, SLAB_ for the spanning-layer bridge-fraction sweep.
    pre = sys.argv[3] if len(sys.argv) > 3 else 'LAM_'
    rows = {}
    for f in sorted(glob.glob(os.path.join(d, 'Job-%s*-ut?.odb' % pre))):
        m = re.match(r'Job-%s(\w+)-(ut[xz])\.odb' % pre, os.path.basename(f))
        if not m:
            continue
        tag, mode = m.group(1), m.group(2)
        # A killed solve leaves an .odb that opens but carries no steps, and
        # reading it either raises here or -- worse -- returns a number from a
        # partly written frame. The .sta is the authority on whether the job
        # finished, so gate on it rather than on the .odb existing.
        sta = f[:-4] + '.sta'
        try:
            done = 'COMPLETED SUCCESSFULLY' in open(sta).read()
        except IOError:
            done = False
        if not done:
            print('skipping %s: solve did not complete'
                  % os.path.basename(f))
            continue
        rows.setdefault(tag, {})[mode] = modulus(f, mode, L)

    print('%-10s %10s %10s %10s' % ('case', 'E_x GPa', 'E_z GPa', 'E_z/E_x'))
    order = ['needle', 'plate2', 'plate4', 'plate8']
    for tag in [t for t in order if t in rows] + sorted(
            t for t in rows if t not in order):
        r = rows[tag]
        ex, ez = r.get('utx', float('nan')), r.get('utz', float('nan'))
        print('%-10s %10.3f %10.3f %10.3f' % (tag, ex / 1e9, ez / 1e9, ez / ex))
    print()
    print('measured base modulus (Kujala, mean of four beams): 1.27 GPa')


if __name__ == '__main__':
    main()
