# -*- coding: utf-8 -*-
"""Extract the torsional rigidity from a second-order torsion ODB.

    abaqus python torsion_extract.py <odb> <L> <run_id> <out_csv>

The twist rate alpha is prescribed at RP_K and its reaction is the conjugate
generalised force, so

    K_rve = R(RP_K) / alpha

with units of N.m^3. For a HOMOGENEOUS periodic cell the warping satisfies
grad^2 w = 0 under periodic boundary conditions, hence w = const and

    K_exact = G L^5 / 6

exactly, which is the patch test this extraction is verified against
(analysis/check_torsion_patch.py). Unlike the bending route there is no
cube-versus-plate reference to choose: the closed form is unambiguous.
"""
import os
import sys

from odbAccess import openOdb


def main():
    if len(sys.argv) < 5:
        raise SystemExit(__doc__)
    odb_path, L, run_id, out_csv = sys.argv[1], float(sys.argv[2]), sys.argv[3], sys.argv[4]

    odb = openOdb(path=odb_path, readOnly=True)

    step = None
    for s in odb.steps.keys():
        if 'TORSION' in s.upper():
            step = odb.steps[s]
            break
    if step is None:
        step = odb.steps[odb.steps.keys()[-1]]
    frame = step.frames[-1]

    rpk = None
    for nm in odb.rootAssembly.nodeSets.keys():
        if 'RP_K' in nm.upper():
            rpk = odb.rootAssembly.nodeSets[nm]
            break
    if rpk is None:
        odb.close()
        raise SystemExit('RP_K not found in %s' % odb_path)

    R = frame.fieldOutputs['RF'].getSubset(region=rpk).values[0].data[0]
    alpha = frame.fieldOutputs['U'].getSubset(region=rpk).values[0].data[0]
    K = abs(R / alpha) if alpha else 0.0

    # meshed volumes, so the realised phase fractions travel with the rigidity
    v_solid = 0.0
    try:
        evol = frame.fieldOutputs['EVOL']
        for v in evol.values:
            v_solid += v.data
    except KeyError:
        v_solid = float('nan')

    new = not os.path.exists(out_csv)
    f = open(out_csv, 'a')
    if new:
        f.write('run_id,L,alpha,R_RPK,K_rve,V_solid\n')
    f.write('%s,%.6g,%.6g,%.6g,%.6g,%.6g\n' % (run_id, L, alpha, R, K, v_solid))
    f.close()

    print('%-16s L=%.3f  alpha=%.4g  K_rve=%.6e' % (run_id, L, alpha, K))
    odb.close()


if __name__ == '__main__':
    main()
