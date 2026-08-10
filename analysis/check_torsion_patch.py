# -*- coding: utf-8 -*-
"""Patch test for the torsion mode.

For a HOMOGENEOUS periodic cell under a uniform twist alpha about z, the warping
w(x,y) satisfies grad^2 w = 0 with periodic boundary conditions, so w is constant
and the warping vanishes identically. The strains are then gamma_xz = -alpha*y,
gamma_yz = +alpha*x, and the strain energy is

    U = 1/2 G alpha^2 INT (x^2+y^2) dV = G alpha^2 L^5 / 12,

so the generalised force conjugate to alpha is R = dU/dalpha = G alpha L^5 / 6
and the torsional rigidity is

    K = R / alpha = G L^5 / 6      (exact, no free parameter)

Run with abaqus python.
"""
from odbAccess import openOdb
import sys

odb_path = sys.argv[1] if len(sys.argv) > 1 else 'Job-TORH_L240-tor.odb'
L = float(sys.argv[2]) if len(sys.argv) > 2 else 0.24
E = float(sys.argv[3]) if len(sys.argv) > 3 else 9.3922e9
nu = float(sys.argv[4]) if len(sys.argv) > 4 else 0.33
alpha = float(sys.argv[5]) if len(sys.argv) > 5 else 0.11

odb = openOdb(path=odb_path, readOnly=True)
step = odb.steps[odb.steps.keys()[-1]]
frame = step.frames[-1]

rf = frame.fieldOutputs['RF']
u = frame.fieldOutputs['U']

# locate RP_K in the assembly node sets
rpk = None
for nm in odb.rootAssembly.nodeSets.keys():
    if 'RP_K' in nm.upper():
        rpk = odb.rootAssembly.nodeSets[nm]
        break
if rpk is None:
    raise SystemExit('RP_K node set not found; sets = %s'
                     % odb.rootAssembly.nodeSets.keys())

R = rf.getSubset(region=rpk).values[0].data[0]
a_applied = u.getSubset(region=rpk).values[0].data[0]

G = E / (2.0 * (1.0 + nu))
K_num = abs(R / a_applied) if a_applied else 0.0
K_exact = G * L ** 5 / 6.0

print('applied twist alpha        : %.6g   (deck asked for %.6g)' % (a_applied, alpha))
print('reaction at RP_K           : %.6g' % R)
print('K numerical  = R/alpha     : %.6e' % K_num)
print('K exact      = G L^5 / 6   : %.6e' % K_exact)
print('ratio numerical/exact      : %.6f  (%.2f%% off)'
      % (K_num / K_exact, 100 * (K_num / K_exact - 1)))
odb.close()
