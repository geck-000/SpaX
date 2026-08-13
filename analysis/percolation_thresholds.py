"""Which of Pringle's three thresholds should switch the layered description on?

He reports percolation of the pore space in three directions:

  4.6 +- 0.7 %   vertically
  9   +- 2   %   within the layer planes
  14  +- 4   %   across the layer planes

The closure currently switches the bridge factor on at phi_c = 0.05, i.e. the
vertical threshold (Golden's rule of fives). But the vertical threshold is the
condition for brine to DRAIN, not for a lamellar plane to become a plane of
weakness. The mechanical condition is that brine spans a plane, which is the
IN-PLANE threshold, 9%.

Columnar ice carries its lamellar substructure at every depth -- the plate
spacing is set at the growth interface. What changes with warming is not
whether the lamellae exist but whether the brine inside one is continuous. Cold:
the brine in a plane is necked into discrete pockets, the ice across the plane
is continuous, and the material is a particle composite. Warm: the pockets in a
plane merge, the plane becomes a continuous brine sheet, and load must cross it
through ice bridges.

That is exactly in-plane percolation, and it is measured.
"""
import sys
import numpy as np
sys.path.insert(0, r'C:/Users/stirpeg2/SpaX/SpaX/.claude/worktrees/skeletal-eringen-weibull/analysis')
import ez_closure as ez

PHI_VERT, PHI_INPLANE, PHI_ACROSS = 0.046, 0.09, 0.14

print('Pringle thresholds and what each should control')
print('  %.3f  vertical    -> brine can drain (drained vs undrained limit)'
      % PHI_VERT)
print('  %.3f  in-plane    -> a lamellar plane becomes a plane of weakness'
      % PHI_INPLANE)
print('  %.3f  across-plane-> brine finds a path across; bridges severed'
      % PHI_ACROSS)
print('  %.3f  phi_0 (Assur, assumed): plane holds no ice at all' % ez.PHI_0)
print()
print('Note the across-plane threshold sits just below phi_0, which is the')
print('same event described two ways. That is a consistency check on phi_0')
print('that does not depend on Assur: 0.14 +- 0.04 brackets 0.20 at its top.')
print()

T = np.array([-19.1, -17.3, -15.4, -13.6, -11.8, -10.0, -8.2, -6.3, -4.5, -2.7])
S = np.array([7.0, 5.5, 4.8, 4.5, 4.3, 4.3, 4.5, 5.0, 6.0, 8.0])
phi = ez.brine_volume(T, S)
z = (np.arange(10) + 0.5) / 10.0

print('reference column: which slices are layered under each threshold?')
print('  %-6s %8s %10s %10s' % ('z/H', 'phi', 'phi_c=0.05', 'phi_c=0.09'))
for i in range(10):
    print('  %-6.2f %8.4f %10s %10s'
          % (z[i], phi[i],
             'layered' if phi[i] >= 0.05 else 'pockets',
             'layered' if phi[i] >= PHI_INPLANE else 'pockets'))
print('  -> %d of 10 slices at 0.05, %d of 10 at 0.09'
      % (int((phi >= 0.05).sum()), int((phi >= PHI_INPLANE).sum())))


def E_of(p, phi_c, phi_0=ez.PHI_0, n=ez.N_MID):
    Ep = ez.E_ICE * (1 - 1.65 * p)
    b = np.clip(1 - np.sqrt(np.clip(p, 0, phi_0) / phi_0), 0, 1)
    ne = n * (ez.A0_REF_MM / ez.A0_MM) ** ez.SPACING_EXP
    w = np.clip((p - phi_c) / (phi_0 - phi_c), 0, 1)
    return np.maximum(Ep * b ** (ne * w), ez.E_FLOOR)


print('\neffect on the closure at the column brine fractions')
print('  %-8s %12s %12s' % ('phi', 'E(phi_c=.05)', 'E(phi_c=.09)'))
for p in (0.05, 0.07, 0.09, 0.12, 0.15, 0.19):
    print('  %-8.3f %12.3f %12.3f'
          % (p, float(E_of(p, 0.05)), float(E_of(p, PHI_INPLANE))))

for pc, lab in ((0.05, 'phi_c = 0.05'), (PHI_INPLANE, 'phi_c = 0.09')):
    E = E_of(phi, pc)
    h = 0.1
    z0 = float(np.sum(E * h * z) / np.sum(E * h))
    print('\n%s : E_base %.3f  alpha %.4f  z0/H %.4f'
          % (lab, E[-1], E[-1] / E[0], z0))
print('\nKujala measured: E_base 0.86-1.56 GPa, alpha 0.12-0.19, z0/H 0.37-0.39')
