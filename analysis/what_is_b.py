r"""What b is, and the fact that we have been using it to mean two things.

DEFINITION. b is the ice area fraction of a lamellar plane: the share of the
plane between two ice platelets that is still ice rather than brine. Its
complement 1-b is how much of that plane the brine covers. In the generator it
is the total cross-section of the bridge cylinders divided by the cell face,
and it enters the geometry through

    phi_layer = t (1 - b) / a0

with t the layer thickness and a0 the plate spacing. Given phi and a0, choosing
b fixes t: a small b means a THIN layer covering nearly the whole plane, a large
b means a THICKER layer covering only part of it.

So b is not a fudge factor. It is the morphology parameter at fixed brine
volume -- whether the brine in a lamellar plane is a continuous sheet or a set
of discrete pockets -- and micro-CT can measure it directly, since it is just
the in-plane ice fraction of an imaged lamellar plane.

THE PROBLEM. We have used b at two values that describe different ice:

  the CELLS we solved      b = 0.03-0.15, brine covering 85-97% of the plane,
                           i.e. thin nearly-continuous sheets
  the CLOSURE we adopted   b = 0.52-0.78 from Assur, brine covering 22-48%,
                           i.e. discrete pockets in a mostly-ice plane

Those are a factor of seventeen apart in b, and the closure extrapolates
FE behaviour measured at one across to the other on an assumed power law. Both
morphologies are physically real -- Light et al. and Pringle et al. both report
brine sitting in isolated pockets in cold ice and connected sheets in warm ice
-- but they are not the same microstructure and should not be conflated.

WHICH MEANS b SHOULD DEPEND ON TEMPERATURE. That is what the microscopy says,
and it is the b(phi, z) that was suggested earlier and that I argued against.
The objection then was to fitting such a function to the moduli it is meant to
predict, and that objection stands. Taking it from micro-CT does not have that
problem, and is the right way to close the model.
"""
import numpy as np

import layered_law as law


def thickness_ratio(phi, b):
    """t/a0 implied by a given brine fraction and in-plane ice fraction."""
    return phi / np.maximum(1.0 - b, 1e-9)


def main():
    print('b = ice area fraction of the lamellar plane;  1-b = brine coverage')
    print('phi = t (1-b) / a0, so b decides how the same brine is arranged.\n')
    print('%8s %8s %14s %26s' % ('phi', 'b', 't/a0', 'morphology'))
    for phi in (0.10, 0.15, 0.227):
        for b in (0.03, 0.15, law.assur_b(phi)):
            t = thickness_ratio(phi, b)
            if t >= 1.0:
                note = 'IMPOSSIBLE, layer thicker than the spacing'
            elif b < 0.2:
                note = 'thin continuous sheet'
            elif b < 0.45:
                note = 'sheet with wide bridges'
            else:
                note = 'discrete pockets in an ice plane'
            print('%8.3f %8.3f %14.3f %26s' % (phi, b, t, note))
        print()

    print('The cells we solved and the closure we adopted sit at opposite ends')
    print('of that table. Both are real sea-ice morphologies -- cold ice holds')
    print('pockets, warm ice holds sheets -- but a power law fitted at b = 0.03')
    print('and evaluated at b = 0.6 is an extrapolation across a factor of %.0f,'
          % (law.assur_b(0.15) / 0.03))
    print('not an interpolation, and the two describe different ice.')

    print('\nWHAT WOULD SETTLE IT')
    print('  b is directly measurable: image a lamellar plane and take the ice')
    print('  area fraction. Pringle et al. resolved near-parallel intracrystalline')
    print('  brine layers between -18 and -3 C, which is the temperature range')
    print('  the column spans, so the data to fix b(T) already exist in the')
    print('  literature we cite. Until then the exponent measured by')
    print('  rve_bracket_bridge, which spans b = 0.02 to 0.50, is what tells us')
    print('  whether one power law even holds across that range.')


if __name__ == '__main__':
    main()
