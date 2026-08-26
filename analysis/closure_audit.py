r"""Is the adopted closure physical? Auditing it piece by piece.

    E_layer(phi) = E_pocket(phi) * b(phi)^n,   b = 1 - sqrt(phi),   n = 2

Four ingredients, and they do not all stand up equally.

1. E_pocket. An FE knockdown over the column cells, R^2 = 0.999, reproduced by
   Mori-Tanaka to under 1%. Sound.

2. b = 1 - sqrt(phi). Assur's load-bearing area. Worth checking for internal
   consistency rather than taking on authority: if brine covers a fraction
   (1-b) of a plane and the layers repeat at spacing a0 with thickness t, then
   phi = (1-b) t / a0, and Assur's form gives t = a0 sqrt(phi). So the closure
   says the brine layer is as thick relative to the plate spacing as it is wide
   relative to the plane -- equiaxed pockets growing in two dimensions, which
   is exactly Assur's picture. Self-consistent.

3. Layers at every depth, no morphology switch. Columnar ice carries its
   lamellar substructure throughout; temperature changes the brine in it, not
   whether it exists. Sound, and it removes a knee no measurement shows.

4. The exponent n = 2. THIS IS THE ONE THAT FAILS, and this script is mostly
   about why. Gibson and Ashby's open-cell square law is a LOW-DENSITY
   asymptote, derived for slender ligaments that bend, and quoted as valid for
   relative densities up to roughly 0.3. Assur's b over the column runs 0.52 to
   0.68. At those values the plane is a perforated solid with the ice in the
   majority, not an open network, bending does not dominate, and the scaling is
   nearer linear. Adopting n = 2 there is using a formula outside its stated
   range because it lands close to the measurement.
"""
import numpy as np

import layered_law as law

GA_VALID = 0.3          # Gibson & Ashby: open-cell scaling below about this
E_MEAS_GOGO = 1.421     # root-corrected
MATRIX_FACTOR_RESULT = 3.628


def main():
    print('CONSISTENCY OF ASSUR b  (does it describe a real geometry?)')
    print('%8s %9s %14s %16s' % ('phi', 'b', 't/a0 implied', 'phi recovered'))
    for phi in (0.05, 0.10, 0.15, 0.227):
        b = law.assur_b(phi)
        t_a0 = np.sqrt(phi)            # from phi = (1-b) t/a0 with 1-b = sqrt(phi)
        print('%8.3f %9.3f %14.3f %16.4f'
              % (phi, b, t_a0, (1 - b) * t_a0))
    print('  recovers phi exactly, so the geometry closes on itself.')

    print('\nVALIDITY OF THE b^2 EXPONENT')
    print('  Gibson & Ashby open-cell scaling holds below a relative density')
    print('  of about %.2f. Assur b over this column:' % GA_VALID)
    for phi in (0.05, 0.10, 0.15, 0.227, 0.35):
        b = law.assur_b(phi)
        ok = 'within' if b <= GA_VALID else 'OUTSIDE, by %.1fx' % (b / GA_VALID)
        print('     phi = %.3f -> b = %.3f   %s' % (phi, b, ok))
    print('\n  Every value the column uses is outside the range the square law')
    print('  is derived for, by roughly a factor of two. At b ~ 0.5-0.7 the')
    print('  plane is a perforated solid with ice in the majority; bending does')
    print('  not dominate and the scaling is nearer linear -- which is also what')
    print('  our own cells report, at b^0.85.')

    print('\n  For the square law to be in range, b would have to fall below')
    print('  %.2f, i.e. phi above %.2f -- the skeletal layer, not the column.'
          % (GA_VALID, (1 - GA_VALID) ** 2))

    print('\nWHAT THE HONEST EXPONENT COSTS')
    z = np.linspace(1e-3, 1, 400)
    phi_g = (0.29315 * (z * 32) ** 2 - 5.124 * (z * 32) + 85.977) / 1000.0

    def flexural(E):
        z0 = np.trapz(E * z, z) / np.trapz(E, z)
        return float(12.0 * np.trapz(E * (z - z0) ** 2, z))

    for n, tag in ((2.0, 'b^2 Gibson-Ashby (out of range)'),
                   (1.0, 'b^1 area (in range here)'),
                   (0.85, 'b^0.85 our own cells')):
        ef = flexural(law.layered(phi_g, n))
        print('  %-34s %6.3f GPa   %.2fx measured' % (tag, ef, ef / E_MEAS_GOGO))
    print('  %-34s %6.3f GPa   %.2fx measured'
          % ('the 0.49 matrix factor it replaces', MATRIX_FACTOR_RESULT,
             MATRIX_FACTOR_RESULT / E_MEAS_GOGO))

    print('\nVERDICT')
    print('  Three of the four ingredients are physical and stay. The exponent')
    print('  is not: b^2 is applied a factor of two outside the density range')
    print('  Gibson and Ashby derive it for, and it was adopted because it')
    print('  landed closest to the data. With the exponent its own geometry')
    print('  supports, the closure does WORSE than the fitted factor it was')
    print('  meant to retire.')
    print('\n  That is not a reason to abandon the layered model, whose other')
    print('  parts are sound and which fixed the profile shape. It is a reason')
    print('  not to quote b^2, and to let rve_bracket_nbridges measure what the')
    print('  geometry actually does instead of asserting it.')


if __name__ == '__main__':
    main()
