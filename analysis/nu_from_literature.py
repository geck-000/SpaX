r"""Can nu be taken from the literature we already cite, and does it agree?

Pringle et al. is the natural source for a bridge density -- his X-ray images
resolve intracrystalline brine layers directly -- but the paper is paywalled
and the abstract reports percolation thresholds and a porosity range, not
counts. Nothing here is taken from those images; that would need the paper.

Light, Maykut and Grenfell DO publish a number density: brine inclusions
averaging about 24 per mm^3 at -15 C, with dimensions from under 0.01 mm to
nearly 10 mm. If those inclusions sit in lamellar planes spaced a0 apart, the
areal density within one plane is n_v * a0, and the ice bridges separating them
are of the same order in number, roughly one per pocket spacing.

That gives an independent estimate of nu, and it does not agree with what our
comparisons imply. Setting the disagreement out is more useful than choosing
one of the two.
"""
import numpy as np

CELL_MODEL = 0.50            # cell edge in model units
CELL_MM = 3.0                # ... and in millimetres, from the paper
MM_PER_UNIT = CELL_MM / CELL_MODEL

LIGHT_NV = 24.0              # brine inclusions per mm^3, Light et al. at -15 C
PLATE_MM = (0.5, 1.0)        # lamellar spacing of columnar sea ice

# nu implied by our work, in model units (bridges per unit area of plane)
NU_CELLS = 8.0               # n_bridges = 2 at L = 0.5
NU_GOGO = 17.0               # from the beam rigidity
NU_MARCH = (3.6, 16.3)       # inverted pointwise down his profile

GOGO_MEAS = 1.421


def to_mm2(nu_model):
    """Model-unit areal density -> per mm^2."""
    return nu_model / MM_PER_UNIT ** 2


def main():
    print('SCALE. The cell edge is %.2f model units and about %.0f mm, so one'
          % (CELL_MODEL, CELL_MM))
    print('model unit is %.0f mm and one model unit of area is %.0f mm^2.\n'
          % (MM_PER_UNIT, MM_PER_UNIT ** 2))

    print('nu IMPLIED BY OUR COMPARISONS')
    print('%-26s %12s %12s' % ('', 'per unit^2', 'per mm^2'))
    print('%-26s %12.1f %12.3f' % ('our cells (n_bridges=2)', NU_CELLS,
                                   to_mm2(NU_CELLS)))
    print('%-26s %12.1f %12.3f' % ('Gogolaze beam', NU_GOGO, to_mm2(NU_GOGO)))
    print('%-26s %12s %12s'
          % ('Marchenko profile', '%.1f-%.1f' % NU_MARCH,
             '%.3f-%.3f' % (to_mm2(NU_MARCH[0]), to_mm2(NU_MARCH[1]))))

    print('\nnu FROM LIGHT ET AL., independent of any of that')
    print('  %.0f inclusions per mm^3; in planes spaced %.1f-%.1f mm apart the'
          % (LIGHT_NV, *PLATE_MM))
    lo, hi = LIGHT_NV * PLATE_MM[0], LIGHT_NV * PLATE_MM[1]
    print('  areal density within one plane is %.0f-%.0f per mm^2.' % (lo, hi))

    ratio_lo = lo / to_mm2(NU_GOGO)
    ratio_hi = hi / to_mm2(NU_MARCH[0])
    print('\n  That is %.0f to %.0f times the densities our comparisons imply.'
          % (ratio_lo, ratio_hi))

    print('\nWHAT THAT DOES TO THE MODULUS')
    print('  E goes as sqrt(b nu), so raising nu to Light\'s value at fixed b')
    print('  stiffens the layered cell by sqrt of that ratio:')
    for nu_real in (lo, hi):
        f = np.sqrt(nu_real / to_mm2(NU_GOGO))
        print('     nu = %4.0f per mm^2 -> Gogolaze beam becomes %.2f GPa (%.1fx'
              % (nu_real, GOGO_MEAS * f, f))
        print('        measured), where the fit needed %.2f.' % GOGO_MEAS)

    print('\n  So a realistic bridge count makes the layered model FAR TOO')
    print('  STIFF, not too soft. Agreement at nu = 0.5 per mm^2 was bought by')
    print('  a bridge density two orders below what the microscopy reports.')

    print('\n  Since only the product b*nu is constrained, the same modulus at')
    print('  Light\'s density needs b smaller by the same factor:')
    for nu_real in (lo, hi):
        b_needed = 0.15 * to_mm2(NU_GOGO) / nu_real
        print('     nu = %4.0f per mm^2 -> b = %.4f, i.e. brine covering %.2f%%'
              % (nu_real, b_needed, 100 * (1 - b_needed)))
        print('        of the lamellar plane, against Assur\'s 0.34 at that phi.')

    print('\nCAVEATS, and they matter here')
    print('  Light counts BRINE INCLUSIONS of every size, from 0.01 mm upward,')
    print('  at one temperature, in first-year Arctic ice. Most of that number')
    print('  is in the smallest inclusions, which need not each interrupt the')
    print('  load path, and ice bridges are not in one-to-one correspondence')
    print('  with brine pockets. So the estimate is an upper bound on the')
    print('  load-bearing bridge density rather than a measurement of it.')
    print('  Even allowing an order of magnitude for that, the gap does not')
    print('  close.')


if __name__ == '__main__':
    main()
