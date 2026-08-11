r"""What N is, why it is not a material property, and what that explains.

In the decks N is `n_bridges`: how many ice bridges pierce one lamellar plane
IN THE CELL. That is not a property of the material, because it counts bridges
in a box whose size we chose. The material property is the bridge number
DENSITY, nu = N / L^2, bridges per unit area of lamellar plane -- or
equivalently the mean bridge spacing 1/sqrt(nu).

Working the constriction argument through with that distinction changes what it
predicts. A circular contact of radius r has spreading compliance ~ 1/(E r);
N of them in parallel over a face of area L^2 give L^2/(E r N) per unit area;
a layer of thickness t has compliance t/E_layer per unit area. So

    E_layer ~ E_ice * t * r * N / L^2  =  E_ice * t * r * nu

and with r = sqrt(b / (pi nu)) this is E_layer ~ E_ice * t * sqrt(b nu). Held
at fixed DENSITY that is independent of cell size, as a homogenisation must be.
Held at fixed COUNT it goes as 1/L, because the density is quietly falling as
the box grows.

WHICH IS WHAT WE DID. The cell-size sweep held n_bridges at two while L ran
from 0.25 to 0.625, so the bridge density fell by a factor of 6.25 across it.
The drained modulus fell with it, and I read that as the cell failing to
homogenise. It was not: the microstructure was changing under the test, for the
second time in this campaign and for the same reason as the layer-count sweep
before it -- a parameter that is not intensive being held fixed.
"""
import numpy as np

# cell-size sweep: pitch 0.125 held, n_bridges = 2 held, L varied
L_SWEEP = np.array([0.250, 0.375, 0.500, 0.625])
E_UNCORR = np.array([0.734, 0.503, 0.346, 0.260])
E_CORR = np.array([0.816, 0.625, 0.544, 0.504])
N_FIXED = 2.0

CELL_MM = 3.0            # model cell edge in physical units, from the paper
GOGO_N = 4.3             # bridges per plane implied by the Gogolaze beam


def main():
    print('N in the decks is bridges per PLANE IN THE CELL.')
    print('The material property is nu = N / L^2, bridges per unit area.\n')

    print('CONSTRICTION, done properly:')
    print('  E_layer ~ E_ice * t * r * nu,  r = sqrt(b/(pi nu))')
    print('          ~ E_ice * t * sqrt(b nu)')
    print('  fixed DENSITY -> independent of cell size, as required')
    print('  fixed COUNT   -> goes as 1/L, since density falls as L^2\n')

    print('THE CELL-SIZE SWEEP HELD COUNT, NOT DENSITY')
    print('%8s %8s %10s %12s %12s'
          % ('L', 'nu', 'E uncorr', 'vs 1/L', 'E corr'))
    nu = N_FIXED / L_SWEEP ** 2
    pred = E_UNCORR[0] * L_SWEEP[0] / L_SWEEP
    for L, n, e, p, c in zip(L_SWEEP, nu, E_UNCORR, pred, E_CORR):
        print('%8.3f %8.1f %10.3f %12.3f %12.3f' % (L, n, e, p, c))
    pu = np.polyfit(np.log(L_SWEEP), np.log(E_UNCORR), 1)[0]
    pc = np.polyfit(np.log(L_SWEEP), np.log(E_CORR), 1)[0]
    print('\n  measured slope, uncorrelated bridges : L^%.2f' % pu)
    print('  measured slope, aligned bridges      : L^%.2f' % pc)
    print('  constriction at fixed COUNT predicts : L^-1.00')
    print('\n  The uncorrelated sweep sits within %.0f%% of that, so what looked'
          % (100 * abs(pu + 1.0)))
    print('  like a homogenisation failure was the bridge density falling by')
    print('  %.1fx across the sweep.' % (nu[0] / nu[-1]))
    print('  Aligning bridges shallowed the slope because it removed the extra')
    print('  lateral hop, but it could not fix a microstructure that was still')
    print('  changing.')

    print('\nWHAT N MEANS PHYSICALLY')
    area_mm2 = CELL_MM ** 2
    nu_gogo = GOGO_N / area_mm2
    print('  The cell edge maps to about %.0f mm, so a plane is %.0f mm^2.'
          % (CELL_MM, area_mm2))
    print('  Gogolaze implied N = %.1f, i.e. nu = %.2f bridges per mm^2,'
          % (GOGO_N, nu_gogo))
    print('  a mean bridge spacing of %.2f mm.' % (1.0 / np.sqrt(nu_gogo)))
    print('  Sea-ice plate spacing is 0.5-1 mm, so bridges would sit roughly')
    print('  one to three plate spacings apart within a lamellar plane. That is')
    print('  the right order, and unlike an exponent it is countable in a')
    print('  micro-CT slice.')

    print('\nWHAT IT MEANS FOR US')
    print('  1. The drained cells were never shown NOT to homogenise. The test')
    print('     that said so held bridge count fixed. Redone at fixed density')
    print('     -- n_bridges scaled as L^2 -- it should come out flat.')
    print('  2. Every layered modulus we have quoted carries an implicit nu set')
    print('     by n_bridges = 2 at L = 0.5, i.e. 8 per unit area. Nothing')
    print('     chose that.')
    print('  3. The closure needs both b and nu. Reporting b alone leaves the')
    print('     modulus undetermined by a factor sqrt(nu).')


if __name__ == '__main__':
    main()
