r"""What phi_0 is, and how well 0.20 is actually supported.

DEFINITION. phi_0 is the brine volume fraction at which the plane of weakness
between two ice platelets becomes entirely brine, so no ice crosses it. Below
it the plane carries load through the ice bridges left between the pockets;
at it the bridges vanish and the plane is severed. It is the point where sea
ice loses tensile strength ACROSS the lamellae, which is why the same constant
appears in Assur's strength relation, sigma = sigma_0 (1 - sqrt(v/v_0)).

It is not the same as the percolation thresholds already in the paper. Sea ice
has several connectivity events at different brine fractions and they answer
different questions:

    phi ~ 0.05   brine connects VERTICALLY, so the ice becomes permeable and
                 drains. Golden's rule of fives. Governs whether a flexural
                 test samples the drained or undrained response.
    phi ~ 0.09   brine connects WITHIN a lamellar plane. Pringle's p_c,parallel.
    phi ~ 0.14   brine connects ACROSS lamellae. Pringle's p_c,perpendicular.
    phi ~ 0.20   the lamellar plane holds no ice at all. phi_0, and the end of
                 the layered description.

WHY 0.20. Two arguments, and they are of different quality.

The first is that it is the value the Assur strength literature uses, around
0.202. That is authority, not derivation, and this script does not lean on it.

The second is a calculation from Pringle's in-plane threshold, and it is worth
setting out because it is what convinced me. If brine spans a lamellar plane at
a bulk fraction of 0.09, then AT that fraction its two-dimensional area
coverage of the plane has reached the continuum-percolation value for discs,
about 0.676. Assur's geometry says that coverage is sqrt(phi/phi_0). Setting
sqrt(0.09/phi_0) = 0.676 gives phi_0 = 0.197.

The weakness in that chain is the 0.676: it is the threshold for randomly
placed OVERLAPPING discs, and real brine pockets are neither randomly placed
nor free to overlap. The sensitivity of phi_0 to that choice is computed below
rather than left implicit.
"""
import numpy as np

PC_PARALLEL = 0.09          # Pringle et al., in-plane brine percolation
PC_PAR_ERR = 0.02
PERC_2D = {
    'overlapping discs (continuum)': 0.676,
    'random close-packed discs': 0.55,
    'square array, touching': np.pi / 4,
    'aligned lanes (lower bound)': 0.50,
}


def phi0_from(pc, coverage):
    """Assur coverage = sqrt(phi/phi_0), so phi_0 = pc / coverage^2."""
    return pc / coverage ** 2


def main():
    print('phi_0 = brine fraction at which the lamellar plane holds no ice.')
    print('Assur strength literature uses about 0.202.\n')

    print('DERIVED from Pringle in-plane percolation p_c = %.2f +- %.2f'
          % (PC_PARALLEL, PC_PAR_ERR))
    print('%-34s %10s %10s' % ('2D threshold assumed', 'coverage', 'phi_0'))
    for name, cov in PERC_2D.items():
        print('%-34s %10.3f %10.3f' % (name, cov, phi0_from(PC_PARALLEL, cov)))

    print('\nsensitivity to Pringle\'s own error bar, at coverage 0.676:')
    for pc in (PC_PARALLEL - PC_PAR_ERR, PC_PARALLEL, PC_PARALLEL + PC_PAR_ERR):
        print('   p_c = %.2f -> phi_0 = %.3f' % (pc, phi0_from(pc, 0.676)))

    lo = phi0_from(PC_PARALLEL - PC_PAR_ERR, 0.676)
    hi = phi0_from(PC_PARALLEL + PC_PAR_ERR, 0.55)
    print('\n  combining both uncertainties: phi_0 = %.2f to %.2f' % (lo, hi))
    print('  so 0.20 is the middle of a range spanning roughly %.0f%%.'
          % (100 * (hi - lo) / 0.2))

    print('\nWHAT phi_0 CHANGES IN THE MODEL')
    print('  b = 1 - sqrt(phi/phi_0), so phi_0 sets how fast the plane empties.')
    print('%8s %10s %10s %10s' % ('phi', 'b(0.15)', 'b(0.20)', 'b(0.25)'))
    for phi in (0.05, 0.10, 0.15, 0.20):
        row = []
        for p0 in (0.15, 0.20, 0.25):
            row.append(max(0.0, 1.0 - np.sqrt(min(phi, p0) / p0)))
        print('%8.3f %10.3f %10.3f %10.3f' % (phi, *row))

    print('\n  The model is sensitive to it: at phi = 0.15 the bridge fraction')
    print('  runs 0.00 to 0.23 over that range of phi_0, and the modulus goes')
    print('  as b to a power near one, so the basal knockdown moves with it.')
    print('  phi_0 therefore needs measuring, not adopting. It is directly')
    print('  observable -- the brine fraction at which micro-CT shows a lamellar')
    print('  plane with no ice crossing it -- and Pringle\'s images span the')
    print('  range in which it should occur.')


if __name__ == '__main__':
    main()
