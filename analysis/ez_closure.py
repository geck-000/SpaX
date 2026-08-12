r"""A usable E(z) for sea ice, with every parameter stated and its status given.

This is a CALIBRATED closure, not a prediction. Four of its ingredients are
measured and two are not, and the two that are not carry most of its
uncertainty. Using it is reasonable; quoting it without the band is not.

    E(phi) = E_pocket(phi) * b(phi)^(n_eff * w(phi))
    E_pocket(phi) = 9.37 (1 - 1.65 phi)          [GPa]
    b(phi) = 1 - sqrt(phi / phi_0)               zero at and above phi_0
    n_eff  = n * (a0_ref/a0)^0.69
    w(phi) = clip((phi - phi_c)/(phi_0 - phi_c), 0, 1)

The bridge factor is switched on GRADUALLY between the percolation threshold,
where brine first spans a layer plane, and phi_0, where the plane holds no ice.
Two limits of w bracket this and both were used at some point:

  w = 1 everywhere   applies b in cold ice where no lamellar plane exists,
                     knocking E_top to 6.3 GPa against a measured 7.18-8.60.
                     Wrong.
  w = Heaviside      a step at phi_c. Reproduces Kujala's endpoint ratio well
                     (alpha 0.129 against a measured 0.12-0.19) but leaves
                     E(phi) discontinuous, so moduli between 4.6 and 8.6 GPa
                     have no corresponding brine fraction and the closure
                     cannot be inverted there.

The ramp adopted here is continuous and monotone, so it inverts everywhere, at
the cost of a milder basal knockdown (alpha 0.235). The transition width is not
measured; the ramp and the step bracket it.

INGREDIENT           VALUE            STATUS
E_pocket             above            MEASURED. R^2 = 0.999 over the column
                                      cells, and Mori-Tanaka reproduces it to
                                      under 1% with nothing fitted.
spacing exponent     0.69             MEASURED. Layer-count sweep at fixed
                                      porosity, bridge fraction and count.
a0                   0.35 mm          MEASURED. Pringle et al.: ice lamellae
                                      200-500 um. Midpoint taken.
a0_ref               0.75 mm          the finest spacing our own cells solved
n                    0.53 (0.49-0.59) CALIBRATED, but the three sources agree
                                      closely: our cells report 0.50 for
                                      constriction, Marchenko's profile asks
                                      0.49 and the Gogolaze beam 0.59. At the
                                      measured spacing this becomes an
                                      EFFECTIVE exponent of 0.83-1.00, since
                                      finer lamellae mean more constrictions.
phi_0                0.20 (0.15-0.36) ASSUMED. Assur's constant, consistent
                                      with Pringle's in-plane percolation
                                      threshold but not measured. Second
                                      largest uncertainty.

WHAT IT IS NOT VALID FOR. Above phi_0 the lamellar plane holds no ice and the
closure returns the floor: the material there is skeletal, a dendritic solid in
seawater, and needs its own description. Nothing here is calibrated against
failure or strength data, which sit below the isotropic bounds and so contain
damage this elastic closure does not model.
"""
import numpy as np

E_ICE = 9.37
PHI_0 = 0.20
PHI_C = 0.05            # Golden's rule of fives: where the lamellar plane begins
N_MID, N_LO, N_HI = 0.53, 0.49, 0.59
A0_MM, A0_REF_MM, SPACING_EXP = 0.35, 0.75, 0.69
E_FLOOR = 0.05          # GPa, nominal skeletal residual; see caveat above


def brine_volume(T, S):
    """Frankenstein and Garner, validated against Pringle to +3.2%."""
    T = np.minimum(np.asarray(T, dtype=float), -0.5)
    return np.asarray(S, dtype=float) * (-49.185 / T + 0.532) / 1000.0


def E_of_phi(phi, n=N_MID, phi_0=PHI_0, a0_mm=A0_MM, floor=E_FLOOR):
    """Transverse Young's modulus, GPa, for drained columnar sea ice.

    The spacing enters through the EXPONENT, not as a prefactor. A prefactor
    would penalise the cold end where b tends to one, and at b = 1 the plane is
    continuous ice with no constriction to penalise -- the modulus must return
    E_pocket exactly however fine the lamellae are. Raising the exponent
    instead has that limit built in, and reproduces the measured spacing
    dependence where b is small, which is the regime it was measured in.
    """
    phi = np.asarray(phi, dtype=float)
    E_pocket = E_ICE * (1.0 - 1.65 * phi)

    # The bridge factor applies only where the lamellar plane it describes
    # actually exists. Below the percolation threshold the brine sits in
    # isolated pockets: there is no plane, b has no referent, and evaluating it
    # anyway charges a mechanism that is not present. At the cold surface, where
    # phi is about 0.022, b comes out at 0.67 and would knock the modulus down
    # by 30% for nothing -- taking E_top to 6.3 GPa against a measured
    # 7.18-8.60 (Kujala).
    #
    # The mechanism is switched on by a weight rather than a step. w rises from
    # 0 at phi_c, where the brine has only just begun to span the layer plane,
    # to 1 at phi_0, where the plane holds no ice at all. A step at phi_c is the
    # w -> Heaviside limit of this and was used earlier; it reproduces Kujala's
    # endpoint ratio more closely but leaves E(phi) discontinuous, and hence a
    # band of moduli with no corresponding brine fraction. The weight keeps the
    # closure continuous and invertible at the cost of a milder basal knockdown.
    b = np.clip(1.0 - np.sqrt(np.clip(phi, 0.0, phi_0) / phi_0), 0.0, 1.0)
    n_eff = n * (A0_REF_MM / a0_mm) ** SPACING_EXP
    w = np.clip((phi - PHI_C) / (phi_0 - PHI_C), 0.0, 1.0)
    E = E_pocket * b ** (n_eff * w)
    return np.maximum(E, floor)


def E_column(z, T_surf, T_base, S, **kw):
    """E(z) for a linear thermal profile and a salinity profile.

    z runs 0 at the cold surface to 1 at the base. S may be scalar or an array
    matching z, so a C-shaped or monotonic salinity can be supplied.
    """
    z = np.asarray(z, dtype=float)
    T = T_surf + (T_base - T_surf) * z
    return E_of_phi(brine_volume(T, S), **kw)


def E_band(phi, **kw):
    """Modulus with the exponent band: (low, mid, high)."""
    kw.pop('n', None)
    return (E_of_phi(phi, n=N_HI, **kw), E_of_phi(phi, n=N_MID, **kw),
            E_of_phi(phi, n=N_LO, **kw))


def flexural(E, z):
    z0 = float(np.trapz(E * z, z) / np.trapz(E, z))
    return float(12.0 * np.trapz(E * (z - z0) ** 2, z) / (z[-1] - z[0]) ** 3)


def main():
    print('E(phi) WITH THE BAND')
    print('%8s %10s %10s %10s %10s'
          % ('phi', 'n=%.2f' % N_HI, 'n=%.2f' % N_MID, 'n=%.2f' % N_LO, 'band'))
    for phi in (0.02, 0.05, 0.08, 0.12, 0.16, 0.19, 0.22):
        lo, mid, hi = (float(v) for v in E_band(phi))
        print('%8.3f %10.3f %10.3f %10.3f %9.0f%%'
              % (phi, lo, mid, hi, 100 * (hi / lo - 1)))

    print('\nWORKED EXAMPLE: 1 m first-year column, -20 C surface, S = 6 ppt')
    z = np.linspace(0, 1, 200)
    for tag, n in (('n = %.2f' % N_HI, N_HI), ('n = %.2f' % N_MID, N_MID),
                   ('n = %.2f' % N_LO, N_LO)):
        E = E_column(z, -20.0, -1.8, 6.0, n=n)
        print('  %-10s E_top %5.2f  E_base %5.2f  alpha %.3f  E_flex %5.2f GPa'
              % (tag, E[0], E[-1], E[-1] / E[0], flexural(E, z)))

    print('\nAGAINST THE THREE COMPARISON CASES')
    zz = np.linspace(1e-3, 1.0, 400)
    zc = zz * 32.0
    phi_g = (0.29315 * zc ** 2 - 5.124 * zc + 85.977) / 1000.0
    lo, mid, hi = (flexural(E_of_phi(phi_g, n=n), zz)
                   for n in (N_HI, N_MID, N_LO))
    print('  Gogolaze beam   %.2f-%.2f GPa (mid %.2f) vs measured 0.79-1.42'
          % (lo, hi, mid))

    phi_m = (np.log(7.23 / (4.4 * (1 - 0.62 * zz ** 0.6))) / 4.2) ** 2
    for n, tag in ((N_HI, 'n=%.2f' % N_HI), (N_MID, 'n=%.2f' % N_MID),
                   (0.375, 'n=0.375')):
        E = E_of_phi(phi_m, n=n)
        print('  Marchenko %s  alpha %.3f vs his 0.384' % (tag, E[-1] / E[0]))

    print('\n  Report the band, not the midpoint alone. At the porosities that')
    print('  matter it is a factor of two, and phi_0 would widen it further.')


if __name__ == '__main__':
    main()
