r"""A usable E(z) for sea ice, with every parameter stated and its status given.

This is a CALIBRATED closure, not a prediction. Four of its ingredients are
measured and two are not, and the two that are not carry most of its
uncertainty. Using it is reasonable; quoting it without the band is not.

    E(phi) = E_pocket(phi) * b(phi)^(n_eff * w(phi))
    E_pocket(phi) = 9.37 (1 - 1.65 phi)          [GPa]
    b(phi) = 1 - sqrt(phi / phi_0)               zero at and above phi_0
    n_eff  = n * (a0_ref/a0)^0.69
    w(phi) = clip((phi - phi_layer)/(phi_0 - phi_layer), 0, 1)

The branch sits at the IN-PLANE percolation threshold, 0.09, not the vertical
one at 0.046. Columnar ice carries its lamellar substructure at every depth --
the plate spacing is fixed at the growth interface -- so nothing about the
lamellae switches on with warming. What switches on is the brine within a plane
becoming continuous, which is in-plane percolation and is measured. Vertical
percolation governs drainage instead, and the two were previously conflated.

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

# Pringle et al. measure percolation of the pore space in three directions, and
# they control three different things. Using one number for two of them, as an
# earlier version of this file did, conflates drainage with morphology.
#
#   PHI_DRAIN  0.046  vertical percolation. Brine can leave: the drained limit
#                     applies at and above this. Golden's rule of fives.
#   PHI_LAYER  0.09   in-plane percolation. Brine spans a lamellar plane, so
#                     the plane becomes a plane of weakness and load must cross
#                     it through ice bridges. THIS is what switches the layered
#                     description on -- not vertical percolation, which is a
#                     statement about permeability.
#   PHI_CROSS  0.14   across-plane percolation. Brine finds a path across the
#                     planes, so the bridge array no longer blocks. phi_0 must
#                     exceed this, which is the only measurement bearing on
#                     phi_0 at all and bounds it from below.
PHI_DRAIN, PHI_LAYER, PHI_CROSS = 0.046, 0.09, 0.14
# Pringle quotes the vertical threshold as 4.6 +- 0.7 %. The drainage factor is
# ramped over that stated uncertainty rather than stepped at the central value,
# which keeps E(phi) continuous and, unlike the layered ramp of Eq. (5), uses a
# width that is measured rather than chosen.
PHI_DRAIN_SD = 0.007
PHI_0 = 0.20            # Assur; consistent with PHI_CROSS < phi_0 as required

# The ordering matters and settles a question the closure used to assert.
# PHI_LAYER > PHI_DRAIN, so wherever brine spans a layer plane it has already
# percolated vertically and can drain. The layered branch is therefore ALWAYS
# the drained one -- which is worth stating, because Section 4.4.3 measures
# drainage at up to 19x in this morphology, the largest single factor in the
# study.
PHI_C = PHI_LAYER       # retained name: the threshold the closure branches at

# Undrained/drained stiffness ratio for the POCKET morphology, measured in
# Section 4.4 by releasing the brine bulk modulus at fixed geometry. Small
# because an isolated pocket is barely confined; the same release is worth up
# to 19x once the brine spans the cell, which is the whole point of that
# section. Applied above PHI_DRAIN only.
DRAIN_FACTOR = 1.04
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

    # Eq. (2) was calibrated on the column cells, which seal the brine at
    # K = 2.2 GPa, so E_pocket as measured is the UNDRAINED pocket law. Above
    # PHI_DRAIN the brine percolates vertically and the pore pressure relaxes,
    # and Section 4.4 measures that release as worth 1.04x in this morphology.
    # It is the smallest of the four transitions by an order of magnitude, and
    # was previously absorbed rather than represented; it is applied here so
    # that each measured threshold carries the mechanism that belongs to it.
    # Ramped over the threshold's own measurement uncertainty, 4.6 +- 0.7 %, so
    # E(phi) stays continuous. Stepping at the central value left moduli
    # between 8.33 and 8.66 GPa with no preimage, which cost the invertibility
    # the ramped layered branch was adopted to secure.
    u = np.clip((phi - (PHI_DRAIN - PHI_DRAIN_SD)) / (2.0 * PHI_DRAIN_SD),
                0.0, 1.0)
    E_pocket = E_ICE * (1.0 - 1.65 * phi) / (1.0 + (DRAIN_FACTOR - 1.0) * u)

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

    # There is deliberately NO separate strut regime above PHI_CROSS.
    #
    # An earlier version added one, E ~ b^2, reasoning that sparse struts carry
    # load by bending rather than by spreading stress into a contact, so Gibson
    # and Ashby's open-cell law should take over once b drops below about 0.3.
    # The layerskel campaign tested that directly -- layered cells at
    # phi = 0.22, 0.28, 0.35 with b imposed at 0.03, 0.07 and 0.13, every one
    # above phi_0 and well inside the b < 0.3 range -- and measured
    #
    #     drained    E_x ~ b^0.839, b^0.866, b^0.875   (R >= 0.998)
    #     undrained  E_x ~ b^0.02,  b^-0.02, b^0.15    (flat)
    #
    # The exponent sits near n_eff and nowhere near 2, and is stable across a
    # factor of 1.6 in brine fraction. Constriction therefore continues to
    # describe the transverse path even when the plane holds only 3% ice, and
    # the undrained flatness persists too, so confinement holds throughout.
    # The b^2 branch is withdrawn on measurement rather than on argument.
    #
    # What does fail above phi_0 is not the law but Assur's b(phi), which forces
    # b to zero there. These cells carry 0.20 to 0.89 GPa against a floor of
    # 0.05, so the floor is premature: what is missing is a measured b(phi) in
    # that range, not a different E(b).
    return np.maximum(E, floor)


def beyond_bridges(phi):
    """True where the bridge geometry the closure assumes is past its measured
    limit.

    Above PHI_CROSS the brine percolates ACROSS the layer planes, which means it
    has broken through the ice platelets and not merely widened the planes
    between them. The load path the closure models -- platelet, bridge, platelet
    -- is no longer the only thing being cut, so b(phi) understates what is
    happening and the modulus here should be read as an upper bound.

    This is deliberately NOT a third branch. A branch needs a model on the far
    side and the skeletal regime has none here: no cell in this study sits above
    PHI_CROSS with resolved bridges, so there is nothing to switch to. It is
    reported as a validity flag instead.

    It matters because the closure does its steepest work inside the flagged
    band -- 3.44 GPa at PHI_CROSS falling to the floor by phi_0 -- and the basal
    slice of every column in this study lies within it.
    """
    return np.asarray(phi, dtype=float) > PHI_CROSS


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
