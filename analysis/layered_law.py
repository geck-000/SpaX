r"""The adopted closure for a bridged brine layer, with every input external.

    E_layer(phi) = E_pocket(phi) * b(phi)^2,   b(phi) = 1 - sqrt(phi)

Both factors come from outside the moduli being predicted.

b(phi) is Assur's load-bearing area fraction, the plane-of-weakness geometry
from which the sea-ice literature's sqrt(v) laws descend. It is not fitted here
and it is not free: given the brine fraction, the fraction of the lamellar plane
still carrying ice follows.

The exponent 2 is Gibson and Ashby's open-cell scaling. A plane held together by
slender ice ligaments between brine pockets is an open-cell structure, and such
a structure's modulus goes as the square of relative density whenever the
ligaments carry load by BENDING rather than by stretching. Slender ligaments
bend; that is the whole content of the assumption.

WHAT IS ASSERTED RATHER THAN MEASURED, stated plainly because it matters. Our
own cells suggest an exponent near 0.85, not 2. That is not evidence against
Gibson and Ashby: the cells carry TWO LARGE circular bridges, and a fat disc
carries load by stretching, for which the area scaling b^1 is correct. Real
skeletal ice has a dense network of thin ligaments, which is the bending-
dominated geometry b^2 describes and which the cells do not resolve at any mesh
we can afford. Adopting b^2 is therefore a sub-grid closure, and it makes a
falsifiable prediction: holding b fixed and SUBDIVIDING it over more, thinner
bridges should drive the measured exponent up from 1 toward 2. If it does not,
the bending interpretation is wrong and this closure should be withdrawn.

The morphology switch is also external. Light, Maykut and Grenfell observed
brine becoming connected above about -5 C, so a stated thermal profile fixes
the depth at which layers replace pockets; it is not chosen to place a knee.
"""
import numpy as np

E_ICE = 9.37
T_MORPH = -5.0          # Light et al. 2003: brine connects above roughly this
PHI_C = 0.05            # Golden et al. 1998, rule of fives


def pocket(phi):
    """Knockdown for isolated pockets, R^2 = 0.999 over the column cells."""
    return E_ICE * (1.0 - 1.65 * np.asarray(phi, dtype=float))


def assur_b(phi):
    """Load-bearing area fraction left in the plane of weakness."""
    return np.clip(1.0 - np.sqrt(np.clip(phi, 0.0, 1.0)), 1e-6, 1.0)


def layered(phi, exponent=2.0):
    """Transverse modulus of a bridged layer, drained."""
    return pocket(phi) * assur_b(phi) ** exponent


def switch_depth(T_surf=-20.0, T_base=-1.8, T_morph=T_MORPH):
    """Depth at which the morphology changes, from the thermal profile alone."""
    if not (min(T_surf, T_base) <= T_morph <= max(T_surf, T_base)):
        return 1.0                       # never reached within the sheet
    return (T_morph - T_surf) / (T_base - T_surf)


def column(z, phi, T_surf=-20.0, T_base=-1.8, exponent=2.0, sharpness=1.0):
    """E(z) with pockets above a morphology switch and layers below.

    RETAINED FOR COMPARISON ONLY; `layered` applied at every depth is the
    adopted form. Switching the layer geometry on at an isotherm puts a knee in
    E(z) that no measured profile shows: on Marchenko's own porosity this gives
    curvature +0.198 against his -0.073, where the unswitched closure gives
    -0.033, the right sign.

    The switch was also poorly motivated. Columnar ice carries its lamellar
    substructure at every depth, set by the growth process; temperature changes
    the brine volume in the lamellae, not whether they exist. Light et al.'s
    pocket-to-sheet observation concerns connectivity WITHIN a layer and does
    not license switching the layer off. And no switch is needed, since
    b = 1 - sqrt(phi) tends to one as phi falls, so a cold layer stops
    softening of its own accord.
    """
    z = np.asarray(z, dtype=float)
    zs = switch_depth(T_surf, T_base)
    if zs >= 1.0:
        return pocket(phi)
    w = np.clip((z - zs) / max(1.0 - zs, 1e-9), 0.0, 1.0) ** sharpness
    Ep, El = pocket(phi), layered(phi, exponent)
    return np.exp((1.0 - w) * np.log(Ep) + w * np.log(np.maximum(El, 1e-9)))


def drained_here(phi):
    """Permeable ice drains far faster than a flexural test lasts."""
    return np.asarray(phi, dtype=float) > PHI_C


if __name__ == '__main__':
    print('adopted closure: E = E_pocket(phi) * (1 - sqrt(phi))^2')
    print('%8s %9s %11s %11s' % ('phi', 'Assur b', 'pocket', 'layered'))
    for p in (0.05, 0.10, 0.15, 0.20, 0.227, 0.30):
        print('%8.3f %9.3f %11.3f %11.3f'
              % (p, assur_b(p), pocket(p), layered(p)))
    print('\nmorphology switch at z/H = %.3f for a -20 to -1.8 C profile'
          % switch_depth())
