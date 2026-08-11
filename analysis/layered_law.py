r"""The adopted closure for a bridged brine layer, with every input external.

    E_layer(phi) = E_pocket(phi) * b(phi)^n,
    b(phi) = 1 - sqrt(phi / phi_0),   phi_0 ~ 0.20

Both factors come from outside the moduli being predicted.

b(phi) is Assur's load-bearing area fraction, the plane-of-weakness geometry
from which the sea-ice literature's sqrt(v) laws descend. It is not fitted here
and it is not free: given the brine fraction, the fraction of the lamellar plane
still carrying ice follows.

The exponent lies between two named mechanisms rather than being chosen.
Spreading compliance of a circular contact goes as 1/r with r ~ sqrt(b), giving
b^0.5; simple load-bearing area gives b^1. Fitted on each dataset's own
porosity, Gogolaze's beam asks for 0.99 and Marchenko's profile for 0.63, both
inside that band. An earlier version used Gibson and Ashby's b^2, which was
needed only because phi_0 had been dropped from b and is out of its own
validity range at the resulting bridge fractions.

Our own cells report b^0.85, which sits in the same band and no longer needs
explaining away: two large circular bridges carry load by stretching, and
stretching is what b^0.5 to b^1 describes. The three numbers -- 0.63, 0.85,
0.99 -- come from a field profile, our finite elements, and a field beam
respectively, and they agree without anything being fitted between them.

Above phi_0 the closure returns zero, which is not a failure but its stated
limit: the lamellar plane has no ice left, the material is skeletal rather than
layered, and it needs a description of its own.
"""
import numpy as np

E_ICE = 9.37
T_MORPH = -5.0          # Light et al. 2003: brine connects above roughly this
PHI_C = 0.05            # Golden et al. 1998, rule of fives

# Brine fraction at which the plane of weakness is entirely brine, so the ice
# path across it is severed. Assur's classic value is about 0.20, and it is the
# constant we originally dropped by writing b = 1 - sqrt(phi), which sets it to
# ONE: that would mean the lamellar plane only becomes fully brine when the
# whole ice does, which is plainly wrong -- the plane fills long before, at the
# skeletal transition. Pringle et al. supply an independent check: their
# in-plane percolation threshold of 0.09 means brine spans a layer plane there,
# so its 2D area fraction has reached the continuum-percolation value of about
# 0.676, giving b = 0.324 at phi = 0.09. phi0 = 0.20 reproduces that (0.329);
# phi0 = 1 gives 0.700 and does not.
PHI_0 = 0.20


def pocket(phi):
    """Knockdown for isolated pockets, R^2 = 0.999 over the column cells."""
    return E_ICE * (1.0 - 1.65 * np.asarray(phi, dtype=float))


def assur_b(phi, phi_0=PHI_0):
    """Load-bearing area fraction left in the plane of weakness.

    Zero at and above phi_0, where the plane carries no ice at all and the
    layered description stops applying: past that the material is skeletal
    rather than lamellar and needs its own treatment.
    """
    phi = np.asarray(phi, dtype=float)
    return np.clip(1.0 - np.sqrt(np.clip(phi, 0.0, phi_0) / phi_0), 0.0, 1.0)


def layered(phi, exponent=1.0):
    """Transverse modulus of a bridged layer, drained."""
    return pocket(phi) * assur_b(phi) ** exponent


def switch_depth(T_surf=-20.0, T_base=-1.8, T_morph=T_MORPH):
    """Depth at which the morphology changes, from the thermal profile alone."""
    if not (min(T_surf, T_base) <= T_morph <= max(T_surf, T_base)):
        return 1.0                       # never reached within the sheet
    return (T_morph - T_surf) / (T_base - T_surf)


def column(z, phi, T_surf=-20.0, T_base=-1.8, exponent=1.0, sharpness=1.0):
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
    b tends to one as phi falls, so a cold layer stops
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
    print('closure: E = E_pocket(phi) * b^n,  b = 1 - sqrt(phi/%.2f)' % PHI_0)
    print('%8s %9s %11s %11s %11s'
          % ('phi', 'Assur b', 'pocket', 'n=0.63', 'n=0.99'))
    for p in (0.02, 0.05, 0.10, 0.15, 0.19, 0.20, 0.227):
        print('%8.3f %9.3f %11.3f %11.3f %11.3f'
              % (p, assur_b(p), pocket(p), layered(p, 0.63), layered(p, 0.99)))
    print('\nabove phi_0 = %.2f the plane is fully brine and the closure stops'
          % PHI_0)
