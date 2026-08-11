r"""Candidate laws for how a bridged plane carries load, and where the switch is.

Two questions, both raised by the same objection: with Assur's physical bridge
fraction the layered cell only reaches 3.4 GPa at the base, so is the exponent
wrong, and is the sharp drop it produces physical at all?

PART 1 -- the exponent. Four laws are physically motivated, not two:

  b^0.5  constriction. Load spreads into and out of a circular contact, and
         the spreading compliance of a disc of radius r goes as 1/r. With
         N discs of total area fraction b, r ~ sqrt(b), so E ~ sqrt(b). This is
         Holm's contact-resistance result carried over to elasticity.
  b^1    area. The bridges simply carry the load in proportion to the section
         they occupy, with no spreading and no bending.
  b^2    bending-dominated ligaments. Gibson & Ashby's open-cell scaling: if
         the load path through the plane is a set of slender ligaments that
         BEND rather than stretch, E goes as the square of relative density.
  b^4    Weeks & Assur's empirical form for sea ice.

PART 2 -- where layers begin. The switch depth was picked at z/H = 0.75, which
is a second fitted quantity and should not be. Light, Maykut and Grenfell
observed the microstructure change directly: brine sits in isolated pockets in
cold ice and becomes connected above about -5 C. That is a TEMPERATURE
criterion, so for a stated thermal profile the switch depth follows rather than
being chosen.
"""
import numpy as np

E_ICE = 9.37
K_BOT, K_TOP = 1.27, 8.05
T_MORPH = -5.0            # Light et al.: connected brine above about -5 C
T_SURF, T_BASE = -20.0, -1.8


def pocket(phi):
    return E_ICE * (1.0 - 1.65 * phi)


def assur_b(phi):
    return 1.0 - np.sqrt(phi)


def weeks_assur(phi):
    return 9.5 * (1.0 - np.sqrt(phi)) ** 4


LAWS = (('b^0.5  constriction', 0.5),
        ('b^1    area / series', 1.0),
        ('b^2    Gibson-Ashby bending', 2.0),
        ('b^4    Weeks & Assur', 4.0))


def main():
    print('PART 1: what each law gives at Assur b, nothing fitted')
    print('(E(b=1) is the pocket value, since a fully bridged plane is no layer)')
    print()
    print('%-30s %9s %9s %9s' % ('law', 'phi=0.10', 'phi=0.15', 'phi=0.227'))
    for name, n in LAWS:
        vals = [pocket(p) * assur_b(p) ** n for p in (0.10, 0.15, 0.227)]
        print('%-30s %9.3f %9.3f %9.3f' % (name, *vals))
    print('%-30s %9.3f %9.3f %9.3f' % ('Weeks & Assur, their own curve',
                                       *[weeks_assur(p) for p in (0.10, 0.15, 0.227)]))
    print('%-30s %9s %9s %9.2f' % ('Kujala base measurement', '-', '-', K_BOT))

    print('\nAt the base, against the measured %.2f GPa:' % K_BOT)
    for name, n in LAWS:
        v = pocket(0.227) * assur_b(0.227) ** n
        print('  %-30s %6.3f GPa   %5.2fx' % (name, v, v / K_BOT))
    print('\n  The bending scaling lands closest, and it is not an arbitrary')
    print('  choice: a plane held together by slender ice ligaments between')
    print('  brine pockets is exactly the open-cell geometry Gibson and Ashby')
    print('  derive b^2 for, and bending dominates whenever the ligaments are')
    print('  slender. Our cells suggested b^0.85, but that rests on one drained')
    print('  point and an extrapolated endpoint, and the sweep will say.')

    print('\n' + '=' * 64)
    print('PART 2: is the sharp drop physical, and where does it belong?')
    print('=' * 64)
    z = np.linspace(0, 1, 400)
    T = T_SURF + (T_BASE - T_SURF) * z
    z_switch = (T_MORPH - T_SURF) / (T_BASE - T_SURF)
    print('linear thermal profile %.1f C at the surface to %.1f C at the base'
          % (T_SURF, T_BASE))
    print('brine becomes connected above %.1f C (Light et al. 2003)' % T_MORPH)
    print('  -> morphology switches at z/H = %.3f, not at the 0.75 assumed'
          % z_switch)
    print('  -> the layered zone is the bottom %.0f%% of the thickness'
          % (100 * (1 - z_switch)))

    print('\nHow sharp is the drop, and why:')
    phi = np.interp(z, [0, .29, .63, .79, .96, 1.0],
                    [0.104, 0.086, 0.128, 0.168, 0.227, 0.227])
    for name, n in (('b^2 Gibson-Ashby', 2.0), ('b^1 area', 1.0)):
        w = np.clip((z - z_switch) / (1 - z_switch), 0, 1)
        El = pocket(phi) * assur_b(phi) ** n
        E = np.exp((1 - w) * np.log(pocket(phi)) + w * np.log(El))
        i = np.searchsorted(z, z_switch)
        print('  %-18s E(switch) %.2f -> E(base) %.2f GPa over %.0f%% of H'
              % (name, E[i], E[-1], 100 * (1 - z_switch)))

    print('\n  Two things make a steep basal fall physical rather than an')
    print('  artefact of switching a model on. Brine volume itself diverges')
    print('  near the freezing point, since phi ~ S/|T| and |T| -> 1.8 C, so')
    print('  phi rises fastest exactly where the morphology also changes. And')
    print('  the skeletal layer is a real structure with a real top, not a')
    print('  gradual blend.')
    print('\n  What is NOT physical is a step. The transition should be spread')
    print('  over the depth where T runs from %.0f C to %.1f C, which is the'
          % (T_MORPH, T_BASE))
    print('  bottom %.0f%% here, and both the porosity and the morphology vary'
          % (100 * (1 - z_switch)))
    print('  continuously through it.')
    print('\n  It is worth stating plainly that this shape disagrees with the')
    print('  empirical forms. Kerr-Palmer as fitted by Marchenko, E = E0(1 -')
    print('  (1-alpha) d^0.6), falls FASTEST AT THE SURFACE and flattens with')
    print('  depth -- the opposite curvature to a flat column with a soft base.')
    print('  Matching alpha while getting the curvature backwards is not')
    print('  agreement, and that tension is not resolved by any exponent here.')


if __name__ == '__main__':
    main()
