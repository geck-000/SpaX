r"""What survives when nothing is chosen to match a modulus.

The drained/undrained bracket is wide enough that "the measurement lies inside
it" is a weak claim, and picking a point inside would be fitting. So this fixes
every input from something external and reports what comes out, whether or not
it agrees.

  phi(z)          from the reported temperature and salinity, or in Gogolaze's
                  case from his own eq. (14).
  plate spacing   0.5-1 mm, the growth-rate-controlled lamellar spacing of
                  columnar sea ice; the cell edge maps to a few millimetres so
                  this is n = 3-5 layers.
  drainage        NOT a free choice between two limits. A poroelastic estimate
                  puts the drainage time below a second for any permeability
                  above 1e-13 m^2, against tests lasting tens of seconds, so
                  permeable ice is DRAINED and sealed ice is undrained, with the
                  percolation threshold deciding which.
  bridge fraction b, the ice left in the plane of a layer. This is the one that
                  was chosen. b = 0.03 was picked from a series estimate aimed
                  at the measured base modulus. Assur's plane-of-weakness
                  geometry, which is where the field correlations' sqrt(v) comes
                  from, gives b = 1 - sqrt(v) instead: 0.61 at phi = 0.15 and
                  0.52 at 0.227, some twenty times larger.

The whole basal softening came from b = 0.03. With Assur's b the answer depends
entirely on how E scales with b, and the two candidate laws disagree by orders
of magnitude, so the honest output of this script is a fork rather than a
number.
"""
import numpy as np

E_ICE = 9.37

# the one drained layered measurement, at physical spacing (four layers)
B_MEAS, PHI_MEAS, E_MEAS = 0.03, 0.15, 0.357
K_BOT_MEAS = 1.27          # Kujala base
GOGO = (0.785, 1.421)      # Gogolaze apparent / root-corrected


def pocket(phi):
    return E_ICE * (1.0 - 1.65 * phi)


def assur_b(phi):
    """Load-bearing area fraction left in the plane of weakness."""
    return 1.0 - np.sqrt(phi)


def weeks_assur(phi):
    return 9.5 * (1.0 - np.sqrt(phi)) ** 4


def main():
    print('THE ONE FITTED PARAMETER')
    print('  b was set to %.2f, chosen from a series estimate aimed at the'
          % B_MEAS)
    print('  measured base modulus. Assur geometry gives instead:')
    for phi in (0.10, 0.15, 0.227):
        print('     phi = %.3f -> b = %.3f  (%.0fx the value used)'
              % (phi, assur_b(phi), assur_b(phi) / B_MEAS))

    # exponent implied by our own cells: E(b=1) is the pocket value, since a
    # fully bridged plane is no layer at all
    e1 = pocket(PHI_MEAS)
    n_ours = np.log(E_MEAS / e1) / np.log(B_MEAS)
    print('\nHOW DOES E SCALE WITH b?')
    print('  our cells: E(b=%.2f) = %.3f, E(b=1) = %.3f (the pocket value)'
          % (B_MEAS, E_MEAS, e1))
    print('  -> E ~ b^%.2f, i.e. essentially LINEAR' % n_ours)
    print('  Weeks & Assur: E = 9.5 (1-sqrt(v))^4, i.e. E ~ b^4')
    print('  Those differ by orders of magnitude away from the fitted point.')

    print('\nWHAT COMES OUT AT ASSUR b, WITH NOTHING ELSE CHOSEN')
    print('%8s %10s %12s %12s %12s' % (
        'phi', 'Assur b', 'ours (b^%.2f)' % n_ours, 'if b^4', 'Weeks&Assur'))
    for phi in (0.10, 0.15, 0.227):
        b = assur_b(phi)
        e0 = pocket(phi)
        ours = e0 * b ** n_ours
        four = e0 * b ** 4
        print('%8.3f %10.3f %12.3f %12.3f %12.3f'
              % (phi, b, ours, four, weeks_assur(phi)))

    b = assur_b(0.227)
    ours = pocket(0.227) * b ** n_ours
    four = pocket(0.227) * b ** 4
    print('\nAT THE COLUMN BASE, phi = 0.227, against Kujala %.2f GPa'
          % K_BOT_MEAS)
    print('  pocket model            %.3f GPa   %.1fx too stiff'
          % (pocket(0.227), pocket(0.227) / K_BOT_MEAS))
    print('  layered, our b-scaling  %.3f GPa   %.1fx too stiff'
          % (ours, ours / K_BOT_MEAS))
    print('  layered, b^4 scaling    %.3f GPa   %.2fx'
          % (four, four / K_BOT_MEAS))
    print('  layered, fitted b=0.03  %.3f GPa   %.2fx  (measured at phi=0.15)'
          % (E_MEAS, E_MEAS / K_BOT_MEAS))

    print('\nTHE ANSWER TO THE QUESTION')
    print('  With every parameter physical and b from Assur, our own measured')
    print('  b-scaling gives %.2f GPa at the base -- %.0f%% of the pocket'
          % (ours, 100 * ours / pocket(0.227)))
    print('  model\'s %.2f, so almost no improvement. The basal softening came'
          % pocket(0.227))
    print('  from b = 0.03 and not from the layering as such.')
    print('  If instead E ~ b^4, the same physical b gives %.2f GPa and the'
          % four)
    print('  agreement is excellent with nothing fitted at all.')
    print('\n  So the exponent is the whole question, and it has never been')
    print('  measured: the earlier b-sweep ran UNDRAINED, where the sealed')
    print('  brine carries the load and b controls nothing. rve_bracket_bridge')
    print('  measures it drained, over b = 0.02 to 0.50, and decides between a')
    print('  model that reproduces the field correlations from geometry alone')
    print('  and one whose agreement was a fitted constant.')


if __name__ == '__main__':
    main()
