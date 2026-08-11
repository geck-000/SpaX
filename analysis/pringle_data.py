r"""What Pringle et al. (2009) actually contains, having now read it.

The paper was consulted to fix nu(phi), the ice-bridge number density in a
lamellar plane. IT DOES NOT REPORT THAT, and two things in it bear against the
way we have been using the literature.

WHAT IS THERE AND IS USABLE

  Ice lamella thickness, paragraph [20]: "the thickness of ice lamellae between
  brine layers is typically in the range 200-500 um". That is the plate
  spacing, and it is a direct measurement at the resolution of the imaging.
  We had been assuming 0.5-1 mm, which is two to five times too coarse.

  Table 1 gives segmented porosity against temperature at S = 9.4 +- 0.5 ppt,
  from 2.24% at -18 C to 11.67% at -3 C. That is a measured phi(T) for a known
  salinity and can be used directly.

  Percolation thresholds, Table 2: vertical 4.6 +- 0.7%, parallel to the layers
  9 +- 2%, perpendicular 14 +- 4%.

WHAT IS NOT THERE

  No ice-bridge count, and no in-plane ice geometry. The 'necks' the paper
  discusses at length are BRINE necks connecting adjacent brine layers, which
  is the opposite phase to our bridges. So nu cannot be taken from this paper.

  Worse for the idea, paragraph [6] records that published inclusion number
  densities disagree by a factor of twenty-four -- Perovich and Gow give 1.0 to
  4.5 per mm^3, Light et al. 24 per mm^3 -- and attributes it to "a power law
  scaling of inclusion number density with length (highlighting the effect of
  imaging resolution)". A quantity that scales with the resolution used to
  measure it is not a material property, so the Light figure should not have
  been carried into an estimate of nu at all.

AND THE PAPER ARGUES AGAINST OUR CLOSURE'S BASIS

  Paragraph [19]: the images show "a pore space much more complicated than
  suggested by simple models of parallel ice lamellae and parallel brine sheets
  or tubes [e.g., Assur, 1960]".

  Paragraph [55]: "The classic model of Assur [1960] ... does not indicate when
  pore segregation occurs and does not allow for transitions in pore
  connectivity. Anisotropy is prescribed through a simple geometric result that
  does not capture the cross-layer connectivity we have addressed."

That is the authors of the measurement we were leaning on saying the geometry
we adopted from Assur is not what their images show.
"""
import numpy as np

# Pringle Table 1: segmented porosity against temperature at S = 9.4 ppt
T_C = np.array([-3.0, -4.0, -5.0, -6.0, -7.0, -8.0, -12.0, -15.0, -18.0])
PHI_PCT = np.array([11.67, 8.81, 7.12, 6.00, 5.22, 4.64, 3.06, 2.54, 2.24])
S_PPT = 9.4

LAMELLA_UM = (200.0, 500.0)     # paragraph [20]
CELL_MM = 3.0                   # our cell edge in physical units


def main():
    print('MEASURED PLATE SPACING, and what it does to our cells')
    lo, hi = LAMELLA_UM[0] / 1000.0, LAMELLA_UM[1] / 1000.0
    print('  Pringle: ice lamellae %.0f-%.0f um thick, i.e. a0 = %.2f-%.2f mm'
          % (*LAMELLA_UM, lo, hi))
    print('  we assumed 0.5-1.0 mm, so the real spacing is %.1f-%.1f times'
          % (0.5 / hi, 1.0 / lo))
    print('  finer than assumed.')
    print('  In a %.0f mm cell that is n_slabs = %.0f to %.0f, against the'
          % (CELL_MM, CELL_MM / hi, CELL_MM / lo))
    print('  1 to 5 layers every campaign has actually run.')

    print('\n  Consequence: our layered cells are too COARSE, not too fine. The')
    print('  layer-count sweep found the modulus falling as spacing tightened,')
    print('  so the physical spacing lies beyond the soft end of what we ran.')

    print('\nMEASURED phi(T) AT KNOWN SALINITY (Table 1, S = %.1f ppt)' % S_PPT)
    print('%8s %10s %14s' % ('T (C)', 'phi', 'Frank-Garner'))
    fg = S_PPT * (-49.185 / T_C + 0.532) / 1000.0
    for t, p, f in zip(T_C, PHI_PCT / 100.0, fg):
        print('%8.0f %10.4f %14.4f' % (t, p, f))
    err = 100 * np.mean((fg - PHI_PCT / 100.0) / (PHI_PCT / 100.0))
    print('  Frankenstein-Garner runs %+.0f%% against his segmented values.'
          % err)
    print('\n  That is NOT a discrepancy, and reading it as one would have been')
    print('  a mistake. His samples are grown from a 50:50 InstantOcean-CsCl')
    print('  mixture, the CsCl added for X-ray contrast, and his Appendix A')
    print('  states the resulting brine volumes are "approximately 25%% less')
    print('  than corresponding brine volumes for pure water ice of the same')
    print('  salinity ... due to the larger mass of Cs+ ions". Correcting for')
    print('  that doping:')
    natural = PHI_PCT / 100.0 / 0.75
    err2 = 100 * np.mean((fg - natural) / natural)
    print('     Frankenstein-Garner vs natural-equivalent porosity: %+.1f%%'
          % err2)
    print('  So the brine relation the whole column rests on agrees with his')
    print('  tomography to a few per cent once his contrast agent is accounted')
    print('  for. That is a real and independent check on it.')

    print('\nWHAT CANNOT BE TAKEN FROM THIS PAPER')
    print('  nu, the ice-bridge density. The paper reports brine necks between')
    print('  layers, which are the other phase. And it records that published')
    print('  inclusion densities differ by 24x purely through imaging')
    print('  resolution, so no literature nu is a material property.')

    print('\n  It also states directly that the pore space is "much more')
    print('  complicated than suggested by simple models of parallel ice')
    print('  lamellae and parallel brine sheets or tubes [Assur, 1960]", and')
    print('  that Assur "does not allow for transitions in pore connectivity".')
    print('  Our closure is built on exactly that geometry.')


if __name__ == '__main__':
    main()
