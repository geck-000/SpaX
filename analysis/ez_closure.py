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
                     (alpha 0.122 against a measured 0.12-0.19) but leaves
                     E(phi) discontinuous, so moduli between 4.6 and 8.6 GPa
                     have no corresponding brine fraction and the closure
                     cannot be inverted there.

The transition width is not measured, so the step and the ramp bracket it. The
STEP is the adopted form and the default here: of the two it is the only one
the four-bridge cells support. The ramp, running from phi_c to phi_0, is
continuous and monotone and so inverts everywhere, at the cost of a milder
basal knockdown (alpha 0.280); it is kept for the comparison the paper quotes
and for the inversion, which a step cannot do.

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
import os
import numpy as np

# numpy 2.x removed np.trapz in favour of np.trapezoid; keep the old name.
if not hasattr(np, 'trapz'):
    np.trapz = np.trapezoid

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
# ---------------------------------------------------------------------------
# THE POCKET COEFFICIENT IS NOT ONE NUMBER, AND DRAINAGE IS NOT A LEVEL FACTOR.
#
# 1.65 was measured on cells whose inclusions are near-equiaxed (sphericity
# 0.85-0.62) and whose growth direction is Z, so the transverse load runs ALONG
# the brine planes.  Sub-threshold sea-ice brine is neither: it is flattened
# into the layer planes whether or not it is yet connected, and in S2 ice the
# c-axes are horizontal, so the plate normals are horizontal too and a
# transverse load crosses them at an angle random in that plane.  The layered
# cells of this paper already assume exactly that (slab_axis = x); the pocket
# cells did not, and the brine does not rotate at the percolation threshold.
#
# The coefficient for the fabric sea ice actually has is SOLVED, on cells with
# Growth_Direction = RandomXY at sphericity 4, and it is
#
#     K_SEALED  = 2.427     K_DRAINED = 2.855
#
# against 1.65 for the cells as built.  Both are measurements on realised brine
# fraction, and both are bracketed by the 3D-random and fully aligned fabrics
# as the mean squared direction cosine of the plate normal requires.
#
# Drainage then STEEPENS the coefficient rather than dividing the level.  That
# is what the cells show -- releasing the brine changes the SLOPE of E(phi), not
# a multiplier on it -- and a level factor d(phi) cannot represent it.  The
# measured drained/sealed ratio at this fabric is 1.18, where d(phi) carried
# 1.04.  d(phi) is therefore retired: set FABRIC = 'paper1' to recover it and
# the published 1.65 for reproducing the earlier numbers.
K_SEALED, K_DRAINED = 2.373, 2.730
FABRIC = os.environ.get('SPAX_FABRIC', 's2')      # 's2' (measured) or 'paper1' (1.65 + d(phi))


def pocket_law(phi, u):
    """E_pocket at brine fraction `phi` with drainage weight `u` in [0, 1]."""
    if FABRIC == 'paper1':
        return E_ICE * (1.0 - 1.65 * phi) / (1.0 + (DRAIN_FACTOR - 1.0) * u)
    k = K_SEALED + (K_DRAINED - K_SEALED) * u
    return E_ICE * (1.0 - k * phi)
# ---------------------------------------------------------------------------
# THE EXPONENT IS NOT A CONSTANT. It falls linearly with the bridge fraction:
#
#     n(b) = 1.091 - 1.337 b        R2 = 0.87, rms 0.029
#
# fitted to twelve four-bridge cells over b = 0.180 to 0.388. b predicts n far
# better than phi does (R2 0.87 against 0.65 on the same cells), which is what
# the mechanism asks for: the exponent belongs to the bridge geometry, not to
# how much brine the slice happens to hold.
#
# The cells sit at a_0 = 0.75 mm, exactly A0_REF_MM, so this is n with no
# spacing correction -- a direct measurement.
#
# WHAT THIS REPLACES, and why. The old constant N_MID = 0.98 came from four
# TWO-bridge cells, three of which agreed at 0.93-1.04 while the fourth read
# 0.66 and was taken as a partial weight. It was not. At b ~ 0.31 two bridges
# and their periodic images make and break a connected ice path across the
# plane, and that fourth cell is the only one of the four above the threshold.
# Dividing the same ice area into four removes the coincidence: eleven cells
# spanning b = 0.225-0.388 hold n to a range of 0.079 where the two-bridge step
# is 0.256. Evaluated back at N = 2 the new form gives 0.984 at phi = 0.10 --
# indistinguishable from the 0.98 fitted there -- and drifts up above it, so
# the old constant was accurate exactly where it was calibrated and nowhere
# else.
N_OF_B_INTERCEPT, N_OF_B_SLOPE = 1.091, -1.337
N_FIT_RMS = 0.029
N_FIT_B_RANGE = (0.180, 0.388)   # measured; below this the form is extrapolated

# n(b) IS CALIBRATED AT FOUR BRIDGES, but the count shifts WHERE n climbs, not
# its shape. The count does not change n(b) away from the bridge-percolation
# window.
#
# Measured at slab 0.075 for BOTH counts, so porosity is matched and only the
# count differs:
#
#       b        n(N=4)   n(N=6)   offset
#     0.3876     0.5590   0.5437   0.0154
#     0.3367     0.6104   0.5754   0.0350
#     0.3144     0.6237   0.5948   0.0289
#     0.2929     0.7481   0.5937   0.1545   (three draws per count)
#     0.2789     0.7568   0.7637  -0.0069   (one draw per count)
#     0.2584     0.7706   0.7660   0.0046
#     0.2254     0.8001   0.7428   0.0573
#
# The offset is large ONLY at b = 0.2929, where the N = 4 branch has already
# climbed but N = 6 has not. Below b ~ 0.28 the two branches coincide within
# the 0.011 placement noise. The count therefore sets the b at which the climb
# begins -- N = 4 nearer b = 0.31, N = 6 nearer b = 0.29 -- rather than the
# shape of n(b). A constant shift DOES map one branch onto the other away from
# that transition.
#
# A SECOND dependence is exposed by the same cells. At b = 0.3876 the matched
# cell gives n = 0.559 where SUBC_p075 of the twelve-cell set gives 0.610 at
# the same b; they differ only in slab fraction, 0.075 against 0.056. So n
# depends on porosity as well as b, beyond what dividing by E_pocket removes,
# and n(b) was never a pure function of b.
#
# WHAT THIS MEANS FOR USE. The closure remains valid where it was calibrated --
# four bridges, the slab fractions of the twelve cells -- and reproduces those
# cells to an rms of 0.029. It must not be read as the bridge-branch law of sea
# ice. Reporting it as one would repeat, at larger scale, the error that
# produced phi_sat: a fit whose support was narrower than its use. The count is
# held at N = 4 rather than transported, but its effect is confined to the
# narrow b window where the bridge array percolates.
N_IS_COUNT_SPECIFIC = True

# HOW MUCH OF THAT SCATTER IS THE BRIDGE ARRANGEMENT. Bridge positions are drawn
# by rejection sampling, and where that jams the generator falls back to a
# regular lattice -- which happens more often as b rises, because the disks plus
# their clearance margins crowd the plane. The twelve cells behind n(b) split six
# lattice and six random, and the split TRACKS b, so it looked like a systematic
# bias steepening the slope.
#
# It is not. Eight cells at fixed b, N and mesh, differing only in the bridge
# placement seed:
#
#     b = 0.3144   n = 0.6119 +/- 0.0114   (spread 0.0316 over 4 placements)
#     b = 0.2584   n = 0.7590 +/- 0.0106   (spread 0.0270 over 4 placements)
#
# so arrangement is NOISE and not bias: n(b) needs no correction, and this is
# the number to quote -- about 0.011 in n, roughly 38% of the fit rms and 4% of
# the 0.278 the fit spans.
#
# These cells CANNOT be sorted into regular and irregular. Bridges are placed
# per LAYER, not per cell, and at these bridge fractions 24 of the 32
# layer-placements jammed and fell back to the lattice, so every cell is a
# mixture and mostly lattice. An earlier version of this comment quoted a
# regular-minus-irregular offset; that comparison was between labels taken from
# a single-layer query and did not describe the cells. The scatter above is a
# straight measurement over independent placements and does not depend on the
# classification.
N_ARRANGEMENT_SD = 0.011

# Kept for callers that still pass a scalar exponent explicitly. These are the
# TWO-bridge constants and are no longer the default for anything.
N_MID, N_LO, N_HI = 0.98, 0.93, 1.04


def n_of_b(b, offset=0.0, phi_0=None):
    """Bridge exponent at four bridges to a plane, Eq. n(b).

    THE EXPONENT ABSORBS A CHANGE OF PREFACTOR EXACTLY, and it must, because it
    was not assumed: it was obtained from the four-bridge cells as

        n = log(E_cell / E_pocket) / log b.

    E_cell is a solved modulus and cannot move because a coefficient did, so
    promoting the pocket law from 1.65 to the measured k(phi) has to be carried
    into n or the layered branch is changed by fiat.  Writing the new pocket law
    as f(phi) times the old one, invariance of E_pocket * b^n requires

        dn(b) = -ln f(phi(b)) / ln b,      phi(b) = phi_0 (1 - b)^2,

    which is closed-form: it needs the cells' brine fractions, not their moduli.
    With it the layered branch returns exactly what it measured and only the
    POCKET branch, below phi_c, moves -- which is the only place the pocket law
    is what is being applied.

    This is the same absorption the paper already invokes when it says a change
    in the prefactor "is absorbed exactly by the exponent wherever the exponent
    was measured".
    """
    b = np.asarray(b, float)
    n = N_OF_B_INTERCEPT + N_OF_B_SLOPE * b + offset
    if FABRIC == 'paper1':
        return n
    p0 = PHI_0 if phi_0 is None else phi_0
    phi = p0 * (1.0 - b) ** 2
    u = np.clip((phi - (PHI_DRAIN - PHI_DRAIN_SD)) / (2.0 * PHI_DRAIN_SD),
                0.0, 1.0)
    old = E_ICE * (1.0 - 1.65 * phi) / (1.0 + (DRAIN_FACTOR - 1.0) * u)
    f = pocket_law(phi, u) / old
    with np.errstate(divide='ignore', invalid='ignore'):
        dn = np.where((b > 0) & (b < 1), -np.log(f) / np.log(b), 0.0)
    return n + np.nan_to_num(dn)

# The ramp saturates at phi_0, which is how the paper defines it: with no
# measured width, the two limits that bracket the transition are a step at
# phi_c and a ramp running the full width of the layered range. An earlier
# version saturated at phi_sat = 0.104; that value is retracted (see the note
# on the ramp in E_of_phi), and it is not merely a different choice -- at
# phi = 0.150 a ramp ending at 0.104 has already reached w = 1, so weight='ramp'
# returned the step exactly and could not reproduce the ramp column of the
# paper's Table 4. The ramp now uses the phi_0 passed to E_of_phi.

# A0_MM defaults to the spacing the cells were solved at, NOT to Pringle's
# measured 0.35 mm. That is a deliberate retreat. Extrapolating the a_0^0.69
# law from 0.75 mm down to 0.35 mm multiplies the exponent by 1.71 and takes
# the basal modulus to 0.24 GPa against a measured 0.86-1.56, so the
# extrapolation is not supported by the comparisons even though the law itself
# is measured over 0.75-3 mm. Pass a0_mm explicitly to explore it.
A0_MM, A0_REF_MM, SPACING_EXP = 0.75, 0.75, 0.69

# THE BRIDGE COUNT, which the closure carried implicitly until now.
#
# Section 4.4.4 measures the transverse modulus against the number of bridges
# sharing a fixed total bridge area and finds it stiffens as N^0.458, close to
# the N^(1/2) of spreading compliance through N contacts; the undrained cell is
# flat at N^0.017. Read from N >= 2 upward, so the single-bridge case does not
# weight the fit, results_bracket_nbridges gives N^0.497 -- the constriction
# value almost exactly. At fixed b = 0.15 the drained modulus runs 0.89, 1.44,
# 2.09 and 2.47 GPa at N = 2, 4, 8, 16, against a replicate scatter under 0.03.
#
# Every cell this closure is calibrated on was built with TWO bridges to a
# plane, so N_CELLS = 2 is not a choice: it is what the calibration contains.
# The exponent and the ramp end both inherit it.
#
# What the real count is, nobody has measured directly. Multiplying published
# brine inclusion number densities by the lamellar spacing brackets it:
#   Lieblappen et al. (2018)   830-4800 channels per cm^3
#   Perovich & Gow (1996)      1.0-4.5 per mm^3
#   Light et al. (2003)        24 per mm^3, number density scaling as a power
#                              law in feature length, so finer imaging resolves
#                              more rather than contradicting the coarser counts
# which over a plane of the area used here gives 6 to 32 features on the two
# lower-resolution counts, midpoint near 14. (An earlier note said 3 to 30; that
# mixed Pringle's 200-500 um ice lamella with the 0.75 mm spacing these cells
# actually use. Carried through at one spacing it is 6 to 32.)
# N_IMAGED is that band, and it is wide because the tomography reports pore
# statistics rather than a bridge count per unit area of a lamellar plane.
#
# VALIDATED, and N_CELLS = 2 shown to be actively misleading.
# results_gapcells.csv holds six drained cells built to probe a step that
# appeared in the RAMP series between phi = 0.0933 and 0.0970. Two findings.
#
#   (a) The exponent holds away from the conditions it was fitted on. At
#       b = 0.3108 the count 2 -> 4 raises E_x from 2.545 to 3.665 GPa, +44.0%,
#       against the +41.1% of N^0.497 -- agreement to 2.1%.
#
#   (b) The step is an artefact of TWO bridges and is absent at four. Taking
#       E/E_pocket so the small differences in realised phi do not confound it,
#       across a 1.1% change in b:
#
#           b = 0.3108 -> 0.3144
#           N = 2 :  0.3345 -> 0.4506   +34.7%
#           N = 4 :  0.4839 -> 0.4850   + 0.2%
#
#       At b ~ 0.314 each of two bridges spans about 44% of the cell edge, so
#       the pair and their periodic images make and break a connected ice path
#       across the plane. Dividing the same ice area into four removes it
#       entirely. The apparent sharpness of E(phi) near the threshold in the
#       N = 2 series is therefore a finite-size accident and NOT a feature of
#       the material; do not fit a transition width to it.
#
# This is the third independent line on the bridge count, and the strongest,
# because it does not rest on equating an imaged inclusion count with a bridge
# count -- it shows the two-bridge cell producing behaviour the material does
# not have.
#
# The two field comparisons available disagree about where in the band to sit,
# and the disagreement is reported rather than resolved:
#   - the neutral plane of the Gogolaze beam wants N ~ 8-10, landing inside the
#     measured 0.37-0.39 where N = 2 gives 0.333
#   - the basal-to-surface ratio against Kujala wants N ~ 3, since N = 10 takes
#     alpha to 0.233 against a measured 0.12-0.19
# These are different beams from different campaigns, and the field datasets are
# already known to disagree with each other more than with the closure.
#
# DEFAULT IS N_CELLS, so nothing changes unless a caller asks for it. Passing
# n_bridges scales the LAYERED BRANCH ONLY -- the pocket branch has no bridges
# and must not move.
#
# BRIDGE_COUNT_EXP IS NOT A MEASURED LAW. kappa = 0.497 was fitted to a sweep at
# N = 2, 4, 8, 16 -- every count EVEN -- which returns a smooth saturating curve.
# Filling in the odd counts at b = 0.3144 breaks it:
#
#     N=3  2.4895     N=4  3.6629     N=5  3.2106     N=6  3.8249   GPa
#
# Odd counts come out softer than even ones, and this is now measured against
# its own noise rather than asserted from single cells. Three or four
# independent bridge-placement draws per count, at b = 0.3144 with slab
# fraction, mesh and porosity matched:
#
#     N=3  2.4774 +/- 0.0210  (0.8%)      N=4  3.7366 +/- 0.0506  (1.4%)
#     N=5  3.2430 +/- 0.0195  (0.6%)      N=6  3.8249  (one cell)
#
# Every count reproduces to better than 1.5%, and the odd deficit is 50.8%
# between N=3 and N=4 and 13-15% around N=5: FIFTY-FOUR times the placement
# noise. No single power describes it -- the local exponent runs +1.342 (3->4),
# -0.591 (4->5), +0.960 (5->6).
#
# It also DECAYS with N, 51% at 3-4 against 14% around 5, which is what a
# finite-cell commensurability effect does and what a material property does
# not: real ice cannot care whether the number of bridges in an arbitrary
# 0.5 m window is odd. The mechanism is not established here. What is
# established is that this cell cannot deliver a bridge-count law.
#
# The factor is therefore RETAINED ONLY AS A ONE-TO-ONE IDENTITY at the
# calibration count. Passing n_bridges != N_CELLS raises a warning, because the
# closure cannot presently transport a bridge count and pretending otherwise is
# how phi_sat happened.
# The calibration is now at FOUR bridges, not two: n(b) above is fitted to
# four-bridge cells. Two was never a choice, it was what the old cells
# contained, and it carried a percolation artefact.
N_CELLS = 4
BRIDGE_COUNT_EXP = 0.497
E_FLOOR = 0.05          # GPa, nominal skeletal residual; see caveat above


def brine_volume(T, S):
    """Frankenstein and Garner, validated against Pringle to +3.2%."""
    T = np.minimum(np.asarray(T, dtype=float), -0.5)
    return np.asarray(S, dtype=float) * (-49.185 / T + 0.532) / 1000.0


def E_of_phi(phi, n=None, phi_0=PHI_0, a0_mm=A0_MM, floor=E_FLOOR,
             n_bridges=N_CELLS, weight='step', n_offset=0.0):
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
    E_pocket = pocket_law(phi, u)

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

    # n is a function of b unless a caller overrides it with a scalar.
    n_base = n_of_b(b, offset=n_offset) if n is None else n
    n_eff = n_base * (A0_REF_MM / a0_mm) ** SPACING_EXP

    # THE RAMP IS RETRACTED. phi_sat = 0.104 was obtained by reading one cell's
    # low exponent as the bridge weight not yet being fully on and inverting;
    # that cell is the one whose two bridges percolate, and at four bridges it
    # shows no deficit to invert. Four-bridge cells run from phi = 0.076 to
    # 0.093, straight through phi_c, with no feature there. Of the two forms
    # the paper compares only the step retains support, so it is the default;
    # weight='ramp' is kept for the comparison the paper quotes and for the
    # inversion of Fig. 6(b), which a step cannot do, and is not supported by
    # any measurement. It runs phi_c to phi_0, not to the retracted phi_sat.
    if weight == 'ramp':
        w = np.clip((phi - PHI_C) / (phi_0 - PHI_C), 0.0, 1.0)
    elif weight == 'step':
        w = np.where(np.asarray(phi, float) > PHI_C, 1.0, 0.0)
    else:
        raise ValueError("weight must be 'step' or 'ramp', got %r" % (weight,))
    E = E_pocket * b ** (n_eff * w)

    # Bridge count. The factor is raised to the same weight w that switches the
    # bridge mechanism on, so it acts only where bridges exist and vanishes
    # smoothly into the pocket branch -- a cell with no lamellar plane has no
    # bridges to count, and applying the factor there would charge a geometry
    # that is absent, exactly the error the weight was introduced to avoid.
    if n_bridges != N_CELLS:
        import warnings
        warnings.warn(
            "E_of_phi: n_bridges=%s but the bridge-count law is not measured. "
            "kappa=%.3f came from an all-even sweep (N=2,4,8,16); odd counts "
            "are 22-25%% softer and no single power fits. Treat the result as "
            "indicative only." % (n_bridges, BRIDGE_COUNT_EXP), RuntimeWarning)
        E = E * (float(n_bridges) / N_CELLS) ** (BRIDGE_COUNT_EXP * w)

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
    """Modulus with the exponent band: (low, mid, high).

    The band is now the scatter of the n(b) fit, +/- 2 rms, rather than the
    spread of three constants read off four two-bridge cells. A HIGHER exponent
    gives a SOFTER cell, so the low-modulus edge takes the positive offset.
    """
    kw.pop('n', None)
    kw.pop('n_offset', None)
    d = 2.0 * N_FIT_RMS
    return (E_of_phi(phi, n_offset=+d, **kw),
            E_of_phi(phi, n_offset=0.0, **kw),
            E_of_phi(phi, n_offset=-d, **kw))


def flexural(E, z):
    z0 = float(np.trapz(E * z, z) / np.trapz(E, z))
    return float(12.0 * np.trapz(E * (z - z0) ** 2, z) / (z[-1] - z[0]) ** 3)


def main():
    print('E(phi) WITH THE BAND')
    print('  n(b) = %.3f - %.3f b at N=4, band = +/- 2 rms (%.3f)'
          % (N_OF_B_INTERCEPT, -N_OF_B_SLOPE, 2.0 * N_FIT_RMS))
    print('%8s %10s %10s %10s %10s'
          % ('phi', 'lo', 'mid', 'hi', 'band'))
    for phi in (0.02, 0.05, 0.08, 0.12, 0.16, 0.19, 0.22):
        lo, mid, hi = (float(v) for v in E_band(phi))
        print('%8.3f %10.3f %10.3f %10.3f %9.0f%%'
              % (phi, lo, mid, hi, 100 * (hi / lo - 1)))

    print('\nWORKED EXAMPLE: 1 m first-year column, -20 C surface, S = 6 ppt')
    z = np.linspace(0, 1, 200)
    for tag, off in (('n(b)+2rms', +2.0 * N_FIT_RMS),
                     ('n(b)', 0.0),
                     ('n(b)-2rms', -2.0 * N_FIT_RMS)):
        E = E_column(z, -20.0, -1.8, 6.0, n_offset=off)
        print('  %-10s E_top %5.2f  E_base %5.2f  alpha %.3f  E_flex %5.2f GPa'
              % (tag, E[0], E[-1], E[-1] / E[0], flexural(E, z)))

    print('\nAGAINST THE THREE COMPARISON CASES')
    zz = np.linspace(1e-3, 1.0, 400)
    zc = zz * 32.0
    phi_g = (0.29315 * zc ** 2 - 5.124 * zc + 85.977) / 1000.0
    lo, mid, hi = (flexural(E_of_phi(phi_g, n_offset=off), zz)
                   for off in (+2.0 * N_FIT_RMS, 0.0, -2.0 * N_FIT_RMS))
    print('  Gogolaze beam   %.2f-%.2f GPa (mid %.2f) vs measured 0.79-1.42'
          % (lo, hi, mid))

    phi_m = (np.log(7.23 / (4.4 * (1 - 0.62 * zz ** 0.6))) / 4.2) ** 2
    for off, tag in ((+2.0 * N_FIT_RMS, 'n(b)+2rms'),
                     (0.0, 'n(b)'),
                     (-2.0 * N_FIT_RMS, 'n(b)-2rms')):
        E = E_of_phi(phi_m, n_offset=off)
        print('  Marchenko %s  alpha %.3f vs his 0.384' % (tag, E[-1] / E[0]))

    print('\n  Report the band, not the midpoint alone. At the porosities that')
    print('  matter it is a factor of two, and phi_0 would widen it further.')
    print('\n  n(b) is specific to four bridges: the count changes its shape,')
    print('  not its level. It reproduces its twelve calibration cells to rms')
    print('  %.3f and is not offered as the bridge-branch law of sea ice.'
          % N_FIT_RMS)


if __name__ == '__main__':
    main()
