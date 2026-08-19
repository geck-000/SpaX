#!/usr/bin/env python3
r"""Layered cells at the in-plane percolation transition, where the closure is
carried by one cell.

Eq. (5)'s weight ramps the bridge mechanism on from phi_c = 0.09 to
phi_sat = 0.104, and Section 4.5.2 obtains that endpoint by reading the low
exponent of a single cell as a partial weight:

    n(phi) = 0.66, 0.93, 1.04, 0.97  at  phi = 0.099, 0.119, 0.139, 0.167

The upper three average 0.98; the lowest is taken to be 0.66/0.98 = 0.67 of the
way on, and inverting Eq. (6) puts phi_sat just above 0.10. That is one sample
inside a window 0.014 wide, and nothing at all below phi_c -- the condition
that would have constrained the low end, LCOL p060, is the one that would not
mesh. Seed scatter on the surviving cell is small (n = 0.641, 0.656, 0.686, so
phi_sat = 0.103-0.104), so what is unmeasured is not the value but the SHAPE:
one point cannot distinguish the linear ramp from the step, and those two are
quoted in Section 4.5.1 as bracketing a transition width that is not measured.

The window matters out of proportion to its width. Across it E(phi) moves by a
factor of three -- at phi = 0.099 the closure returns 2.29 GPa with the ramp
ending at 0.095 and 6.47 GPa with it ending at 0.16 -- and the invertibility
adopted in Section 4.5.1 makes every measured modulus in that band a statement
about phi. The reference column of Table 1 steps from 0.068 to 0.150 and skips
the window entirely, which is why nothing downstream of it has noticed; the
Gogolaze beam does not. Seven of its twelve slices fall between 0.086 and
0.125, and its predicted apparent modulus runs 2.45 GPa at phi_sat = 0.095 to
3.61 GPa at 0.20, with the neutral axis moving 0.328H to 0.378H.

Three decks, in the order their conclusions depend on each other.

RAMPB -- rve_layerb.csv, already written by make_layer_decks.py and never
    solved. b swept independently of phi at two fixed phi. Every n in this work
    is extracted as ln(E/E_pocket)/ln(b) from a cell whose b was set by Assur at
    the deck's SLAB brine fraction while phi was read back as the realised total,
    which includes the pocket population and runs about 0.019 higher. The four
    LCOL cells therefore sit off Assur's curve by construction: at phi = 0.099
    the cell carries b = 0.368 where Assur asks for 0.296. That is harmless if n
    does not depend on b and fatal if it does, and nothing has tested which.
    Solve this first.

RAMP -- four new brine fractions through the window, at realised phi = 0.092,
    0.096, 0.104 and 0.110, joining the existing 0.099. Five points across
    [0.09, 0.11] measure n(phi) rather than inverting one deficit into a
    linear form, and separate the ramp from the step by measurement.

SUBC -- realised phi = 0.075, 0.082 and 0.088, below phi_c. The closure sets
    w = 0 there on the geometric argument that no plane spans the material, and
    Section 4.5.1 shows what applying b unweighted would cost: 6.31 GPa against
    a measured 7.18-8.60 at the cold surface. The argument is sound and the
    turn-off is still an assertion. A cell with layers built in at phi = 0.075
    either returns n near zero, and the weight is confirmed where it matters
    most, or it does not, and phi_c is in the wrong place. This deck also
    carries LCOL p060, which failed to mesh: at its thin layers and wide bridges
    a pocket straddling a plane leaves a sliver Gmsh cannot repair. Every cell
    in both new decks is at or below that layer thickness, so all of them take
    the wider min_distance described below.

Two departures from make_layer_decks.py, both stated because they make the new
cells not quite twins of the old.

    1. phi is TARGETED ON THE REALISED TOTAL. The pocket population adds about
       0.019 of brine on top of the slab fraction -- measured as +0.0191,
       +0.0195, +0.0189, +0.0172 across the four LCOL cells -- so nominal
       slab_vof is set to the target less that offset. LCOL's own nominal
       values were the slab fractions and its realised phi came out where it
       came out; hitting the window needs the offset carried explicitly.
    2. b IS ASSUR AT THE TARGET REALISED phi, not at the slab fraction. This
       puts the cells on the curve the closure evaluates, which is what makes
       them a test of it. It is also why RAMPB comes first: pooling these with
       the LCOL cells is only legitimate if n is b-independent.

    python3 make_ramp_decks.py <outdir>
"""
import csv
import os
import sys

from make_layer_decks import (BASE, COLS, K, assur_b, mesh_for, thickness)

# Brine the pocket population adds on top of the slab fraction, from the four
# LCOL cells (+0.0191, +0.0195, +0.0189, +0.0172 at nominal 0.08 to 0.15).
# Flat to a thousandth across the range, as it should be -- the packing is the
# same in every cell of the deck.
POCKET_OFFSET = 0.019

PHI_C = 0.09
# min_distance is opened up across the whole campaign, and the criterion is the
# layer thickness rather than the brine fraction. LCOL meshed at 0.002 down to
# t = 0.0158 and failed at t = 0.0137, a pocket straddling a plane leaving a
# sliver facet Gmsh reports as an overlapping boundary and cannot repair. Every
# cell here runs t = 0.0114 to 0.0153, so all of them sit at or below the
# thickness that failed and all of them get the wider gap. Keeping it uniform
# also means the pocket packing is identical across the campaign, which matters
# because the offset below is what puts phi where it is aimed.
THIN_T = 0.0155


def deck(path, targets, tag, seeds=(1, 2, 3), states=('drn', 'und')):
    """One row per (target phi, drainage state, seed).

    targets are REALISED brine fractions; slab_vof is backed out through
    POCKET_OFFSET and b is Assur's value at the target itself.
    """
    rows = []
    with open(path, 'w', newline='', encoding='utf8') as fh:
        w = csv.writer(fh)
        w.writerow(COLS)
        for phi in targets:
            slab = round(phi - POCKET_OFFSET, 4)
            b = assur_b(phi)
            t = thickness(slab, b)
            lm = mesh_for(slab, b)
            md = '0.005' if t < THIN_T else '0.002'
            for state in states:
                for s in seeds:
                    r = dict(BASE)
                    r['run_id'] = '%s_p%03d_%s_s%d' % (
                        tag, round(phi * 1000), state, s)
                    r['slab_vof'] = '%.4f' % slab
                    r['bridge_fraction'] = '%.4f' % b
                    r['K_inclusion'] = K[state]
                    r['L_mesh'] = '%.4f' % lm
                    r['min_distance'] = md
                    w.writerow([r[c] for c in COLS])
            rows.append((phi, slab, b, t, lm, md))
    print('wrote %s  (%d cells)' % (path, len(rows) * len(states) * len(seeds)))
    print('  %-8s %-8s %-8s %-9s %-8s %-8s %s'
          % ('phi', 'slab', 'b Assur', 't', 'L_mesh', 'min_d', 'el across'))
    for phi, slab, b, t, lm, md in rows:
        print('  %-8.3f %-8.4f %-8.4f %-9.5f %-8.5f %-8s %.1f'
              % (phi, slab, b, t, lm, md, t / lm))


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else '.'
    os.makedirs(out, exist_ok=True)

    print('the window: phi_c = %.3f to phi_sat = 0.104, one existing cell at '
          '0.099\n' % PHI_C)

    # Inside and just above the ramp. 0.104 is the endpoint itself, where the
    # weight is claimed to be fully on and n should have reached the plateau;
    # 0.110 is the first point past it and tests that the plateau is flat
    # rather than still climbing toward the 0.119 cell.
    deck(os.path.join(out, 'rve_rampn.csv'),
         (0.092, 0.096, 0.104, 0.110), 'RAMP')
    print()

    # Below the threshold, where the closure asserts w = 0. 0.088 sits just
    # under phi_c, 0.075 is LCOL p060's realised fraction.
    deck(os.path.join(out, 'rve_subc.csv'),
         (0.075, 0.082, 0.088), 'SUBC')

    print('\nRAMPB is params/rve_layerb.csv, already written and never solved;'
          '\nno new deck is needed for it.')
    print('\nWhat each deck can return that would change the closure:')
    print('  RAMPB  n drifting with b at fixed phi -- every exponent in')
    print('         Section 4.5.2 is then contaminated by the b-phi mismatch')
    print('  RAMP   n(phi) convex or stepped rather than linear -- Eq. (6)')
    print('         gets a measured width instead of two bracketing limits')
    print('  SUBC   n materially above zero below phi_c -- the branch does')
    print('         not turn off where Section 4.5.1 puts it')
    return 0


if __name__ == '__main__':
    sys.exit(main())
