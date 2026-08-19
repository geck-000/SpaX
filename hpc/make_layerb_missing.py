#!/usr/bin/env python3
r"""The LAYERB cells that out_lb does not already hold, and why each is missing.

An abandoned generation pass left 31 of the 48 LAYERB cells meshed on scratch as
complete utx+utz pairs. Those are reused rather than rebuilt -- generation is the
slow half of this campaign, and rebuilding them would buy nothing. This script
writes the deck for the other 17, and it does not treat them as one group,
because the logs say they failed for two different reasons and only one of them
is a knob.

TIMEOUT (5 cells: p080 b010 drn_s3 / und_s1 / und_s2, p080 b020 und_s2 / und_s3).
    "mesh subprocess exceeded 900s (degenerate geometry -- stuck in mesher)".
    These are the finest cells in the deck: b = 0.10 leaves the thinnest layer,
    mesh_for clamps L_mesh to the 0.005 floor, and the cell carries a million
    elements to the edge. Nothing is degenerate about them -- the mesher simply
    needed longer than SPAX_MESH_TIMEOUT allows by default. Raising it is the
    whole fix, and it costs nothing when a cell does not need it.

GEOMETRY (12 cells: every b = 0.45 cell at both brine fractions).
    "Invalid boundary mesh (overlapping facets)" and "periodic pair(s)
    unmatched", after six packings with sliver rejection already engaged. That
    is not a timeout and will not yield to one. At b = 0.45 two bridges carry
    45% of the plane between them, so each is a wide patch, and a pocket landing
    against one leaves the sliver facet Gmsh reports and cannot repair. The same
    fix that rescued LCOL p060 applies: hold the pockets off the planes with a
    wider min_distance.

    It may still not be enough, and the campaign is designed so that it need not
    be. LAYERB exists to test whether n depends on b, and the b values that
    matter are the ones the other decks actually sit at: the LCOL cells run
    b = 0.134 to 0.368 and the new RAMP and SUBC cells b = 0.258 to 0.388. Those
    are bracketed by b = 0.10, 0.20 and 0.30, all of which meshed. Losing b =
    0.45 costs the upper bracket point and leaves a factor of three in b still
    swept, which is enough to answer the question that gates the pooling.

    So these 12 are worth one attempt and not a campaign of retries. If they
    fail again, the honest report is that the sweep spans 0.10 to 0.30.

    python3 make_layerb_missing.py <have_list> <outdir>

<have_list> is one run_id per line, the cells already generated -- produced on
the cluster by

    ls out_lb/*-utx.inp | sed 's|.*/Job-||; s|-utx.inp||' | sort > have_lb.txt

Deriving it from the filesystem rather than hardcoding it means a partial rerun
of this deck narrows the next one automatically.
"""
import csv
import os
import sys

# Pockets held off the layer planes, as for LCOL p060 and the SUBC deck. Applied
# only to the cells whose failure was geometric; the timeout cells meshed their
# geometry fine and changing their packing would make them differ from the 31
# reused cells for no reason.
WIDE_MIN_DISTANCE = '0.005'


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    have_path = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else '.'
    here = os.path.dirname(os.path.abspath(__file__))
    deck = os.path.join(os.path.dirname(here), 'params', 'rve_layerb.csv')
    if not os.path.exists(deck):
        deck = 'rve_layerb.csv'

    have = set()
    if os.path.exists(have_path):
        with open(have_path) as fh:
            have = set(l.strip() for l in fh if l.strip())

    with open(deck, newline='', encoding='utf8') as fh:
        rows = list(csv.DictReader(fh))
        cols = list(rows[0].keys())

    missing = [r for r in rows if r['run_id'] not in have]
    if not missing:
        print('nothing missing: all %d cells are already generated' % len(rows))
        return 0

    n_geom = 0
    for r in missing:
        # b = 0.45 is the geometric failure; everything else was a timeout.
        if float(r['bridge_fraction']) >= 0.40:
            r['min_distance'] = WIDE_MIN_DISTANCE
            n_geom += 1

    p = os.path.join(out, 'rve_layerb_missing.csv')
    with open(p, 'w', newline='', encoding='utf8') as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in missing:
            w.writerow(r)

    print('deck holds %d cells; out_lb already has %d'
          % (len(rows), len(rows) - len(missing)))
    print('wrote %s  (%d cells: %d geometric, %d timeout)'
          % (p, len(missing), n_geom, len(missing) - n_geom))
    print()
    print('  %-24s %8s %8s %10s %s' % ('run_id', 'b', 'L_mesh', 'min_dist',
                                       'why it is missing'))
    for r in missing:
        geom = float(r['bridge_fraction']) >= 0.40
        print('  %-24s %8s %8s %10s %s'
              % (r['run_id'], r['bridge_fraction'], r['L_mesh'],
                 r['min_distance'], 'overlapping facets' if geom
                 else 'mesher exceeded 900s'))
    print()
    print('Generate with SPAX_MESH_TIMEOUT well above the 900s default -- the')
    print('five timeout cells are the reason the deck is not already complete.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
