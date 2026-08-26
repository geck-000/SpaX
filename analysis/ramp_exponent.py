"""Does the bridge branch really turn off below phi_c, and did the step survive?

Eq. (5) used to switch the bridge mechanism on with a weight ramping linearly
from phi_c = 0.09 to phi_sat = 0.104, with that endpoint fixed from one cell:
the four LCOL cells gave n = 0.66, 0.93, 1.04, 0.97, the lowest read as a
partial weight w = 0.66/0.98 which inverted put phi_sat just above 0.10.

THE RAMP IS RETRACTED. phi_sat came from the one cell whose two bridges
percolate at b ~ 0.31, where the pair and their periodic images make and break
a connected ice path. At four bridges that deficit is gone: seven N=4 cells
spanning b = 0.388 to 0.311 hold n to a range of 0.079 where the two-bridge
step is 0.256, and the phi_sat cell itself returns n = 0.583 with nothing left
to invert. Four of those cells run straight through phi_c from 0.076 to 0.093
with no feature there. Of the two forms the paper compares, only the step
retains support. Nothing in this script's ramp fit should therefore be read as
a measurement of a transition width.

This script now does three things, and the first gates the other two.

1. GATE. n is extracted at the b the cell was BUILT with, while phi is read back
   as the realised total. Those two are not consistent in the LCOL deck: b was
   Assur's value at the slab fraction, and the pocket population adds ~0.019 on
   top, so the cell at phi = 0.099 carries b = 0.368 where Assur asks for 0.296.
   Pooling cells built on different conventions is only legitimate if n does not
   depend on b at fixed phi. rve_layerb.csv was written to test exactly that and
   has never been solved. If its results are present the gate is evaluated; if
   they are absent the script says so and treats the pooling as unverified
   rather than quietly assuming it.

2. THE STEP. With the four-bridge cells at phi = 0.076 to 0.093 the window
   around phi_c holds no feature to measure: the step-vs-ramp comparison is
   resolved in the paper by direct measurement (the step survives, the ramp
   does not), so this section reports the N=4 data and no longer inverts
   phi_sat.

3. THE TURN-OFF. Below phi_c the closure sets w = 0 on a geometric argument.
   The SUBC cells have layers built in at phi = 0.075 to 0.088 and measure what
   a layered cell actually returns there. n near zero confirms the weight where
   it matters most; n materially above zero puts phi_c in the wrong place.

Runs on whatever is present. With no new results it reports the current
state, which is the baseline the campaign is against.

    python3 analysis/ramp_exponent.py
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

E_ICE = 9.37
PHI_0 = 0.20
PHI_C = 0.09
DRAIN_FACTOR = 1.04          # pocket law as calibrated is undrained; Sec. 4.4
PHI_SAT_PUBLISHED = 0.104

DECKS = [('rve_layercol', 'results_layercol'),
         ('rve_rampn', 'results_rampn'),
         ('rve_subc', 'results_subc')]


def assur_b(phi, phi_0=PHI_0):
    return np.clip(1.0 - np.sqrt(np.clip(phi, 0.0, phi_0) / phi_0), 0.0, 1.0)


def e_pocket(phi):
    """The drained pocket law, GPa. Eq. (2) divided by the drainage factor."""
    return E_ICE * (1.0 - 1.65 * np.asarray(phi, float)) / DRAIN_FACTOR


def load(deck, result):
    """Join a result set to its deck so each cell carries the b it was built at.

    b is not in the result file and cannot be recovered from it: the generator
    takes it as an input and the postprocessor reports moduli. Without the join
    every exponent here would be computed against an assumed b, which is the
    thing this script exists to avoid.
    """
    rp = os.path.join(ROOT, 'results', result + '.csv')
    dp = os.path.join(ROOT, 'params', deck + '.csv')
    if not (os.path.exists(rp) and os.path.exists(dp)):
        return None
    r = pd.read_csv(rp)
    d = pd.read_csv(dp)[['run_id', 'slab_vof', 'bridge_fraction', 'L_mesh']]
    m = r.merge(d, on='run_id', suffixes=('', '_deck'))
    m = m.dropna(subset=['E_x', 'phi_inclusion'])
    m['drained'] = m.run_id.str.contains('_drn_')
    m['deck'] = deck
    return m


# Which b convention a deck was built on. The distinction is not cosmetic: n is
# extracted as ln(E/E_pocket)/ln(b), so a b that is too large makes ln(b) less
# negative and DEPRESSES n. LCOL took b from Assur at the slab fraction while phi
# was read back as the realised total, ~0.019 higher, so every LCOL cell carries
# a b above the Assur curve and an n biased low. At LCOL_p080 that is b = 0.368
# against the 0.296 the curve asks for -- and that cell is the one phi_sat was
# derived from, its low n read as the weight not yet being fully on.
#
# RAMP and SUBC take b from Assur at the realised target, so they sit ON the
# curve. Pooling the two conventions is only legitimate if n does not depend on
# b, which is what the LAYERB gate decides. Until it does, they are reported
# apart.
ON_CURVE = ('rve_rampn', 'rve_subc')


def exponents(df):
    """n per cell, and the same aggregated over the seed replicates."""
    d = df[df.drained].copy()
    d['phi'] = d.phi_inclusion
    d['b'] = d.bridge_fraction
    d['E'] = d.E_x / 1e9
    d['n'] = np.log(d.E / e_pocket(d.phi)) / np.log(d.b)
    key = d.run_id.str.rsplit('_s', n=1).str[0]
    g = d.groupby(key.values).agg(
        deck=('deck', 'first'), phi=('phi', 'mean'), b=('b', 'first'),
        E=('E', 'mean'), E_sd=('E', lambda s: s.std(ddof=0)),
        n=('n', 'mean'), n_sd=('n', lambda s: s.std(ddof=0)),
        seeds=('n', 'size'))
    return d, g.sort_values('phi')


def gate(verbose=True):
    """Does n depend on b at fixed phi? Returns None if the cells are unsolved."""
    lb = load('rve_layerb', 'results_layerb')
    print('=' * 72)
    print('1. GATE -- is n independent of b at fixed phi?')
    print('=' * 72)
    if lb is None or lb.empty:
        print('results_layerb.csv absent: rve_layerb.csv (48 cells) has never')
        print('been solved. The LCOL cells carry b = Assur(slab) while phi is')
        print('the realised total ~0.019 higher, so at phi = 0.099 the cell')
        print('has b = 0.368 against Assur\'s 0.296 -- pooling decks built on')
        print('different conventions is UNVERIFIED until this deck runs.')
        return None
    _, g = exponents(lb)
    g = g.assign(phi_nom=np.round(g.phi, 2))
    ok = True
    for p, s in g.groupby(g.index.str.extract(r'p(\d{3})')[0].values):
        rng = s.n.max() - s.n.min()
        print('  slab phi %s: b = %s' % (p, np.round(sorted(s.b), 3)))
        print('               n = %s   spread %.3f'
              % (np.round(s.n.values, 3), rng))
        if rng > 0.10:
            ok = False
    print('\n  n %s independent of b to within 0.10'
          % ('IS' if ok else 'is NOT'))
    if not ok:
        print('  -> every exponent in Section 4.5.2 is contaminated by the')
        print('     b-phi mismatch, and the decks cannot be pooled as they are')
    return ok


def ramp(g):
    """Report the four-bridge cells through phi_c, with the retraction noted."""
    print('\n' + '=' * 72)
    print('2. THE STEP -- n(phi) around phi_c, no width to invert')
    print('=' * 72)
    print('  phi_sat = %.3f is RETRACTED. It was read off the one cell whose'
          % PHI_SAT_PUBLISHED)
    print('  two bridges percolate at b ~ 0.31; at four bridges that cell shows')
    print('  no deficit to invert and no feature marks the window. Only the')
    print('  step retains support, so nothing here fits a transition width.')
    print('  %-22s %7s %7s %8s %8s %7s %7s'
          % ('cell', 'phi', 'b', 'E GPa', 'sd', 'n', 'n sd'))
    for k, r in g.iterrows():
        print('  %-22s %7.4f %7.4f %8.3f %8.4f %7.3f %7.4f'
              % (k, r.phi, r.b, r.E, r.E_sd, r.n, r.n_sd))

    plateau = g[g.phi > 0.11]
    if not plateau.empty:
        print('\n  plateau (phi > 0.11): n = %.3f over %d conditions, range %.2f-%.2f'
              % (plateau.n.mean(), len(plateau),
                 plateau.n.min(), plateau.n.max()))
    win = g[(g.phi >= PHI_C) & (g.phi <= 0.115)]
    if not win.empty:
        print('  cells in [%.3f, 0.115]: %d (reported; no width is fitted)'
              % (PHI_C, len(win)))


def turnoff(g):
    print('\n' + '=' * 72)
    print('3. THE TURN-OFF -- what a layered cell returns below phi_c')
    print('=' * 72)
    sub = g[g.phi < PHI_C]
    if sub.empty:
        print('  no cells below phi_c = %.3f. The closure sets w = 0 there on' % PHI_C)
        print('  the geometric argument of Section 4.5.1 and nothing measures')
        print('  it; rve_subc.csv is the deck that would.')
        return
    plateau = g[g.phi > 0.11]
    n_pl = plateau.n.mean() if not plateau.empty else np.nan
    print('  %-22s %7s %7s %8s %7s %7s' % ('cell', 'phi', 'b', 'E GPa', 'n', 'w'))
    for k, r in sub.iterrows():
        print('  %-22s %7.4f %7.4f %8.3f %7.3f %7.3f'
              % (k, r.phi, r.b, r.E, r.n, r.n / n_pl))
    w = (sub.n / n_pl).values
    if np.nanmax(w) < 0.15:
        print('\n  w stays under 0.15 below phi_c: the weight is confirmed')
        print('  where Section 4.5.1 needed it, at the cold end.')
    else:
        print('\n  w reaches %.2f below phi_c. The bridge mechanism has not' % np.nanmax(w))
        print('  turned off where the closure says it has, and phi_c is in')
        print('  the wrong place or the weight has the wrong support.')


def control():
    """One increment against ten, on bit-identical meshes.

    The campaign solves nlgeom-OFF cells in a single increment because the
    extractor reads the last frame only and the response is proportional to the
    imposed displacement. That is an argument about linearity, and arguments
    about linearity are exactly the ones worth checking against a solve. The
    RAMPC decks are copies of the phi = 0.104 decks with the increment line
    edited back, so the mesh, the packing and the periodic equations are
    identical and the increment size is the only difference between them.
    """
    print('\n' + '=' * 72)
    print('CONTROL -- one increment against ten')
    print('=' * 72)
    a = load('rve_rampn', 'results_rampn')
    c = load('rve_rampctl', 'results_rampctl')
    if a is None or c is None or c.empty:
        print('  results_rampctl.csv absent; the single-increment solve is')
        print('  argued and not yet measured.')
        return
    a = a[a.run_id.str.startswith('RAMP_p104')].set_index(
        a.run_id.str.replace('RAMP_', '', regex=False))
    c = c.set_index(c.run_id.str.replace('RAMPC_', '', regex=False))
    both = sorted(set(a.index) & set(c.index))
    if not both:
        print('  no matching run_ids between the two sets')
        return
    print('  %-18s %12s %12s %10s' % ('cell', '1 inc', '10 inc', 'rel diff'))
    worst = 0.0
    for k in both:
        for col in ('E_x', 'E_z'):
            e1, e10 = float(a.loc[k, col]), float(c.loc[k, col])
            r = abs(e1 - e10) / e10
            worst = max(worst, r)
            print('  %-18s %12.6g %12.6g %10.2e' % (k + ':' + col, e1, e10, r))
    print('\n  largest relative difference: %.2e' % worst)
    if worst < 1e-6:
        print('  the two are the same solve; the nine discarded frames were')
        print('  discardable and the saving costs nothing.')
    else:
        print('  NOT identical. Something in these cells is not linear, and')
        print('  the single-increment results are not interchangeable with')
        print('  the LCOL cells they are pooled with.')


def consequences():
    """The downstream stake is the n(b) band now, not the retracted window."""
    try:
        import ez_closure as ez
    except ImportError:
        return
    print('\n' + '=' * 72)
    print('WHAT THE BAND IS WORTH')
    print('=' * 72)
    probes = (0.090, 0.095, 0.099, 0.104, 0.110, 0.120)
    print('  E(phi), GPa, n(b) and its +/- 2 rms band')
    print('  %7s %10s %10s %10s' % ('phi', 'lo', 'mid', 'hi'))
    for p in probes:
        lo, mid, hi = (float(v) for v in ez.E_band(p))
        print('  %7.3f %10.3f %10.3f %10.3f' % (p, lo, mid, hi))


def main():
    frames = [f for f in (load(d, r) for d, r in DECKS) if f is not None]
    if not frames:
        print('no layered results found under results/')
        return 1
    have = sorted(set(pd.concat(frames).deck))
    print('decks present: %s\n' % ', '.join(have))
    gate()
    allf = pd.concat(frames, ignore_index=True)
    on = allf[allf.deck.isin(ON_CURVE)]
    off = allf[~allf.deck.isin(ON_CURVE)]

    print(chr(10) + '#' * 72)
    print('# CELLS ON THE ASSUR CURVE (b from the realised phi): RAMP + SUBC')
    print('#' * 72)
    if not on.empty:
        _, g_on = exponents(on)
        ramp(g_on)
        turnoff(g_on)

    if not off.empty:
        print(chr(10) + '#' * 72)
        print('# CELLS OFF THE CURVE (b from the slab fraction): LCOL')
        print('# Reported apart, not pooled. Their n is biased LOW by the')
        print('# b convention -- see ON_CURVE above -- and the size of that')
        print('# bias is what the LAYERB gate has to bound before these can')
        print('# be combined with the set above.')
        print('#' * 72)
        _, g_off = exponents(off)
        print('  %-22s %7s %7s %8s %7s' % ('cell', 'phi', 'b', 'E GPa', 'n'))
        for k, r in g_off.iterrows():
            print('  %-22s %7.4f %7.4f %8.3f %7.3f' % (k, r.phi, r.b, r.E, r.n))
    g = g_on if not on.empty else g_off
    control()
    consequences()
    return 0


if __name__ == '__main__':
    sys.exit(main())
