"""How wide is the ramp, and does the bridge branch really turn off below it?

Eq. (5) switches the bridge mechanism on with a weight that ramps linearly from
phi_c = 0.09 to phi_sat = 0.104, and Section 4.5.2 fixes that endpoint from one
cell. The four LCOL cells give

    n = ln(E_x / E_pocket(phi)) / ln(b)  =  0.66, 0.93, 1.04, 0.97
                                     at  phi = 0.099, 0.119, 0.139, 0.167

the upper three averaging 0.98 and the lowest read as a partial weight,
w = 0.66/0.98, which inverted puts phi_sat just above 0.10. Seed scatter on
that cell is small, so the number is stable; what one point cannot do is
distinguish the linear ramp from the step, and Section 4.5.1 quotes those two
as bracketing a transition width it does not measure.

This script does three things, and the first gates the other two.

1. GATE. n is extracted at the b the cell was BUILT with, while phi is read back
   as the realised total. Those two are not consistent in the LCOL deck: b was
   Assur's value at the slab fraction, and the pocket population adds ~0.019 on
   top, so the cell at phi = 0.099 carries b = 0.368 where Assur asks for 0.296.
   Pooling cells built on different conventions is only legitimate if n does not
   depend on b at fixed phi. rve_layerb.csv was written to test exactly that and
   has never been solved. If its results are present the gate is evaluated; if
   they are absent the script says so and treats the pooling as unverified
   rather than quietly assuming it.

2. THE RAMP. With the new cells at phi = 0.092, 0.096, 0.104, 0.110 the window
   holds five points instead of one. Three forms are fitted -- linear ramp with
   phi_sat free, step at phi_c, and a free power -- and compared on residual.
   The comparison can fail: if the residuals are indistinguishable the width
   stays unmeasured and the honest report is that five points were not enough.

3. THE TURN-OFF. Below phi_c the closure sets w = 0 on a geometric argument.
   The SUBC cells have layers built in at phi = 0.075 to 0.088 and measure what
   a layered cell actually returns there. n near zero confirms the weight where
   it matters most; n materially above zero puts phi_c in the wrong place.

Runs on whatever is present. With no new results it reports the current
one-point state, which is the baseline the campaign is against.

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
    """Fit the weight through the window and compare the candidate forms."""
    print('\n' + '=' * 72)
    print('2. THE RAMP -- n(phi) between phi_c and the plateau')
    print('=' * 72)
    print('  %-22s %7s %7s %8s %8s %7s %7s'
          % ('cell', 'phi', 'b', 'E GPa', 'sd', 'n', 'n sd'))
    for k, r in g.iterrows():
        print('  %-22s %7.4f %7.4f %8.3f %8.4f %7.3f %7.4f'
              % (k, r.phi, r.b, r.E, r.E_sd, r.n, r.n_sd))

    plateau = g[g.phi > 0.11]
    if plateau.empty:
        print('\n  no cells above phi = 0.11; the plateau is what defines the')
        print('  weight, so nothing can be fitted')
        return
    n_pl = plateau.n.mean()
    print('\n  plateau (phi > 0.11): n = %.3f over %d conditions, range %.2f-%.2f'
          % (n_pl, len(plateau), plateau.n.min(), plateau.n.max()))

    win = g[(g.phi >= PHI_C) & (g.phi <= 0.115)]
    print('  cells in [%.3f, 0.115]: %d' % (PHI_C, len(win)))
    if len(win) < 2:
        print('\n  ONE POINT (or none). phi_sat is an inversion, not a fit:')
        for k, r in win.iterrows():
            w = r.n / n_pl
            print('    %s: w = %.3f/%.3f = %.3f -> phi_sat = %.4f'
                  % (k, r.n, n_pl, w, PHI_C + (r.phi - PHI_C) / w))
        print('  Nothing here distinguishes the linear ramp from the step.')
        return

    w_obs = (win.n / n_pl).values
    phi = win.phi.values

    def sse_linear(phi_sat):
        w = np.clip((phi - PHI_C) / (phi_sat - PHI_C), 0.0, 1.0)
        return float(np.sum((w_obs - w) ** 2))

    grid = np.linspace(0.091, 0.20, 2000)
    s = np.array([sse_linear(x) for x in grid])
    phi_sat = float(grid[int(np.argmin(s))])
    sse_lin = float(s.min())
    sse_step = float(np.sum((w_obs - 1.0) ** 2))

    # A free power on the same support: w = ((phi-phi_c)/(phi_sat-phi_c))^q.
    # q < 1 is a fast take-up (closer to the step), q > 1 a slow one. One extra
    # parameter against four or five points, so it is read as a direction and
    # not as a measurement of q.
    best = (np.inf, None, None)
    for x in np.linspace(0.095, 0.25, 400):
        u = np.clip((phi - PHI_C) / (x - PHI_C), 0.0, 1.0)
        for q in np.linspace(0.2, 3.0, 200):
            v = float(np.sum((w_obs - u ** q) ** 2))
            if v < best[0]:
                best = (v, x, q)

    print('\n  observed weight w = n/n_plateau:')
    for p, w in zip(phi, w_obs):
        print('    phi = %.4f   w = %.3f' % (p, w))
    print('\n  %-28s %10s %10s' % ('form', 'phi_sat', 'SSE'))
    print('  %-28s %10.4f %10.5f' % ('linear ramp (fitted)', phi_sat, sse_lin))
    print('  %-28s %10s %10.5f' % ('step at phi_c', '-', sse_step))
    print('  %-28s %10.4f %10.5f  q = %.2f'
          % ('free power', best[1], best[0], best[2]))
    print('\n  published value: phi_sat = %.3f' % PHI_SAT_PUBLISHED)
    if sse_step < sse_lin:
        print('  the STEP fits better: the transition is sharper than Eq. (6)')
    elif sse_lin < 0.5 * sse_step:
        print('  the ramp is resolved and beats the step by %.1fx in residual'
              % (sse_step / max(sse_lin, 1e-12)))
    else:
        print('  the two forms are not separated by these points; the width')
        print('  stays bracketed rather than measured')


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


def consequences(phi_sat_values=(0.095, PHI_SAT_PUBLISHED, 0.12, 0.16)):
    """What the window is worth downstream, so the campaign has a stake."""
    try:
        import ez_closure as ez
    except ImportError:
        return
    print('\n' + '=' * 72)
    print('WHAT THE WINDOW IS WORTH')
    print('=' * 72)
    keep = ez.PHI_SAT
    probes = (0.090, 0.095, 0.099, 0.104, 0.110, 0.120)
    print('  E(phi), GPa, against the ramp endpoint')
    print('  %7s' % 'phi', ''.join('%10s' % ('sat=%.3f' % s)
                                   for s in phi_sat_values))
    for p in probes:
        row = []
        for s in phi_sat_values:
            ez.PHI_SAT = s
            row.append(ez.E_of_phi(p))
        print('  %7.3f' % p, ''.join('%10.3f' % v for v in row))

    # Gogolaze beam 3: their eq. (14) brine profile over a 0.32 m section.
    H, B = 0.32, 0.60
    z = np.linspace(0, H, 400)
    zc = z * 100.0
    vb = (0.29315 * zc ** 2 - 5.124 * zc + 85.977) / 1000.0
    print('\n  Gogolaze beam 3: phi = %.4f to %.4f, %.0f%% of the depth above'
          ' phi_c' % (vb.min(), vb.max(), 100 * np.mean(vb >= PHI_C)))
    print('  %10s %10s %12s %10s' % ('phi_sat', 'z_NA/H', 'EI (MN m2)', 'E_app GPa'))
    for s in phi_sat_values:
        ez.PHI_SAT = s
        E = np.array([ez.E_of_phi(p) for p in vb]) * 1e9
        A = np.trapz(E, z) * B
        zn = np.trapz(E * z, z) * B / A
        EI = np.trapz(E * (z - zn) ** 2, z) * B
        print('  %10.3f %10.3f %12.4f %10.3f'
              % (s, zn / H, EI / 1e6, EI / (B * H ** 3 / 12) / 1e9))
    ez.PHI_SAT = keep


def main():
    frames = [f for f in (load(d, r) for d, r in DECKS) if f is not None]
    if not frames:
        print('no layered results found under results/')
        return 1
    have = sorted(set(pd.concat(frames).deck))
    print('decks present: %s\n' % ', '.join(have))
    gate()
    _, g = exponents(pd.concat(frames, ignore_index=True))
    ramp(g)
    turnoff(g)
    consequences()
    return 0


if __name__ == '__main__':
    sys.exit(main())
