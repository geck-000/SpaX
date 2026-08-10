# -*- coding: utf-8 -*-
"""Does a layered brine topology give E(z) the steepness the profiles show?

The pocket cells are too flat. Across the Gogolaze column phi_soft moves from
0.085 to 0.134 and the computed E_x moves from 7.98 to 7.15 GPa -- a 10% spread
against measurements that fall by a factor of several toward the base. The
knockdown law behind that, E/E_ice = 1 - 1.65 phi, is linear and gentle because
the matrix percolates around every pocket, so no pocket shape can steepen it.

A layer topology is steep for free. If the brine sits in cell-spanning planes
pierced by ice bridges of area fraction b, the transverse modulus is a series
sum dominated by the bridges,

    E_x = E_ice / [ (1 - t) + t / b ],   t = layer thickness fraction,

and t/b blows up as b falls. The question is whether b falls with phi_b the way
it has to, WITHOUT b being fitted to the moduli it is supposed to predict.

The closure used here is Assur's plane-of-weakness geometry, which is where the
sea-ice literature's sqrt(phi) law comes from: brine occupies the plane between
ice platelets and the load-bearing area left over is

    b = 1 - sqrt(phi_b / phi_0),   phi_0 ~ 0.20.

phi_0 is Assur's constant, taken from the strength literature and not adjusted
here. That it was calibrated for STRENGTH and is being applied to MODULUS is a
real assumption and is flagged in the output, not hidden.
"""
import csv
import os
import re
import sys

import numpy as np

E_ICE = 9.37
PHI_0 = 0.20


def bridge_fraction(phi, phi_0=PHI_0):
    """Assur load-bearing area fraction left in the plane of weakness."""
    return np.clip(1.0 - np.sqrt(np.clip(phi, 0.0, phi_0) / phi_0), 1e-4, 1.0)


def E_layer(phi, phi_0=PHI_0, E_ice=E_ICE):
    """Transverse modulus of a cell whose brine is all in bridged layers."""
    b = bridge_fraction(phi, phi_0)
    t = np.clip(phi / np.maximum(1.0 - b, 1e-6), 0.0, 0.98)
    return E_ice / ((1.0 - t) + t / b), b, t


def E_pocket(phi, E_ice=E_ICE):
    """The knockdown law the pocket cells obey (R^2 = 0.999 over 50 cells)."""
    return E_ice * (1.0 - 1.65 * phi)


def load_column(path):
    depth, phi, Ex = [], [], []
    for r in csv.DictReader(open(path)):
        m = re.search(r'z(\d+)', r['run_id'])
        try:
            e = float(r['E_x'])
        except (ValueError, KeyError):
            continue
        if not m:
            continue
        depth.append(int(m.group(1)) / 1000.0)
        phi.append(float(r['phi_soft_total']))
        Ex.append(e / 1e9)
    return np.array(depth), np.array(phi), np.array(Ex)


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        'results', 'results_gogo_column.csv')
    d, phi, ex = load_column(path)
    if not len(d):
        print('no usable rows in %s' % path)
        return 1

    # collapse the seed replicates at each depth
    zs = sorted(set(d))
    print('Gogolaze column: %d depths, %d cells\n' % (len(zs), len(d)))
    print('%8s %9s %11s %11s %9s %9s' % (
        'z (m)', 'phi_soft', 'E pocket', 'E layer', 'b', 't'))
    print('%8s %9s %11s %11s %9s %9s' % (
        '', '', 'FE, GPa', 'pred, GPa', 'Assur', 'frac'))
    pk, lay, pk_fe = [], [], []
    for z in zs:
        s = d == z
        p = phi[s].mean()
        el, b, t = E_layer(p)
        print('%8.3f %9.4f %11.3f %11.3f %9.4f %9.4f'
              % (z, p, ex[s].mean(), el, b, t))
        pk.append(E_pocket(p)); lay.append(el); pk_fe.append(ex[s].mean())
    pk, lay, pk_fe = np.array(pk), np.array(lay), np.array(pk_fe)

    print('\nspread top-to-bottom over the column')
    print('  pocket cells (FE)      %.3f -> %.3f GPa   ratio %.2f'
          % (pk_fe.max(), pk_fe.min(), pk_fe.max() / pk_fe.min()))
    print('  pocket law             %.3f -> %.3f GPa   ratio %.2f'
          % (pk.max(), pk.min(), pk.max() / pk.min()))
    print('  layered prediction     %.3f -> %.3f GPa   ratio %.2f'
          % (lay.max(), lay.min(), lay.max() / lay.min()))

    print('\nsensitivity of E to phi_b, d(lnE)/d(ln phi) at the column mean')
    pm = phi.mean()
    h = 0.02 * pm
    for name, f in (('pocket law', lambda x: E_pocket(x)),
                    ('layered', lambda x: E_layer(x)[0])):
        g = (np.log(f(pm + h)) - np.log(f(pm - h))) / (
            np.log(pm + h) - np.log(pm - h))
        print('  %-20s %7.2f' % (name, g))

    print('\nbridge fractions this column asks for: %.3f (top) to %.3f (base)'
          % (bridge_fraction(phi.min()), bridge_fraction(phi.max())))
    print('the FE sweep covers b = 0.01 to 0.10, so the base is inside the')
    print('bracket and the upper column is not -- those cells need a sweep')
    print('at larger b before the profile can be claimed end to end.')
    print('\nCAVEAT: phi_0 = %.2f is Assur\'s constant from the STRENGTH'
          % PHI_0)
    print('literature. Using it for modulus is an assumption, not a result.')

    invert(np.array(zs), np.array([phi[d == z].mean() for z in zs]))
    return 0


# Kujala et al. (1990) Table 2, the four strain-gauged beams, endpoints meaned.
# Kujala's own reduction assumes E(z) linear and fits two parameters, so the
# straight line interpolated below is their model rather than sampled data.
# The monotonic fall itself is not an artefact of that choice: Marchenko and
# Gogolaze recover the same decreasing profile by different methods, so a model
# that returns a non-monotonic or flat E(z) is failing against measurement and
# not merely against one group's fitting form.
K_TOP, K_BOT = 8.05, 1.27


def invert(zs, phis):
    """What bridge fraction would the measured profile require?

    Fitting b to the moduli would prove nothing -- b was already shown to be a
    knife edge, so some b always exists. The test that can fail is whether the
    REQUIRED b(phi) has the shape of a physical law. If it tracks Assur's
    1 - sqrt(phi/phi_0) over the whole column with one phi_0, the closure is
    doing real work. If it wanders, the layer model is a two-parameter fit
    dressed as a mechanism.
    """
    from scipy.optimize import brentq

    print('\n' + '=' * 66)
    print('INVERSION: bridge fraction the measured profile demands')
    print('=' * 66)
    H = zs.max()
    print('%8s %9s %10s %10s %10s' % (
        'z (m)', 'phi_soft', 'E meas', 'b needed', 'b Assur'))
    need, assur = [], []
    for z, p in zip(zs, phis):
        Em = K_TOP + (K_BOT - K_TOP) * (z / H)

        def resid(b):
            t = min(p / max(1.0 - b, 1e-9), 0.98)
            return E_ICE / ((1.0 - t) + t / b) - Em

        try:
            b = brentq(resid, 1e-6, 0.999999)
        except ValueError:
            b = float('nan')
        ba = bridge_fraction(p)
        need.append(b); assur.append(ba)
        print('%8.3f %9.4f %10.3f %10.4f %10.4f' % (z, p, Em, b, ba))

    need, assur = np.array(need), np.array(assur)
    ok = np.isfinite(need)
    # One free constant: the phi_0 that best matches the required b, found by
    # least squares on b itself. If the closure has the right SHAPE this should
    # track across the column; the residual is what says whether it does.
    def sse(p0):
        return np.nansum((need[ok] - bridge_fraction(phis[ok], p0)) ** 2)
    grid = np.linspace(0.15, 1.5, 400)
    p0 = grid[int(np.argmin([sse(g) for g in grid]))]
    pred = bridge_fraction(phis[ok], p0)
    ss_res = np.sum((need[ok] - pred) ** 2)
    ss_tot = np.sum((need[ok] - need[ok].mean()) ** 2)
    print('\nbest-fit phi_0 = %.3f  (Assur\'s strength value is %.2f)'
          % (p0, PHI_0))
    print('R^2 of Assur shape against required b : %.4f'
          % (1.0 - ss_res / ss_tot if ss_tot > 0 else float('nan')))
    print('required b spans %.4f to %.4f' % (np.nanmin(need), np.nanmax(need)))


if __name__ == '__main__':
    sys.exit(main())
