#!/usr/bin/env python3
"""Figure: where the stress concentrates inside the cell, not just how much.

Sections 3.3 and 4.3.1 reduce the matrix stress field to percentiles and to an
m-norm. Those are the right summaries for ranking cases, but they discard the
spatial information, and the spatial information is what distinguishes the two
morphologies whose ranking depends on the Weibull modulus: gas voids raise the
stress a little almost everywhere, whereas a percolating channel network raises
it enormously in the thin ligaments between channels and hardly at all
elsewhere. That is why the two exchange places at m~4.

Panels (a)-(b) map the field. Elements are shown at their centroids in a slab
through the cell, coloured by the local concentration factor on a shared scale,
so the two cases can be read against each other directly.

Panel (c) gives the volume-weighted exceedance curve, the fraction of matrix
volume above a given concentration. This is the quantity a weakest-link
criterion integrates, and it separates the cases where the percentiles do not:
the curves cross, the channelled case carrying less moderately-stressed volume
but a far longer tail.

Run from results/:  python3 ../analysis/plot_scf_field.py [dump_dir]
"""
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as fs

fs.apply()
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

LABEL = {'BASE': 'warm base (channels + pockets + gas)',
         'CHAN': 'brine channels', 'GAS': 'gas voids',
         'CTRL': 'control (2% brine)', 'POCK': 'brine pockets',
         'ELON': 'elongated pockets'}


def load(d):
    out = {}
    for p in sorted(glob.glob(os.path.join(d, '*.npz'))):
        z = np.load(p, allow_pickle=True)
        if 'cent' not in z:
            continue
        tag = os.path.splitext(os.path.basename(p))[0].replace('WBL_', '').split('_s')[0]
        out[tag] = (z['scf'], z['vol'], z['cent'], float(z['L']))
    return out


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else 'scf_fields'
    data = load(d)
    if not data:
        raise SystemExit('no dumps with centroids in %s' % d)

    show = [t for t in ('BASE', 'CHAN', 'GAS', 'CTRL') if t in data][:2]
    if len(show) < 2:
        show = list(data)[:2]

    fig = plt.figure(figsize=(14.0, 4.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.25], wspace=0.42)

    # Shared colour scale, clipped at a high percentile of the pooled field:
    # the extreme tail is a handful of elements and would otherwise flatten
    # every map to a single colour.
    allscf = np.concatenate([data[t][0] for t in show])
    vmax = float(np.percentile(allscf, 99.5))
    norm = Normalize(vmin=1.0, vmax=vmax)

    for i, tag in enumerate(show):
        scf, vol, cent, L = data[tag]
        ax = fig.add_subplot(gs[0, i])
        # a slab through the middle of the cell, so interior structure is
        # visible rather than the outer surface
        zc = cent[:, 2]
        sel = (zc > 0.45 * L) & (zc < 0.55 * L)
        if sel.sum() < 200:
            sel = np.ones(len(zc), bool)
        o = np.argsort(scf[sel])          # draw hot elements last
        x = cent[sel][:, 0][o]; y = cent[sel][:, 1][o]; c = scf[sel][o]
        sc = ax.scatter(x, y, c=c, s=5, cmap='inferno', norm=norm, linewidths=0)
        ax.set_aspect('equal')
        ax.set_xlabel('$x$'); ax.set_ylabel('$y$' if i == 0 else '')
        ax.set_title('(%s) %s' % ('ab'[i], LABEL.get(tag, tag)), fontsize=13)
        ax.text(0.03, 0.97, 'P99 = %.2f\nmax = %.1f' % (np.percentile(scf, 99), scf.max()),
                transform=ax.transAxes, va='top', fontsize=11,
                bbox=dict(fc='white', ec='0.7', alpha=0.85))
    # Horizontal bar under the two maps. A vertical one sits between panel (b)
    # and panel (c) and collides with the latter's y-axis label.
    cb = fig.colorbar(sc, ax=fig.axes[:2], orientation='horizontal',
                      fraction=0.06, pad=0.18, aspect=45)
    cb.set_label('local SCF $=\\sigma_1/\\bar\\sigma_{11}$')

    # ---- (c) volume-weighted exceedance ---------------------------------
    ax = fig.add_subplot(gs[0, 2])
    COL = {'BASE': fs.VERM, 'CHAN': fs.BLUE, 'GAS': fs.ORANGE, 'CTRL': fs.SKY,
           'POCK': fs.PURPLE, 'ELON': fs.GREEN}
    for tag in ('BASE', 'CHAN', 'GAS', 'CTRL'):
        if tag not in data:
            continue
        scf, vol, cent, L = data[tag]
        o = np.argsort(scf)
        s = scf[o]; w = vol[o]
        frac = 1.0 - np.cumsum(w) / w.sum()      # volume fraction above s
        ax.semilogy(s, np.clip(frac, 1e-6, 1), color=COL.get(tag), lw=2.0,
                    label=LABEL.get(tag, tag))
    ax.axhline(0.01, color=fs.BLACK, ls='--', lw=1.0)
    ax.text(0.98, 0.30, 'P99', transform=ax.transAxes, ha='right', fontsize=11)
    ax.set_xlabel('local SCF'); ax.set_ylabel('matrix volume fraction above SCF')
    ax.set_ylim(1e-5, 1.05)
    ax.set_xlim(1.0, None)          # SCF<1 is unstressed matrix, not of interest
    ax.set_title('(c) volume-weighted exceedance', fontsize=13)
    ax.legend(fontsize=10, loc='upper right')

    for ext in ('png', 'pdf'):
        fig.savefig('scf_field.%s' % ext, dpi=200, bbox_inches='tight')
    print('wrote scf_field.{png,pdf}')
    for tag in data:
        scf, vol, cent, L = data[tag]
        above2 = float(vol[scf > 2].sum() / vol.sum())
        print('  %-5s n=%6d  P99=%5.2f  max=%6.2f  vol frac SCF>2 = %5.1f%%'
              % (tag, len(scf), np.percentile(scf, 99), scf.max(), 100 * above2))


if __name__ == '__main__':
    main()
