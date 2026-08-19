"""Can the layered cells be made smaller, and cheaper to generate?

Generation is the slow half of every layered campaign -- a RAMP cell spent over
an hour in the mesher -- so it is worth asking whether the cells need to be the
size they are. The arithmetic says the saving would be large, and the data says
not to take it.

THE ARITHMETIC. The element count goes as (L/L_mesh)^3, and mesh_for sets
L_mesh from the brine-layer thickness t = slab_vof*L/(n_slabs*(1-b)). Scaling L
and n_slabs together by k leaves t unchanged, so L_mesh is unchanged and the
count falls as k^3. It also leaves a0 = L/n_slabs unchanged, and a0 is the
lamellar spacing -- the quantity the exponent depends on through a0^0.69 and the
reason the cells are solved at 0.75 mm in the first place. So the shrink looks
free: same spacing, same resolution across the layer, an eighth of the elements
at half the edge.

    L=0.500 n=4   610562 elements   1.000x   36 pockets
    L=0.375 n=3   257581 elements   0.422x   15 pockets
    L=0.250 n=2    76320 elements   0.125x    4 pockets

THE DATA. It is not free, and the experiment has already been run.
rve_eringen_layer.csv sweeps exactly this scaling -- L = 0.24, 0.36, 0.48 at
n_slabs = 2, 3, 4, so a0 = 0.12 throughout -- with six seeds at each size, and
rve_eringen_layer_homog.csv is its matched phi=0 control at the same sizes and
the same element size. The control is what makes the sweep readable: Section
4.4.1's own finding is that a size sweep without one measures the extraction
bias as much as the material.

Here the control barely moves, 11.07 to 11.32 GPa across the range, so the bias
is about 2% and the trend is the material. Normalised by it, the layered cell at
n_slabs = 2 comes out 46% stiffer than at n_slabs = 4, and at n_slabs = 3 still
11% stiffer. Two lamellae is not a representative volume for this morphology.
The seed scatter says the same thing from the other side: 3.9% of the mean at
L = 0.24 against 1.0% at L = 0.48, because a cell that small holds four pockets
and the packing stops averaging.

WHAT THIS MEANS FOR THE RAMP CAMPAIGN. Two things, and the second is the one
that settles it.

First, L = 0.48 is itself only approaching the asymptote rather than sitting on
it -- the normalised ratios run 0.358, 0.273, 0.245 and are still falling, by
about 5% at the last step. The campaign's L = 0.5 carries that residual.

Second, and decisively, it carries it in common with everything it is compared
against. The four LCOL cells are L = 0.5 with n_slabs = 4, the reused LAYERB
decks are L = 0.5, and the RAMP and SUBC cells are L = 0.5. n(phi) is a
comparison ACROSS those cells, so a size bias they all share cancels out of it.
Shrinking some of them would convert a cancelling common-mode bias into a
differential one, and a 46% differential is fatal to a measurement whose whole
content is a factor of about 1.3 in n across a narrow window in phi.

So the cells stay at L = 0.5. The saving is real and it is not available here.

CAVEAT worth keeping in view: the sweep above is a BENDING measure, E_bending
from the Eringen decks, and the RAMP cells are uniaxial. A uniaxial size sweep
of the same morphology might well show a smaller effect, and nothing here
measures that. It would not change the decision -- common-mode cancellation is
the argument that does -- but it would be the thing to run before shrinking
cells in some future campaign where comparability is not at stake.

    python3 analysis/layered_cell_size.py
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def cost_table():
    """Elements and pocket count under the a0-preserving scaling."""
    slab, b, r, vof = 0.085, 0.2789, 0.030, 0.0325
    print('The a0-preserving shrink, on the RAMP phi=0.104 geometry')
    print('  %6s %8s %8s %9s %10s %9s %8s'
          % ('L', 'n_slabs', 'a0', 'L_mesh', 'elements', 'rel cost', 'pockets'))
    base = None
    for L, n in ((0.500, 4), (0.375, 3), (0.250, 2)):
        t = slab * L / (n * (1.0 - b))
        lm = min(max(t / 2.5, 0.005), 0.012)
        ne = (L / lm) ** 3
        npk = vof * L ** 3 / ((4.0 / 3.0) * np.pi * r ** 3)
        base = ne if base is None else base
        print('  %6.3f %8d %8.4f %9.5f %10.0f %9.3f %8.1f'
              % (L, n, L / n, lm, ne, ne / base, npk))


def size_sweep():
    """The measured size dependence, normalised by its matched control."""
    lp = os.path.join(ROOT, 'results', 'results_eringen_layer.csv')
    hp = os.path.join(ROOT, 'results', 'results_eringen_layer_homog.csv')
    if not (os.path.exists(lp) and os.path.exists(hp)):
        print('\nresults_eringen_layer{,_homog}.csv not found; sweep not shown')
        return
    d = pd.read_csv(lp)
    h = pd.read_csv(hp)
    d['Eb'] = d.E_bending / 1e9
    h['Eb'] = h.E_bending / 1e9
    lay = d.groupby('L').Eb.agg(['mean', lambda s: s.std(ddof=0), 'size'])
    lay.columns = ['mean', 'sd', 'n']
    ctl = h.groupby('L').Eb.mean()

    print('\nMeasured, a0 = 0.12 held constant, six seeds per size')
    print('  %6s %8s %10s %8s %10s %10s %9s'
          % ('L', 'n_slabs', 'layered', 'scatter', 'control', 'normalised', 'vs L=max'))
    common = [L for L in lay.index if L in ctl.index]
    ref = lay.loc[common[-1], 'mean'] / ctl[common[-1]]
    for L in common:
        m, sd = lay.loc[L, 'mean'], lay.loc[L, 'sd']
        nrm = m / ctl[L]
        print('  %6.2f %8d %10.4f %7.1f%% %10.4f %10.4f %+8.1f%%'
              % (L, int(round(L / 0.12)), m, 100 * sd / m, ctl[L], nrm,
                 100 * (nrm / ref - 1.0)))
    print('\n  control moves only %.1f%% across the range, so the trend is the'
          % (100 * (ctl[common[-1]] / ctl[common[0]] - 1.0)))
    print('  material and not the extraction.')
    print('  Two lamellae is not a representative volume for this morphology.')


def main():
    cost_table()
    size_sweep()
    print('\nConclusion: the cells stay at L = 0.5. The size bias at that edge')
    print('is shared by the LCOL, LAYERB, RAMP and SUBC cells alike, so it')
    print('cancels from n(phi); shrinking any of them would turn a cancelling')
    print('bias into a differential one.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
