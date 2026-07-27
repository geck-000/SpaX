#!/usr/bin/env python3
"""Test whether the channel generator is isotropic in the RVE plane.

The warm-base full-tensor replicates returned E_y/E_x = 1.012 +/- 0.006 with
all five packings on the same side of unity. Two readings are possible: the
channel network prefers an in-plane direction (a generator bias), or one cell
simply holds too few channels for the two in-plane directions to be equivalent
(a resolution limit). They call for different fixes -- a code change, or a
larger cell -- so the question is worth settling directly.

The test bypasses meshing and solving. It generates channel arrangements at the
base-slice parameters and measures the fabric tensor of the near-neighbour
vectors, F = <n (x) n> over pairs closer than 4 r_avg. For an in-plane isotropic
generator E[F_xx - F_yy] = 0; a preferred direction shows up as a non-zero mean.
Inclusions are omitted (empty octree) so that only the channel placement is
under test.

Result at the base parameters (700 runs, two seed streams): the mean fabric
anisotropy is zero within its standard error, while a cell holds only 3-5
channels and the per-realisation fabric scatter is ~0.7. The generator is
unbiased; the in-plane split is a resolution limit, and more packings will not
remove it -- a larger cell would.

Usage: python3 check_channel_isotropy.py [n_runs] [seed]
"""
import contextlib
import importlib.util
import io
import math
import os
import sys

import numpy as np

# base-slice channel parameters, from params/rve_basetensor_seeds.csv
L, VOF, R_AVG, R_STD, SEP, MAX_IT, L_MESH = 0.50, 0.060, 0.022, 0.005, 0.002, 200000, 0.033
CUTOFF = 4.0 * R_AVG


class _NoInclusions:
    """Octree stand-in: the channel-inclusion clearance never binds."""

    def query(self, centre, radius):
        return []


def _load_generator():
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, '..', 'SpaX_Standalone.py')
    spec = importlib.util.spec_from_file_location('spax_standalone', path)
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
    except SystemExit:          # module runs a CLI guard on import
        pass
    return mod.generate_channels


def fabric_anisotropy(P):
    """F_xx - F_yy over min-image near-neighbour unit vectors; None if too few."""
    if len(P) < 2:
        return None
    fxx = fyy = 0.0
    cnt = 0
    for i in range(len(P)):
        dx = P[i, 0] - P[i + 1:, 0]
        dy = P[i, 1] - P[i + 1:, 1]
        dx -= L * np.round(dx / L)
        dy -= L * np.round(dy / L)
        d = np.hypot(dx, dy)
        sel = (d < CUTOFF) & (d > 1e-12)
        if sel.any():
            fxx += ((dx[sel] / d[sel]) ** 2).sum()
            fyy += ((dy[sel] / d[sel]) ** 2).sum()
            cnt += int(sel.sum())
    return (fxx - fyy) / cnt if cnt else None


def main():
    n_runs = int(sys.argv[1]) if len(sys.argv) > 1 else 700
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 777
    gen = _load_generator()
    np.random.seed(seed)

    fab, counts = [], []
    for _ in range(n_runs):
        with contextlib.redirect_stdout(io.StringIO()):     # generator is chatty
            ch = gen(L, VOF, R_AVG, R_STD, SEP, MAX_IT, _NoInclusions(),
                     L_mesh=L_MESH, densify=True)
        # keep primaries only; the returned array also carries periodic copies
        P = np.array([c for c in ch if 0 <= c[0] < L and 0 <= c[1] < L])
        counts.append(len(P))
        a = fabric_anisotropy(P) if len(P) else None
        if a is not None:
            fab.append(a)

    fab = np.array(fab)
    counts = np.array(counts)
    sem = fab.std(ddof=1) / math.sqrt(len(fab))
    print('runs %d (seed %d) | channels per cell: mean %.1f, median %d, range %d-%d'
          % (n_runs, seed, counts.mean(), int(np.median(counts)),
             counts.min(), counts.max()))
    print('usable realisations (>=2 channels): %d' % len(fab))
    print('fabric anisotropy F_xx - F_yy: mean %+.4f  sd %.4f  SEM %.4f  t = %+.2f'
          % (fab.mean(), fab.std(ddof=1), sem, fab.mean() / sem))
    print('  95%% CI on the mean: [%+.4f, %+.4f]'
          % (fab.mean() - 1.96 * sem, fab.mean() + 1.96 * sem))
    print('  P(F_xx > F_yy) = %.3f   (isotropic generator gives 0.500)'
          % np.mean(fab > 0))
    verdict = ('no detectable in-plane bias -- the split is a resolution limit'
               if abs(fab.mean()) < 1.96 * sem else
               'IN-PLANE BIAS DETECTED in the channel generator')
    print('  verdict: %s' % verdict)
    return 0


if __name__ == '__main__':
    sys.exit(main())
