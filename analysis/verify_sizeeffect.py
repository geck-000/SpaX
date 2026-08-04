#!/usr/bin/env python3
"""Recompute every number quoted in Section 4.2.1, and check it against the text.

Why this exists
---------------
Section 4.2.1 rests on two controls applied to the same bending sweep, and the
result depends on which of them is in force:

  1. the extraction bias, measured on homogeneous (phi=0) cells at the same six
     sizes and divided out point by point;
  2. the composition drift, the realised brine and void fractions admitted as
     covariates because the packer does not hit its target equally well in a
     small cell as in a large one.

The manuscript's covariate model is run with BOTH in force. That was not stated
in an earlier draft, and it matters: on the uncorrected modulus the identical
three-covariate model returns a coefficient on 1/L^2 of +0.005 GPa (p=0.78)
rather than +0.026 GPa (p=0.17). Neither is a size effect -- both are null, which
is the paper's conclusion either way -- but a reader reproducing the number from
the raw data would not recover the quoted one and would reasonably conclude the
paper was wrong. This script makes the distinction explicit and mechanical.

It is a verification harness, not an analysis: every quantity it computes is
compared with the value quoted in the manuscript, and it exits non-zero if any
of them has drifted. Run it after re-solving the sweep, or before submission.

    python3 verify_sizeeffect.py            # run from results/
    python3 verify_sizeeffect.py --verbose  # also dump the per-size table
"""
import argparse
import csv
import os
import sys
from collections import defaultdict

import numpy as np
from scipy import stats

D_INCL = 0.08          # mean inclusion diameter, sets L/d

# Every value Section 4.2.1 states, with the tolerance to which it is quoted.
# (label, quoted value, absolute tolerance)
QUOTED = {
    'n_realisations':      (33,     0),
    'baseline_rise_pct':   (4.8,    0.05),
    'channelled_end_pct':  (2.9,    0.05),
    'size_mean_spread_pct': (13.9,  0.05),
    'raw_slope':           (-0.020, 0.0005),
    'raw_slope_se':        (0.013,  0.0005),
    'raw_slope_p':         (0.14,   0.005),
    'corr_slope':          (-0.003, 0.0005),
    'corr_slope_se':       (0.014,  0.0005),
    'corr_slope_p':        (0.82,   0.005),
    'corr_slope_r2':       (0.002,  0.0005),
    'anova_p':             (0.065,  0.0005),
    'cov_min_pct':         (2.3,    0.05),
    'cov_max_pct':         (9.2,    0.05),
    'brine_vs_invL2_p':    (0.30,   0.005),
    'void_smallest':       (0.050,  0.0005),
    'void_largest':        (0.013,  0.0005),
    'void_target':         (0.012,  0.0005),
    'void_vs_invL2_r2':    (0.47,   0.005),
    'cov_corr_invL2':      (0.026,  0.0005),
    'cov_corr_invL2_se':   (0.018,  0.0005),
    'cov_corr_invL2_p':    (0.17,   0.005),
    'cov_corr_brine_p':    (0.005,  0.0005),
    'cov_corr_void_p':     (0.020,  0.0005),
    'cov_raw_invL2':       (0.005,  0.0005),
    'cov_raw_invL2_se':    (0.018,  0.0005),
    'cov_raw_invL2_p':     (0.78,   0.005),
}


def load(path, col='E_bending'):
    """Rows with a usable modulus, grouped by cell edge L."""
    g = defaultdict(list)
    for r in csv.DictReader(open(path)):
        try:
            v = float(r[col])
        except (ValueError, KeyError, TypeError):
            continue
        if v > 0:
            g[float(r['L'])].append(r)
    return g


def ols(y, *cols):
    """Least squares with intercept. Returns (beta, stderr, pvalue) per column,
    intercept first. Written out rather than pulled from statsmodels so the
    repository keeps a scipy-only dependency footprint."""
    X = np.column_stack([np.ones_like(y)] + list(cols))
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = len(y) - X.shape[1]
    s2 = float(resid @ resid) / dof
    se = np.sqrt(np.diag(np.linalg.inv(X.T @ X)) * s2)
    p = 2.0 * (1.0 - stats.t.cdf(np.abs(beta / se), dof))
    return beta, se, p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--verbose', action='store_true')
    ap.add_argument('--power', action='store_true',
                    help='also report the smallest couple-stress length scale '
                         'this sweep could have detected')
    ap.add_argument('--channelled', default='results_eringen.csv')
    ap.add_argument('--homogeneous', default='results_eringen_homog.csv')
    args = ap.parse_args()

    for f in (args.channelled, args.homogeneous):
        if not os.path.isfile(f):
            raise SystemExit('missing %s -- run from results/' % f)

    ch = load(args.channelled)
    ho = load(args.homogeneous)
    Ls = sorted(ch)

    # ---- the extraction baseline ----------------------------------------
    # Homogeneous cells carry no microstructure, so their size dependence is
    # the cube-versus-plate kinematics plus discretisation. Normalised to the
    # widest cell, which is taken as the asymptote.
    href = {L: float(np.mean([float(r['E_bending']) / 1e9 for r in ho[L]]))
            for L in ho}
    ref = href[max(href)]
    rel = {L: href[L] / ref for L in href}

    got = {}
    got['n_realisations'] = sum(len(ch[L]) for L in Ls)
    got['baseline_rise_pct'] = 100.0 * (max(href.values()) / min(href.values()) - 1.0)

    means = np.array([np.mean([float(r['E_bending']) / 1e9 for r in ch[L]])
                      for L in Ls])
    got['channelled_end_pct'] = 100.0 * (means[-1] / means[0] - 1.0)
    got['size_mean_spread_pct'] = 100.0 * (means.max() / means.min() - 1.0)

    # ---- slope on 1/L^2, before and after the baseline ------------------
    for tag, use_corr in (('raw', False), ('corr', True)):
        x, y = [], []
        for L in Ls:
            f = rel[L] if use_corr else 1.0
            for r in ch[L]:
                x.append(1.0 / L ** 2)
                y.append(float(r['E_bending']) / 1e9 / f)
        res = stats.linregress(np.array(x), np.array(y))
        got['%s_slope' % tag] = res.slope
        got['%s_slope_se' % tag] = res.stderr
        got['%s_slope_p' % tag] = res.pvalue
        if tag == 'corr':
            got['corr_slope_r2'] = res.rvalue ** 2

    groups = [np.array([float(r['E_bending']) / 1e9 for r in ch[L]]) / rel[L]
              for L in Ls]
    got['anova_p'] = float(stats.f_oneway(*groups).pvalue)
    covs = [100.0 * np.std(g) / np.mean(g) for g in groups]
    got['cov_min_pct'], got['cov_max_pct'] = min(covs), max(covs)

    # ---- composition actually realised in the mesh ----------------------
    # phi_inclusion and porosity are recovered from the meshed phase volumes;
    # VoF_incl_sphere / VoF_void_sphere in the same CSV are what was REQUESTED
    # of the packer and are identical across sizes, so they cannot be used here.
    x, yraw, brine, void = [], [], [], []
    for L in Ls:
        for r in ch[L]:
            x.append(1.0 / L ** 2)
            yraw.append(float(r['E_bending']) / 1e9)
            brine.append(float(r['phi_inclusion']))
            void.append(float(r['porosity']))
    x = np.array(x); yraw = np.array(yraw)
    brine = np.array(brine); void = np.array(void)
    ycorr = np.array([v / rel[L] for L in Ls for v in
                      [float(r['E_bending']) / 1e9 for r in ch[L]]])

    got['brine_vs_invL2_p'] = stats.linregress(x, brine).pvalue
    rv = stats.linregress(x, void)
    got['void_vs_invL2_r2'] = rv.rvalue ** 2
    got['void_smallest'] = float(void[np.isclose(x, 1.0 / min(Ls) ** 2)].mean())
    got['void_largest'] = float(void[np.isclose(x, 1.0 / max(Ls) ** 2)].mean())
    got['void_target'] = float(ch[Ls[0]][0]['VoF_void_sphere'])

    # ---- the three-covariate model, both ways ---------------------------
    for tag, y in (('corr', ycorr), ('raw', yraw)):
        beta, se, p = ols(y, x, brine, void)
        got['cov_%s_invL2' % tag] = beta[1]
        got['cov_%s_invL2_se' % tag] = se[1]
        got['cov_%s_invL2_p' % tag] = p[1]
        if tag == 'corr':
            got['cov_corr_brine_p'] = p[2]
            got['cov_corr_void_p'] = p[3]

    # ---- report ----------------------------------------------------------
    if args.verbose:
        print('%5s %6s %9s %9s %9s %9s %9s'
              % ('L', 'L/d', 'n', 'E_raw', 'baseline', 'E_corr', 'void'))
        for L in Ls:
            v = np.mean([float(r['porosity']) for r in ch[L]])
            m = np.mean([float(r['E_bending']) / 1e9 for r in ch[L]])
            print('%5.2f %6.0f %9d %9.3f %9.4f %9.3f %9.4f'
                  % (L, L / D_INCL, len(ch[L]), m, rel[L], m / rel[L], v))
        print()

    print('%-24s %12s %12s   %s' % ('quantity', 'computed', 'in the text', ''))
    bad = []
    for k, (want, tol) in QUOTED.items():
        have = got[k]
        ok = abs(have - want) <= tol
        if not ok:
            bad.append((k, have, want))
        print('%-24s %12.4f %12.4f   %s'
              % (k, have, want, 'ok' if ok else '<-- DIFFERS'))

    if args.power:
        # A null is only as strong as the effect it could have caught. MCST
        # puts a couple-stress length l on the slope as E_app = E_inf +
        # 12*mu*l^2 / L^2, so injecting a synthetic l into the observed
        # (corrected) data and re-running the same covariate model says which
        # length scales this sweep can actually exclude, rather than merely
        # failing to find.
        mu = np.mean([float(r['E_eff']) / 1e9 / (2 * (1 + float(r['nu_eff'])))
                      for L in Ls for r in ch[L]
                      if r.get('E_eff') not in (None, '', 'MISSING')])
        Lx = np.array([L for L in Ls for _ in ch[L]])
        beta, se, p = ols(ycorr, x, brine, void)
        dof = len(ycorr) - 4
        tcrit = stats.t.ppf(0.975, dof)
        upper = beta[1] + tcrit * se[1]              # 95% upper confidence limit
        l_max = np.sqrt(max(upper, 0.0) / (12.0 * mu))
        print()
        print('What the null excludes')
        print('  MCST puts the length scale on the slope as 12*mu*l^2, with the')
        print('  effective shear modulus mu = %.2f GPa here, so slope = %.1f*l^2 GPa.'
              % (mu, 12 * mu))
        print('  coefficient on 1/L^2  = %+.4f +- %.4f GPa' % (beta[1], se[1]))
        print('  95%% upper conf. limit = %+.4f GPa (t=%.3f, dof=%d)'
              % (upper, tcrit, dof))
        print('  => l < %.4f model units = %.2f d  (d = mean inclusion diameter)'
              % (l_max, l_max / D_INCL))
        print('  This is the bound the data support: not that l is zero, but')
        print('  that it is below a fraction of one inclusion diameter, which is')
        print('  smaller than any microstructural feature the cell resolves.')
        print()
        print('  Power check -- inject a synthetic l and re-run the same model:')
        print('  %8s %8s %12s %10s' % ('l/d', 'l (units)', 'coef 1/L^2', 'p'))
        detected = None
        for lod in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0):
            l = lod * D_INCL
            y = ycorr + 12.0 * mu * l ** 2 / Lx ** 2
            b2, s2, p2 = ols(y, x, brine, void)
            print('  %8.2f %8.4f %12.4f %10.4f' % (lod, l, b2[1], p2[1]))
            if detected is None and p2[1] < 0.05 and b2[1] > 0:
                detected = lod
        if detected is not None:
            print('  An l of %.1f d ON TOP of what is present would have turned the'
                  % detected)
            print('  coefficient significant. That is a detection threshold, not an')
            print('  exclusion bound -- the bound above is the one to quote.')

    print()
    print('The contrast the section turns on:')
    print('  coefficient on 1/L^2, baseline-corrected : %+.3f +- %.3f GPa (p=%.2f)'
          % (got['cov_corr_invL2'], got['cov_corr_invL2_se'], got['cov_corr_invL2_p']))
    print('  coefficient on 1/L^2, uncorrected        : %+.3f +- %.3f GPa (p=%.2f)'
          % (got['cov_raw_invL2'], got['cov_raw_invL2_se'], got['cov_raw_invL2_p']))
    print('  Both are null. A couple-stress length scale needs a positive,')
    print('  significant coefficient here, and neither model provides one.')

    if bad:
        print()
        print('%d quantity(ies) no longer match the manuscript:' % len(bad))
        for k, have, want in bad:
            print('   %-24s computed %.4f, text says %.4f' % (k, have, want))
        return 1
    print()
    print('All %d quoted quantities reproduce.' % len(QUOTED))
    return 0


if __name__ == '__main__':
    sys.exit(main())
