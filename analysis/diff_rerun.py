# -*- coding: utf-8 -*-
"""What did the geometry fixes actually change, campaign by campaign?

The generator was wrong in several ways at once -- inclusions built as their
bounding spheres, semi-axes that did not preserve volume, an orientation
sampler that returned the same aligned pack whatever concentration was asked
for, and two phase fractions that overshot their targets. Every number in the
paper was computed on cells carrying at least one of those. The re-run fixes
them, but a re-run only helps if we know WHICH conclusions moved: rewriting the
text against 1400 new cells by reading them is not feasible and not reliable.

This compares an old results file against its re-run counterpart, per run_id
and per column, and reports the shift against the replicate scatter of the
campaign itself. A change smaller than the scatter it sits in is not a change,
and the paper's own convention is a population standard deviation (ddof = 0),
so that is what is used here.

    python3 analysis/diff_rerun.py OLD.csv NEW.csv [--quiet]

Columns are matched by name, rows by run_id. Anything that fails to parse as a
number on either side is reported as unusable rather than silently skipped,
because a column that turned into ERROR is itself a finding.
"""
import csv
import os
import sys

import numpy as np

# quantities the paper actually quotes; others are echoed but not ranked
KEY = ('E_x', 'E_y', 'E_z', 'E_eff', 'G_eff', 'nu_eff',
       'E_anisotropy', 'E_z_over_xy', 'phi_soft_total', 'phi_inclusion',
       'porosity', 'K_rve')


def load(path):
    rows = {}
    with open(path, encoding='utf8', errors='replace') as f:
        for r in csv.DictReader(f):
            rid = (r.get('run_id') or '').strip()
            if rid:
                rows[rid] = r
    return rows


def num(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return x if np.isfinite(x) else None


def base(rid):
    """Strip a trailing _sN replicate tag so scatter can be pooled."""
    return rid.rsplit('_s', 1)[0] if '_s' in rid else rid


def scatter(rows, col):
    """Population s.d. within replicate groups, pooled over the campaign."""
    g = {}
    for rid, r in rows.items():
        v = num(r.get(col))
        if v is not None:
            g.setdefault(base(rid), []).append(v)
    sds = [float(np.std(v)) for v in g.values() if len(v) > 1]
    return float(np.mean(sds)) if sds else None


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    old_p, new_p = sys.argv[1], sys.argv[2]
    quiet = '--quiet' in sys.argv
    for p in (old_p, new_p):
        if not os.path.exists(p):
            print('missing: %s' % p)
            return 2
    old, new = load(old_p), load(new_p)

    shared = sorted(set(old) & set(new))
    print('%s  ->  %s' % (os.path.basename(old_p), os.path.basename(new_p)))
    print('rows: %d old, %d new, %d shared' % (len(old), len(new), len(shared)))
    gone = sorted(set(old) - set(new))
    added = sorted(set(new) - set(old))
    if gone:
        print('  only in old (%d): %s%s' % (
            len(gone), ', '.join(gone[:4]), ' ...' if len(gone) > 4 else ''))
    if added:
        print('  only in new (%d): %s%s' % (
            len(added), ', '.join(added[:4]), ' ...' if len(added) > 4 else ''))
    if not shared:
        print('  no shared run_ids -- cannot compare')
        return 1

    cols = [c for c in (old[shared[0]].keys()) if c in new[shared[0]]]
    ranked = []
    for c in cols:
        if c == 'run_id':
            continue
        pairs, bad = [], 0
        for rid in shared:
            a, b = num(old[rid].get(c)), num(new[rid].get(c))
            if a is None or b is None:
                if (old[rid].get(c) or '').strip() or (new[rid].get(c) or '').strip():
                    bad += 1
                continue
            pairs.append((a, b))
        if not pairs:
            continue
        a = np.array([p[0] for p in pairs])
        b = np.array([p[1] for p in pairs])
        with np.errstate(divide='ignore', invalid='ignore'):
            rel = np.where(a != 0, (b - a) / np.abs(a), np.nan)
        mrel = float(np.nanmean(rel)) if np.isfinite(rel).any() else float('nan')
        sd = scatter(new, c)
        # shift measured in units of the campaign's own replicate scatter
        nsig = (abs(float(np.mean(b - a))) / sd) if sd else float('nan')
        ranked.append((c, float(np.mean(a)), float(np.mean(b)),
                       100 * mrel, nsig, bad, len(pairs)))

    ranked.sort(key=lambda r: (-(r[4] if np.isfinite(r[4]) else 0),
                               -abs(r[3] if np.isfinite(r[3]) else 0)))
    print('\n%-18s %12s %12s %9s %9s %6s' % (
        'column', 'old mean', 'new mean', 'change %', 'n_sigma', 'bad'))
    for c, am, bm, pc, ns, bad, n in ranked:
        if quiet and c not in KEY:
            continue
        flag = ''
        if np.isfinite(ns):
            flag = ' <<<' if ns > 3 else ('  <<' if ns > 1 else '')
        print('%-18s %12.5g %12.5g %+9.2f %9s %6d%s' % (
            c, am, bm, pc, ('%.1f' % ns) if np.isfinite(ns) else '-', bad, flag))

    big = [r for r in ranked if r[0] in KEY and np.isfinite(r[4]) and r[4] > 3]
    print('\n%d quoted quantit%s moved by more than 3x the replicate scatter%s'
          % (len(big), 'y' if len(big) == 1 else 'ies',
             ':' if big else '.'))
    for r in big:
        print('   %-16s %.5g -> %.5g  (%+.1f%%, %.1f sigma)'
              % (r[0], r[1], r[2], r[3], r[4]))
    print('\nA shift under one sigma is not a change. Anything marked <<< needs')
    print('the sentence that quotes it rewritten, not just the number swapped.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
