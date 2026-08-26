"""What did the mesher fix change, measured against replicate scatter?

Every published number moves when inclusions stop being meshed as bounding
spheres, because the cells no longer carry 1.6-2.2x the intended soft phase.
The interesting question is not whether a number moved but whether it moved by
more than the packing-to-packing scatter, which is the criterion the papers use
throughout for calling an effect real.

Replicates are recognised by a trailing _sN on the run_id, so a column of ten
slices at ten packings groups into ten ensembles. Scatter is the population
standard deviation (ddof=0), as everywhere else in this project.

    python3 compare_rerun.py <old.csv> <new.csv> [label]

Reports, per group, the ensemble mean before and after, the shift, and the
shift expressed in units of the pooled replicate scatter. A shift below about
one sigma is not resolvable at this ensemble size, however large it looks as a
percentage.
"""
import collections
import csv
import re
import statistics as st
import sys


def load(path):
    """group -> {quantity: [values over replicates]}"""
    g = collections.defaultdict(lambda: collections.defaultdict(list))
    with open(path, encoding='utf8', errors='replace') as fh:
        for r in csv.DictReader(fh):
            rid = (r.get('run_id') or '').strip()
            if not rid:
                continue
            key = re.sub(r'_s\d+$', '', rid)

            def num(k):
                try:
                    v = float(r.get(k))
                    return v if v > 0 else None
                except (TypeError, ValueError):
                    return None

            ex, ez = num('E_x'), num('E_z')
            eeff = num('E_eff') or ex
            if eeff:
                g[key]['E'].append(eeff / 1e9)
            if ex and ez:
                g[key]['aniso'].append(ez / ex)
            phi = num('phi_inclusion')
            if phi:
                g[key]['phi'].append(phi)
    return g


def stats(v):
    if not v:
        return None, None, 0
    return st.mean(v), st.pstdev(v), len(v)


def main():
    old, new = load(sys.argv[1]), load(sys.argv[2])
    label = sys.argv[3] if len(sys.argv) > 3 else 'campaign'
    keys = [k for k in sorted(new) if k in old]
    if not keys:
        print('no groups in common')
        return 1

    print('=' * 92)
    print('%s : old (bounding spheres) vs new (true ellipsoids)' % label)
    print('=' * 92)

    for q, name, unit in (('E', 'effective modulus', 'GPa'),
                          ('aniso', 'anisotropy E_z/E_x', ''),
                          ('phi', 'realised inclusion fraction', '')):
        rows = [(k, stats(old[k][q]), stats(new[k][q])) for k in keys]
        rows = [r for r in rows if r[1][0] is not None and r[2][0] is not None]
        if not rows:
            continue
        print()
        print('%s %s' % (name, ('[%s]' % unit) if unit else ''))
        print('%-16s %9s %9s %9s %9s %8s %8s' %
              ('group', 'old', 'old sd', 'new', 'new sd', 'change', 'sigma'))
        shifts = []
        for k, (mo, so, no), (mn, sn, nn) in rows:
            pooled = ((so or 0) ** 2 + (sn or 0) ** 2) ** 0.5 or float('nan')
            sig = (mn - mo) / pooled if pooled == pooled and pooled else float('nan')
            pct = 100 * (mn / mo - 1) if mo else float('nan')
            shifts.append((pct, sig))
            print('%-16s %9.4f %9.4f %9.4f %9.4f %+7.1f%% %+8.1f'
                  % (k, mo, so or 0, mn, sn or 0, pct, sig))
        mp = st.mean([s[0] for s in shifts])
        ms = st.mean([abs(s[1]) for s in shifts if s[1] == s[1]])
        resolved = sum(1 for s in shifts if s[1] == s[1] and abs(s[1]) > 1.0)
        print('  mean change %+.1f%% ; mean |shift| %.1f sigma ; '
              '%d of %d groups move by more than one sigma'
              % (mp, ms, resolved, len(shifts)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
