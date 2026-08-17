"""One line per pending campaign: how many cells have actually produced a number.

Existence and row count both lie. A postprocessor enumerates the deck and emits
a row per cell whether or not that cell solved, so a file can be freshly written,
carry the full complement of rows, and still be almost empty -- which is exactly
what results_eringen_layer.csv did. The only honest measure is how many rows
carry a real value in the column the campaign exists to produce.
"""
import csv
import os
import subprocess
import sys

WORK = '/scratch/project_2019020/test_rve'

# file, expected cells, the column that must be populated
CAMPAIGNS = [
    ('results_eringen_layer.csv',            18, 'E_bending'),
    ('results_eringen_layer_homog.csv',       5, 'E_bending'),
    ('results_torsion_big_K.csv',             7, 'K_rve'),
    ('results_weibull_layer_scf.csv',        20, None),   # any SCF-like column
]


def quota():
    """Percent of the scratch allocation in use.

    Reported on every poll because the failure it precedes is silent: at the
    limit Abaqus is killed mid-write and leaves a truncated ODB, so the campaign
    looks solved and the extraction dies later with 'database file is corrupt'.
    The cost is a re-solve, and it lands on whatever happened to be writing --
    not necessarily the campaign that filled the disk.
    """
    try:
        out = subprocess.check_output(
            'lfs quota -q -p 602019020 /scratch 2>/dev/null | head -1',
            shell=True).decode().split()
        used = int(out[1].replace('*', ''))
        lim = int(out[3])
        return 100.0 * used / lim, used / 1024.0 / 1024.0, lim / 1024.0 / 1024.0
    except Exception:
        return None, None, None


def populated(path, col):
    if not os.path.exists(path):
        return 0
    try:
        rows = list(csv.DictReader(open(path)))
    except Exception:
        return 0
    if not rows:
        return 0
    if col is None:
        cand = [c for c in rows[0] if 'scf' in c.lower() or 'sigma_max' in c.lower()]
        if not cand:
            return 0
        col = cand[0]
    n = 0
    for r in rows:
        v = (r.get(col) or '').strip()
        if v in ('', 'MISSING', 'nan'):
            continue
        try:
            if float(v) != 0.0:
                n += 1
        except ValueError:
            pass
    return n


def main():
    os.chdir(WORK)
    done = True
    for f, exp, col in CAMPAIGNS:
        n = populated(f, col)
        state = 'complete' if n >= exp else 'waiting'
        if n < exp:
            done = False
        print('%-38s %2d/%-2d %s' % (f, n, exp, state))
    try:
        q = subprocess.check_output(
            'squeue -u $USER -h | wc -l', shell=True).decode().strip()
    except Exception:
        q = '?'
    print('QUEUE %s' % q)
    pct, used, lim = quota()
    if pct is not None:
        flag = 'CRITICAL' if pct >= 90 else ('HIGH' if pct >= 75 else 'ok')
        print('QUOTA %.0f%% (%.0fG of %.0fG) %s' % (pct, used, lim, flag))
    print('ALLDONE' if done else 'PENDING')


if __name__ == '__main__':
    main()
