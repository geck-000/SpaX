"""Did every RVE the deck asks for actually survive to the results table?

The volume audit that gates the re-run checks whether the cells that *were*
built carry the right microstructure. It says nothing about cells that were
never built at all. A cell can be lost quietly: generation exhausts its mesh
retries and skips the row, the solve array skips a deck whose ODB already
exists, or a post-processor writes a table with fewer rows than the deck. The
chain advances in every case, and the results file exists and looks current.

This walks the manifest and, for each campaign, compares four counts that
should agree:

    deck rows  ->  decks generated  ->  ODBs solved  ->  rows in the results

and names the run_ids lost at each step.

    python3 check_completeness.py [manifest] [workdir] [--since EPOCH]

`--since` is what makes this trustworthy. A results table left over from an
earlier campaign has exactly the run_ids the deck asks for, so on row counts
alone it looks complete; only its date betrays it. Pass the time the re-run
started and any table older than that is reported STALE rather than ok.

Exits non-zero if anything is missing or stale, so it can gate a downstream
step.
"""
import csv
import glob
import os
import re
import sys


_SOLVE_RC = {}


def solve_outcomes(work):
    """deck -> exit code of its most recent solve, read from the solve logs.

    An Abaqus killed for memory or walltime still leaves an ODB behind, and it
    is indistinguishable from a good one by existence or by size -- a partly
    written 60 MB file looks exactly like a small complete one. The solve array
    then skips that deck on any resubmit, because its ODB is present, so the
    loss is silent and the post-processor reads whatever is in the file.
    Counting ODBs cannot catch this.

    The solver prunes everything but the deck and the ODB, so the .sta is gone
    and the only durable record is the array task's own log, which names the
    deck and ends with 'Abaqus exit: <rc>'. A task killed by Slurm never
    reaches that line, so a missing exit code is exactly the signature wanted.
    """
    if work in _SOLVE_RC:
        return _SOLVE_RC[work]
    out = {}
    logs = glob.glob(os.path.join(work, 'logs', 'csc_solve_*.out'))
    for path in sorted(logs, key=lambda p: os.path.getmtime(p)):
        deck, rc = None, None
        try:
            with open(path, encoding='utf8', errors='replace') as fh:
                for line in fh:
                    if line.startswith('====='):
                        m = re.match(r'=====\s+(\S+)\s', line)
                        if m and deck is None:
                            deck = m.group(1)
                    elif line.startswith('Abaqus exit:'):
                        try:
                            rc = int(line.split(':', 1)[1].strip())
                        except ValueError:
                            rc = None
        except OSError:
            continue
        if deck:
            out[deck] = rc            # later log wins: mtime-sorted
    _SOLVE_RC[work] = out
    return out


def truncated_odbs(work, ids):
    """ODBs present whose most recent solve did not exit cleanly."""
    rc = solve_outcomes(work)
    bad = []
    for rid in ids:
        for odb in glob.glob(os.path.join(work, 'Job-%s-*.odb' % rid)):
            deck = os.path.basename(odb)[:-4]
            if deck in rc and rc[deck] != 0:
                bad.append('%s (exit %s)' % (deck, rc[deck]))
    return bad


def deck_ids(path):
    if not os.path.isfile(path):
        return None
    with open(path, encoding='utf8', errors='replace') as fh:
        return [r['run_id'] for r in csv.DictReader(fh) if r.get('run_id')]


def results_ids(path):
    if not os.path.isfile(path):
        return None
    out = []
    with open(path, encoding='utf8', errors='replace') as fh:
        for r in csv.DictReader(fh):
            rid = r.get('run_id')
            if not rid:
                continue
            # results carry one row per RVE; tolerate a -mode suffix
            out.append(rid.rsplit('-', 1)[0] if rid.rsplit('-', 1)[-1] in
                       ('utx', 'uty', 'utz', 'ss12', 'ss13', 'ss23', 'ben')
                       else rid)
    return out


def main():
    argv = [a for a in sys.argv[1:]]
    since = None
    if '--since' in argv:
        i = argv.index('--since')
        since = float(argv[i + 1])
        del argv[i:i + 2]
    man = argv[0] if len(argv) > 0 else 'rerun_paper_manifest.tsv'
    work = argv[1] if len(argv) > 1 else '.'

    rows = []
    with open(man, encoding='utf8', errors='replace') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            f = line.split('\t')
            if len(f) >= 6:
                rows.append(f)

    print('%-20s %6s %6s %6s %6s   %s'
          % ('campaign', 'deck', 'gen', 'odb', 'result', 'status'))
    print('-' * 86)

    bad = []
    for f in rows:
        name, deck, outdir, _glob, results = f[0], f[1], f[2], f[3], f[4]
        ids = deck_ids(os.path.join(work, 'params', deck))
        if ids is None:
            print('%-20s %6s %6s %6s %6s   deck missing' % (name, '-', '-', '-', '-'))
            bad.append((name, 'deck missing', []))
            continue

        gen, odb = [], []
        for rid in ids:
            g = (glob.glob(os.path.join(work, outdir, '**', 'Job-%s-*.inp' % rid),
                           recursive=True)
                 or glob.glob(os.path.join(work, 'Job-%s-*.inp' % rid)))
            if g:
                gen.append(rid)
            if glob.glob(os.path.join(work, 'Job-%s-*.odb' % rid)):
                odb.append(rid)

        rpath = os.path.join(work, results)
        res = results_ids(rpath)
        nres = len(set(res) & set(ids)) if res is not None else 0
        stale = (since is not None and os.path.isfile(rpath)
                 and os.path.getmtime(rpath) < since)

        trunc = truncated_odbs(work, ids)

        status = 'ok'
        if res is None:
            status = 'NO RESULTS FILE'
        elif stale:
            status = 'STALE (predates this run)'
        elif nres < len(ids):
            status = 'INCOMPLETE'
        elif trunc:
            status = 'TRUNCATED ODBs (%d)' % len(trunc)
        elif len(odb) < len(ids):
            status = 'results ok, some ODBs absent'

        print('%-20s %6d %6d %6d %6d   %s'
              % (name, len(ids), len(gen), len(odb), nres, status))

        if status != 'ok' and not status.startswith('results ok'):
            if status.startswith('TRUNCATED'):
                bad.append((name, status, trunc))
            else:
                missing = [r for r in ids if res is None or r not in set(res)]
                bad.append((name, status, missing))

    print()
    if not bad:
        print('All campaigns complete: every deck row reached its results table.')
        return 0

    print('INCOMPLETE CAMPAIGNS')
    for name, status, missing in bad:
        print('  %-20s %s' % (name, status))
        for m in missing[:12]:
            ng = 'no deck' if not glob.glob(
                os.path.join(work, 'Job-%s-*.inp' % m)) else 'deck built'
            print('      %-24s (%s)' % (m, ng))
        if len(missing) > 12:
            print('      ... and %d more' % (len(missing) - 12))
    return 1


if __name__ == '__main__':
    sys.exit(main())
