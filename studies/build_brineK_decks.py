"""Study #5 (redo) deck builder -- single-mesh, geometry-controlled brine K(T).

Mesh generation is non-deterministic (inclusions are randomly placed), so two
separate generation runs do NOT share geometry. To isolate the pure K(T) effect we
therefore generate ONE mesh per slice (the fixed-K=2.2 GPa "Kconst" decks, utx+utz)
and stamp the K(T) "Ktemp" twin onto that SAME mesh by rewriting only the brine
*Elastic card -- the second *Elastic block, the near-incompressible inclusion. The
paired decks then share byte-identical geometry, mesh, PBC and matrix; only the
brine bulk modulus differs, exactly as intended.

Brine card follows the solver's own convention (SpaX_Standalone.py:2197):
    E  = 9 K G / (3K + G),   nu = (3K - 2G) / (2(3K + G))
written with plain str() formatting.

In:  out_bkbase/Job-BKC_z{05..95}-{utx,utz}.inp   (fixed K = 2.2 GPa)
Out: out_brineK/ with, per slice, the Kconst deck (copied) and the Ktemp twin.
Usage: python3 build_brineK_decks.py
"""
import os, csv, shutil
from make_ice_studies import G_BRINE

SRC, OUT = 'out_bkbase', 'out_brineK'

def brine_card(K, G):
    E = 9.0 * K * G / (3.0 * K + G)
    nu = (3.0 * K - 2.0 * G) / (2.0 * (3.0 * K + G))
    return '{}, {}\n'.format(E, nu)

def stamp_ktemp(base_path, out_path, new_card):
    """Copy base deck, replacing the brine card (line after the 2nd *Elastic)."""
    with open(base_path) as f:
        lines = f.readlines()
    elastic_idx = [i for i, ln in enumerate(lines) if ln.strip() == '*Elastic']
    if len(elastic_idx) < 2:
        raise RuntimeError('expected 2 *Elastic blocks in %s, found %d'
                           % (base_path, len(elastic_idx)))
    brine_line = elastic_idx[1] + 1          # data line of the 2nd (brine) *Elastic
    lines[brine_line] = new_card
    with open(out_path, 'w') as f:
        f.writelines(lines)

def main():
    if os.path.isdir(OUT):
        shutil.rmtree(OUT)
    os.makedirs(OUT)
    ktemp = {r['run_id'].split('_z')[1]: float(r['K_inclusion'])
             for r in csv.DictReader(open('rve_brineKtemp.csv'))}
    n = 0
    for z, K in sorted(ktemp.items()):
        card = brine_card(K, G_BRINE)
        for mode in ('utx', 'utz'):
            base = os.path.join(SRC, 'Job-BKC_z%s-%s.inp' % (z, mode))
            if not os.path.exists(base):
                raise SystemExit('missing base mesh: %s' % base)
            shutil.copy(base, os.path.join(OUT, 'Job-BKC_z%s-%s.inp' % (z, mode)))
            stamp_ktemp(base, os.path.join(OUT, 'Job-BKT_z%s-%s.inp' % (z, mode)), card)
            n += 2
        print('z%s: K(T)=%.4f GPa  brine card=%s' % (z, K / 1e9, card.strip()))
    print('built %d decks in %s/ (expect 40 = 10 slices x 2 modes x 2 K-cases)' % (n, OUT))

    # verify each Ktemp deck differs from its Kconst base by exactly ONE line
    bad = 0
    for z in sorted(ktemp):
        for mode in ('utx', 'utz'):
            a = open(os.path.join(OUT, 'Job-BKC_z%s-%s.inp' % (z, mode))).readlines()
            b = open(os.path.join(OUT, 'Job-BKT_z%s-%s.inp' % (z, mode))).readlines()
            d = sum(1 for x, y in zip(a, b) if x != y)
            if d != 1 or len(a) != len(b):
                bad += 1
                print('  MISMATCH z%s %s: %d differing lines (len %d/%d)'
                      % (z, mode, d, len(a), len(b)))
    print('verification: %s' % ('all pairs differ by exactly the brine line OK'
                                 if bad == 0 else '%d BAD pairs' % bad))

if __name__ == '__main__':
    main()
