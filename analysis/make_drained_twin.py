# -*- coding: utf-8 -*-
"""Drained twin of an existing deck, on the identical mesh.

Drainage changes exactly one thing in a deck: the inclusion's bulk modulus.
Geometry, packing, mesh and every periodic constraint are untouched. So the
twin does not need generating -- rewriting the one elastic card that carries
the brine turns an undrained deck into its drained counterpart, on the same
nodes and the same elements.

That matters for more than cost. Regenerating would re-pack, and the pair would
then differ by realisation scatter as well as by drainage, which at the few
percent the pocket bracket spans would swamp the effect being measured. Sharing
the mesh removes that entirely: the difference between the two solves is
drainage and nothing else. It is the same paired design the brine-K(T) study
used, taken one step further by not regenerating at all.

    python3 analysis/make_drained_twin.py <dir> [K_drained_Pa]

Writes Job-<id>-drn-<mode>.inp beside each Job-<id>-<mode>.inp it finds.
"""
import glob
import os
import re
import sys

G_BRINE = 440029.33528897085     # Pa, unchanged by drainage
K_DRAINED = 2.2e6                # Pa; brine at 2.2 GPa is the undrained end


def en_from_kg(K, G):
    """Young's modulus and Poisson ratio from bulk and shear."""
    E = 9.0 * K * G / (3.0 * K + G)
    nu = (3.0 * K - 2.0 * G) / (2.0 * (3.0 * K + G))
    return E, nu


def twin(path, K=K_DRAINED, G=G_BRINE):
    """Rewrite the inclusion elastic card; everything else is copied verbatim."""
    with open(path) as f:
        lines = f.readlines()

    E, nu = en_from_kg(K, G)
    out, i, done = [], 0, False
    while i < len(lines):
        out.append(lines[i])
        if re.match(r'\*Material\s*,\s*name\s*=\s*Mat_Inclusion', lines[i], re.I):
            # the card is *Material / *Elastic / E, nu -- replace the third line
            if i + 2 < len(lines) and lines[i + 1].strip().lower().startswith('*elastic'):
                out.append(lines[i + 1])
                out.append('%.10g, %.6f\n' % (E, nu))
                i += 3
                done = True
                continue
        i += 1

    if not done:
        raise ValueError('no Mat_Inclusion elastic card found in %s' % path)

    base = os.path.basename(path)
    m = re.match(r'(Job-.+)-(ut[xyz]|ss\d\d|ben|tor)\.inp$', base)
    if not m:
        raise ValueError('unexpected deck name %s' % base)
    new = os.path.join(os.path.dirname(path),
                       '%s-drn-%s.inp' % (m.group(1), m.group(2)))
    with open(new, 'w') as f:
        f.writelines(out)
    return new, E, nu


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else '.'
    K = float(sys.argv[2]) if len(sys.argv) > 2 else K_DRAINED
    found = [f for f in sorted(glob.glob(os.path.join(d, 'Job-*.inp')))
             if '-drn-' not in os.path.basename(f)]
    if not found:
        print('no decks in %s' % d)
        return 1
    E, nu = en_from_kg(K, G_BRINE)
    print('drained fill: K = %.4g Pa, G = %.4g Pa -> E = %.6g Pa, nu = %.6f'
          % (K, G_BRINE, E, nu))
    n = 0
    for f in found:
        try:
            new, _, _ = twin(f, K)
        except ValueError as e:
            print('  skip: %s' % e)
            continue
        n += 1
    print('wrote %d twin deck(s) in %s' % (n, d))
    return 0


if __name__ == '__main__':
    sys.exit(main())
