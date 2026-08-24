#!/usr/bin/env python3
"""fbares.py -- rewrite a C3D4 deck as F-barES-FEM-T4.

Onishi, Iida & Amaya, Int. J. Comput. Methods 15(7) 1845003 (2018).  The
formulation, its small-strain reduction and the verification ladder are in
elements_ccx/docs/fbar_es_fem_t4.md; this script is only the deck side.

    K = K_dev   per EDGE (U2): V_h Bt^T D_dev Bt          eq. (1), (4), (13)
      + K_vol   per EDGE (U3): (K V_h) tbar^T sbar        eq. (6)-(11), (17)

The base tets stay in the deck retyped to U5 and are run with CCX_U5_ZERO=1,
contributing nothing: U2 and U3 carry the whole stiffness and both read the
tets for geometry through the 'U5' node->element map.  That is a different
arrangement from U5+U6, where U5 carries the deviatoric itself.

WHY THE CONNECTIVITY IS WRITTEN OUT AT ALL.  The elements recompute their own
weights from geometry -- which subsets of a ring form tets is not recoverable
from a node list -- but ccx has to know the stencil to size the element and to
build the matrix structure, so the node SET must be in the deck and must match
what u3vol walks exactly.  u3vol stops and names the element if it meets a
node the connectivity does not carry, so a mismatch is loud.

STENCIL WIDTH IS THE BINDING CONSTRAINT.  Measured on LMESH_m0p0240's soft
phase (117437 tets, 36323 nodes, 184572 edges):

    U2 deviatoric      mean  6.3 nodes/edge   max  14   ->   42 DOF
    U3 volumetric c=1  mean 33.7              max 173   ->  519 DOF
    U3 volumetric c=2  mean 93.6              max 494   -> 1482 DOF

c=2 cannot be expressed at all: ccx carries a user element's node count in the
single byte lakon(8:8) and userelements.f rejects NODES > 255.  This script
refuses it rather than letting the solver truncate.  c=1 needs the element
matrix widened from 150 to 520 DOF along the whole e_c3d_u* path.
"""
import argparse
import os
import sys
from collections import defaultdict

import numpy as np
import scipy.sparse as sp

LETTERS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
# Two suffix characters, LETTERS ONLY.  ccx's built-in element dispatch keys on
# digits in the label (elements.f: label(4:4).eq.'4' -> nope=4, '10' -> 10,
# '20' -> 20), so a type named U214 is claimed by the nope=4 rule before the
# *USER ELEMENT lookup and ccx reads 4 nodes instead of 14.  nodalbbar.py hit
# exactly this and it showed up on only 3 of ~36000 elements.
NAMES = [a + b for a in LETTERS for b in LETTERS]


def iskw(line):
    s = line.lstrip()
    return s.startswith('*') and not s.startswith('**')


def parse(path):
    """nodes, tets, elset membership, elset -> material, material -> (E, nu)."""
    nodes, tets = {}, {}
    elset_of, mat_of, elastic = defaultdict(list), {}, {}
    mode, cur, curmat = None, None, None
    for ln in open(path):
        s = ln.strip()
        if not s or s.startswith('**'):
            continue
        if iskw(ln):
            u = s.upper().replace(' ', '')
            mode = None
            if u.startswith('*NODE') and 'OUTPUT' not in u:
                mode = 'n'
            elif u.startswith('*ELEMENT'):
                if 'TYPE=C3D4' in u:
                    mode = 'e'
                    cur = next((p.split('=')[1] for p in s.split(',')
                                if p.strip().upper().startswith('ELSET=')),
                               'ALL').strip()
            elif u.startswith('*SOLIDSECTION'):
                es = next((p.split('=')[1] for p in s.split(',')
                           if p.strip().upper().startswith('ELSET=')), None)
                mt = next((p.split('=')[1] for p in s.split(',')
                           if p.strip().upper().startswith('MATERIAL=')), None)
                if es and mt:
                    mat_of[es.strip()] = mt.strip()
            elif u.startswith('*MATERIAL'):
                curmat = next((p.split('=')[1] for p in s.split(',')
                               if p.strip().upper().startswith('NAME=')),
                              None)
                if curmat:
                    curmat = curmat.strip()
            elif u.startswith('*ELASTIC'):
                mode = 'el'
            continue
        f = [x.strip() for x in s.split(',') if x.strip()]
        if mode == 'n' and len(f) >= 4:
            nodes[int(f[0])] = (float(f[1]), float(f[2]), float(f[3]))
        elif mode == 'e' and len(f) >= 5:
            tets[int(f[0])] = tuple(int(x) for x in f[1:5])
            elset_of[cur].append(int(f[0]))
        elif mode == 'el' and curmat and len(f) >= 2:
            elastic.setdefault(curmat, (float(f[0]), float(f[1])))
            mode = None
    return nodes, tets, elset_of, mat_of, elastic


def stencils(conn, ncyc, chunk=20000):
    """Edge list and, per edge, the U2 ring and the U3 = E A^c node support.

    Built with boolean sparse products rather than a Python walk: on the real
    cells this is ~500k edges and the walk is the whole runtime.  The PATTERN
    is all that is needed -- the elements recompute the weights themselves --
    and verify_u8_chain.py has already checked that the weighted walk agrees
    with the prototype operator to 1e-15.

    conn is (ne, 4) LOCAL node indices for one material's tets.
    """
    ne = len(conn)
    nn = int(conn.max()) + 1 if ne else 0
    inc = sp.coo_matrix((np.ones(4 * ne, dtype=np.int8),
                         (np.repeat(np.arange(ne), 4), conn.ravel())),
                        shape=(ne, nn)).tocsr()
    inc.data[:] = 1

    # edges of this material, and the elements at each
    ekey = {}
    for e in range(ne):
        t = conn[e]
        for i in range(4):
            for j in range(i + 1, 4):
                a, b = (t[i], t[j]) if t[i] < t[j] else (t[j], t[i])
                ekey.setdefault((int(a), int(b)), []).append(e)
    keys = sorted(ekey)
    rows = np.repeat(np.arange(len(keys)),
                     [len(ekey[k]) for k in keys])
    cols = np.fromiter((e for k in keys for e in ekey[k]), dtype=np.int64)
    E = sp.coo_matrix((np.ones(len(cols), dtype=np.int8), (rows, cols)),
                      shape=(len(keys), ne)).tocsr()
    E.data[:] = 1

    # A = P Q as a PATTERN: elements sharing a node
    A = (inc @ inc.T).tocsr()
    A.data[:] = 1

    ring, supp = [], []
    for lo in range(0, len(keys), chunk):
        blk = E[lo:lo + chunk]
        r = (blk @ inc).tocsr()          # U2: nodes of the edge's own tets
        r.data[:] = 1
        for i in range(r.shape[0]):
            ring.append(r.indices[r.indptr[i]:r.indptr[i + 1]])
        s = blk
        for _ in range(ncyc):
            s = (s @ A).tocsr()
            s.data[:] = 1
        w = (s @ inc).tocsr()            # U3: nodes of the E A^c support
        w.data[:] = 1
        for i in range(w.shape[0]):
            supp.append(w.indices[w.indptr[i]:w.indptr[i + 1]])
    return keys, ring, supp


def card(eid, conn_nodes):
    """One *ELEMENT card, wrapped.

    At most 15 fields per line and NO trailing comma.  textpart in elements.f
    is dimensioned (16); a trailing comma adds a field and overflows it, after
    which ccx reads the continuation as a fresh element and reports
    'element N is already defined'.  ccx continues on node count instead.
    """
    row = [str(eid)] + [str(x) for x in conn_nodes]
    return [','.join(row[c:c + 15]) for c in range(0, len(row), 15)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('src')
    ap.add_argument('dst')
    ap.add_argument('--elset', action='append', default=None,
                    help='element set to treat; repeat for more than one. '
                         'Default: every C3D4 set in the deck.')
    ap.add_argument('--cycles', type=int,
                    default=int(os.environ.get('CCX_FBAR_C', 1)),
                    help='c, the number of cyclic smoothings of J, eq. (6)-(7)')
    ap.add_argument('--solver', default=os.environ.get('SPAX_FBAR_SOLVER',
                                                       'PARDISO'))
    a = ap.parse_args()

    nodes, tets, elset_of, mat_of, elastic = parse(a.src)
    sets = a.elset or sorted(elset_of)
    for es in sets:
        if es not in elset_of:
            raise SystemExit('fbares: no C3D4 elements in ELSET=%s '
                             '(deck has %s)' % (es, sorted(elset_of)))
        if es not in mat_of:
            raise SystemExit('fbares: ELSET=%s has no *SOLID SECTION' % es)

    if a.cycles < 0 or a.cycles > 3:
        raise SystemExit('fbares: --cycles must be 0..3')

    lines = open(a.src).read().splitlines()
    maxel = max(tets) if tets else 0

    # --- build the stencils, per material -------------------------------
    #
    # PER MATERIAL, always.  The smoothing of eqs. (1) and (6)-(8) must not
    # cross a phase boundary: averaging a divergence across a 1000x modulus
    # contrast is not the method, and u6patch.f records that a shared patch
    # took equilibrium_gap from 2.5e-3 to 2.6e-1 on the finer layered cell.
    # An interface edge legitimately gets one smoothing domain per phase.
    u2, u3 = {}, {}
    stats = {}
    for es in sets:
        ids = sorted(elset_of[es])
        loc = {}
        conn = np.empty((len(ids), 4), dtype=np.int64)
        for i, e in enumerate(ids):
            for j, n in enumerate(tets[e]):
                if n not in loc:
                    loc[n] = len(loc)
                conn[i, j] = loc[n]
        back = np.empty(len(loc), dtype=np.int64)
        for n, i in loc.items():
            back[i] = n
        keys, ring, supp = stencils(conn, a.cycles)
        u2[es] = [(int(back[k[0]]), int(back[k[1]]), back[r]) 
                  for k, r in zip(keys, ring)]
        u3[es] = [(int(back[k[0]]), int(back[k[1]]), back[s])
                  for k, s in zip(keys, supp)]
        rs = np.array([len(r) for r in ring])
        ss = np.array([len(s) for s in supp])
        stats[es] = (len(ids), len(loc), len(keys), rs, ss)

    # --- the 255-node wall ----------------------------------------------
    worst = max((s[4].max(), es) for es, s in stats.items())
    if worst[0] > 255:
        raise SystemExit(
            'fbares: the c=%d volumetric stencil reaches %d nodes in '
            'ELSET=%s. ccx stores a user element\'s node count in the single '
            'byte lakon(8:8) and userelements.f rejects NODES > 255, so this '
            'cannot be expressed as an element at all. Use --cycles 1, or '
            'assemble the volumetric term outside the element loop (see '
            'docs/fbar_es_fem_t4.md section 8).'
            % (a.cycles, worst[0], worst[1]))

    # --- type names: one *USER ELEMENT per (kind, elset, node count) -----
    groups = {}
    for kind, table in (('U2', u2), ('U3', u3)):
        for es in sets:
            for na, nb, nl in table[es]:
                # nl ALREADY contains both edge nodes -- it is the node set of
                # the tets at the edge (U2) or of the E A^c support (U3), and
                # both contain the edge itself.  Declaring len(nl)+2 made ccx
                # read two fields past the end of every card and swallow the
                # next element's id, which showed up as duplicate ids and a
                # connectivity carrying its own successor.
                groups.setdefault((kind, es, len(nl)), []).append(
                    (na, nb, nl))
    for kind in ('U2', 'U3'):
        g = [k for k in groups if k[0] == kind]
        if len(g) > len(NAMES):
            raise SystemExit('fbares: %d %s (elset, size) groups, more than '
                             'the %d type names available'
                             % (len(g), kind, len(NAMES)))
    suffix = {}
    for kind in ('U2', 'U3'):
        for i, k in enumerate(sorted(x for x in groups if x[0] == kind)):
            suffix[k] = NAMES[i]

    # --- emit -------------------------------------------------------------
    out, u5decl, done_step = [], False, False
    eid = maxel
    for ln in lines:
        if iskw(ln):
            u = ln.upper().replace(' ', '')
            if u.startswith('*ELEMENT') and any(
                    ('ELSET=' + e).upper() in u for e in sets):
                if not u5decl:
                    out.append('*USER ELEMENT,TYPE=U5,NODES=4,'
                               'INTEGRATIONPOINTS=1,MAXDOF=3')
                    u5decl = True
                out.append(ln.replace('C3D4', 'U5').replace('c3d4', 'U5'))
                continue
            if u.startswith('*STATIC'):
                # The volumetric operator of eq. (17) is NOT symmetric, so the
                # assembled matrix is not either and incomplete-Cholesky PCG
                # does not apply.  PARDISO takes it through the asymmetric
                # path (nasym -> mafillsmas -> mtype=11).
                ln = '*STATIC, SOLVER=' + a.solver
            if u.startswith('*STEP') and not done_step:
                done_step = True
                # ALL *USER ELEMENT declarations first, then all *ELEMENT
                # blocks.  ccx sorts the deck into per-keyword chains
                # (keystart.f: *USER ELEMENT is position 3, *ELEMENT is 4), so
                # interleaving them puts a multi-line element's continuation
                # at a chain boundary, where ccx reads it as a fresh element
                # and reports 'element N is already defined'.
                out.append('** SPAX F-barES-FEM-T4(c=%d): %d edge domains'
                           % (a.cycles, sum(len(u2[e]) for e in sets)))
                for k in sorted(groups):
                    out.append('*USER ELEMENT,TYPE=%s%s,NODES=%d,'
                               'INTEGRATIONPOINTS=1,MAXDOF=3'
                               % (k[0], suffix[k], k[2]))
                for k in sorted(groups):
                    kind, es, sz = k
                    tag = 'FBD' if kind == 'U2' else 'FBV'
                    out.append('*ELEMENT,TYPE=%s%s,ELSET=%s_%s'
                               % (kind, suffix[k], tag, es))
                    for na, nb, nl in groups[k]:
                        eid += 1
                        # konl(1), konl(2) ARE THE EDGE NODES, in both U2 and
                        # U3; u2edge/u3vol identify the edge from them and
                        # find the tets themselves.
                        rest = [int(x) for x in nl if x != na and x != nb]
                        out.extend(card(eid, [na, nb] + sorted(rest)))
                # Their OWN elset + *SOLID SECTION, never the phase's: an
                # element in the phase elset joins its *EL PRINT set, and
                # printoutelem.f would try to integrate a smoothing domain
                # that has no material volume of its own -- and it would
                # corrupt the volume-averaged stress the homogenisation reads.
                for k in sorted(groups):
                    kind, es, sz = k
                    tag = 'FBD' if kind == 'U2' else 'FBV'
                    out.append('*SOLID SECTION,ELSET=%s_%s,MATERIAL=%s'
                               % (tag, es, mat_of[es]))
        out.append(ln)

    open(a.dst, 'w').write('\n'.join(out) + '\n')

    print('fbares: %s -> %s   (c = %d)' % (a.src, a.dst, a.cycles))
    for es in sets:
        nel, nnd, ned, rs, ss = stats[es]
        E, nu = elastic.get(mat_of[es], (float('nan'), float('nan')))
        kg = (1.0 / (3.0 * (1.0 - 2.0 * nu))) / (1.0 / (2.0 * (1.0 + nu)))
        print('  %-14s %7d tets  %6d nodes  %7d edges   material %-14s '
              'K/G = %.0f' % (es, nel, nnd, ned, mat_of[es], kg))
        print('      U2 ring   mean %5.1f  p99 %4d  max %4d nodes '
              '(%d DOF)' % (rs.mean(), np.percentile(rs, 99), rs.max(),
                            3 * (rs.max() + 0)))
        print('      U3 stencil mean %5.1f  p99 %4d  max %4d nodes '
              '(%d DOF)' % (ss.mean(), np.percentile(ss, 99), ss.max(),
                            3 * (ss.max() + 0)))
    ndof = 3 * max(s[4].max() for s in stats.values())
    print('  %d element types; widest element %d DOF' % (len(groups), ndof))
    if ndof > 150:
        print('  NOTE: the e_c3d_u* family and mafillsm.f hold 150 DOF. '
              'Raise them to %d together, or ccx overruns s(150,150) '
              'silently -- patch 0006 documents this.' % (((ndof + 59) // 60) * 60))
    print('  run with CCX_U5_ZERO=1 CCX_FBAR_C=%d -- the base tets must '
          'contribute nothing' % a.cycles)


if __name__ == '__main__':
    main()
