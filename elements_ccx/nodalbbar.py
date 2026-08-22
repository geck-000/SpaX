"""Nodal-averaged B-bar for one phase of a CalculiX deck.

    python3 nodalbbar.py <in-ccx.inp> <out.inp> [--elset Sphere_Only]

Splits the phase's response in two so the bulk modulus never appears in an
element-local operator and cannot lock:

    K = K_dev   per tet, element type U5 (deviatoric only)
      + K_vol   per node, built here out of ccx primitives

The nodal volumetric strain is
    V_a     = sum over the phase's elements at node a of V_e/4
    theta_a = (1/V_a) sum over those elements of (V_e/4) div(u)|_e
            = sum over the 1-ring of  b_n^a . u_n
and its energy is (1/2) K V_a theta_a^2.

Rather than build a patch element, theta is carried as one extra dof on a
dummy node, tied by an *EQUATION and given energy by a grounded SPRING1 --
all existing ccx features, so no change to the matrix-structure code.

The scaling is folded into the equation so every spring is identical:
with  t_a = sqrt(K V_a) theta_a  the energy is exactly (1/2) t_a^2, i.e. unit
stiffness for all of them, and one *SPRING card covers the lot. Emitting one
card per node would otherwise mean tens of thousands of element sets.

    *EQUATION :  t_a  -  sum_n sqrt(K V_a) b_n^a . u_n  =  0

t_a is written first because ccx eliminates the leading term as the dependent
dof, and t_a is the one that must go.
"""
import os
import sys
import numpy as np


def main():
    argv = sys.argv[1:]
    # --elset may be given more than once. Covering only the soft phase gives
    # its interface nodes ONE-SIDED patches, and in a slab 2-3 elements thick
    # most of its nodes are interface nodes -- the averaging then has almost
    # nothing to average over, which is why the brine-only version tracked
    # plain C3D4 under refinement while the same element gained 50x on a
    # homogeneous cantilever. Pass every elset to average across the interface.
    elsets = []
    while '--elset' in argv:
        i = argv.index('--elset')
        elsets.append(argv[i + 1])
        del argv[i:i + 2]
    if not elsets:
        elsets = ['Sphere_Only']
    elset = elsets[0]
    src, dst = argv[0], argv[1]
    lines = open(src).read().split('\n')

    def iskw(l):
        return l.startswith('*') and not l.startswith('**')

    # --- parse nodes, the phase's elements, and the material -------------
    co, els, maxnode = {}, [], 0
    mode = None
    cur_elset = None
    for ln in lines:
        if iskw(ln):
            u = ln.upper().replace(' ', '')
            if u.startswith('*NODE') and 'OUTPUT' not in u and 'PRINT' not in u \
               and 'FILE' not in u:
                mode = 'n'
            elif u.startswith('*ELEMENT') and any(
                    ('ELSET=' + e).upper() in u for e in elsets):
                mode = 'e'
                cur_elset = next(e for e in elsets
                                 if ('ELSET=' + e).upper() in u)
            else:
                mode = None
            continue
        f = [x.strip() for x in ln.split(',') if x.strip()]
        if mode == 'n' and len(f) >= 4:
            n = int(f[0])
            co[n] = np.array([float(f[1]), float(f[2]), float(f[3])])
            maxnode = max(maxnode, n)
        elif mode == 'e' and len(f) >= 5:
            els.append(([int(x) for x in f[1:5]], cur_elset))
    if not els:
        raise SystemExit('nodalbbar: no elements in ELSET=%s' % elset)

    # material of each elset and its bulk modulus. Each patch spans one
    # material, so each needs its own section card and its own K.
    mat_of = {}
    for ln in lines:
        u = ln.upper().replace(' ', '')
        if u.startswith('*SOLIDSECTION'):
            es = mt = None
            for fld in ln.split(','):
                t = fld.strip()
                if t.upper().startswith('ELSET='):
                    es = t.split('=', 1)[1].strip()
                elif t.upper().startswith('MATERIAL='):
                    mt = t.split('=', 1)[1].strip()
            if es and mt:
                mat_of[es] = mt
    missing = [e for e in elsets if e not in mat_of]
    if missing:
        raise SystemExit('nodalbbar: no *SOLID SECTION for %s' % missing)

    def elastic_of(mt):
        for i2, ln2 in enumerate(lines):
            if ln2.strip().upper().replace(' ', '').startswith(
                    '*MATERIAL,NAME=' + mt.upper()):
                for j2 in range(i2 + 1, min(i2 + 6, len(lines))):
                    f2 = [x.strip() for x in lines[j2].split(',') if x.strip()]
                    if len(f2) == 2 and not lines[j2].startswith('*'):
                        return float(f2[0]), float(f2[1])
                break
        return None, None

    Kof = {}
    for e in elsets:
        E, nu = elastic_of(mat_of[e])
        if E is None:
            raise SystemExit('nodalbbar: no elastic card for material %s'
                             % mat_of[e])
        Kof[e] = E / (3.0 * (1.0 - 2.0 * nu))

    # --- nodal volumes and the averaged divergence operator --------------
    Va = {}
    b = {}                       # b[a][n] = 3-vector
    for e, es in els:
        p = np.array([co[n] for n in e])
        J = np.array([p[1] - p[0], p[2] - p[0], p[3] - p[0]]).T
        det = np.linalg.det(J)
        vol = abs(det) / 6.0
        if vol <= 0.0:
            continue
        Ji = np.linalg.inv(J)
        gr = np.zeros((4, 3))
        gr[1:, :] = Ji                      # rows of J^-1 are grad L1..L3
        gr[0, :] = -gr[1:, :].sum(axis=0)   # grads sum to zero
        w = vol / 4.0
        for a in e:
            key = (a, es)
            Va[key] = Va.get(key, 0.0) + w
            ba = b.setdefault(key, {})
            for k, n in enumerate(e):
                ba[n] = ba.get(n, np.zeros(3)) + w * gr[k]
    for key in b:
        for n in b[key]:
            b[key][n] /= Va[key]

    # --- emit ----------------------------------------------------------
    #
    # The volumetric term goes in as U6 PATCH ELEMENTS, not as *EQUATION +
    # SPRING1. That first design was verified correct in every part -- and in
    # ccx it still came out 1.55x too stiff on a confined block, because
    # adjacent patches overlap so each mesh dof appears as an independent term
    # in ~24 equations and ccx's MPC machinery does not survive that. One
    # patch alone reproduced an independent Python assembly exactly; the full
    # set did not. As elements the MPCs are gone.
    #
    # *USER ELEMENT fixes NODES= per TYPE, and ring sizes vary, so one type is
    # declared per distinct ring size.
    keys = sorted(Va, key=lambda k: (k[1], k[0]))
    ring = {k: [k[0]] + sorted(n for n in b[k] if n != k[0]) for k in keys}
    bysize = {}
    for k in keys:
        bysize.setdefault((k[1], len(ring[k])), []).append(k)
    # One *USER ELEMENT type per (material, ring size): NODES= is fixed per
    # type. Two suffix characters give 36^2 names -- the label only has to
    # start 'U6' for the dispatcher, and positions 7-8 carry ndof and nope.
    alnum = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    names = [a + c for a in alnum for c in alnum]
    if len(bysize) > len(names):
        raise SystemExit('nodalbbar: %d (material, ring size) groups, more '
                         'than the U6 type names available' % len(bysize))
    suffix = {g: names[i] for i, g in enumerate(sorted(bysize))}

    out, done_step, u5decl = [], False, False
    eid = 10 ** 8
    for ln in lines:
        if iskw(ln):
            u = ln.upper().replace(' ', '')
            if u.startswith('*ELEMENT') and any(
                    ('ELSET=' + e).upper() in u for e in elsets):
                if not u5decl:
                    out.append('*USER ELEMENT,TYPE=U5,NODES=4,'
                               'INTEGRATIONPOINTS=1,MAXDOF=3')
                    u5decl = True
                out.append(ln.replace('C3D4', 'U5').replace('c3d4', 'U5'))
                continue
            if u.startswith('*STEP') and not done_step:
                done_step = True
                out.append('** SPAX nodal-averaged B-bar: %d patches, '
                           'one material each' % len(keys))
                # All *USER ELEMENT declarations first, then all *ELEMENT
                # blocks. ccx sorts the deck into per-keyword chains
                # (keystart.f: *USER ELEMENT is position 3, *ELEMENT is 4), and
                # interleaving them puts a multi-line element's continuation at
                # a chain boundary, where ccx reads it as a fresh element and
                # reports "element N is already defined".
                for grp in sorted(bysize):
                    es, sz = grp
                    out.append('*USER ELEMENT,TYPE=U6%s,NODES=%d,'
                               'INTEGRATIONPOINTS=1,MAXDOF=3'
                               % (suffix[grp], sz))
                for grp in sorted(bysize):
                    es, sz = grp
                    t = 'U6' + suffix[grp]
                    out.append('*ELEMENT,TYPE=%s,ELSET=BBARP_%s' % (t, es))
                    for k in bysize[grp]:
                        eid += 1
                        row = [str(eid)] + [str(x) for x in ring[k]]
                        # ccx takes at most 16 fields per line, and a continued
                        # card MUST end in a comma. Without it ccx reads the
                        # next line as a fresh element and reports
                        # "element N is already defined" -- the first field of
                        # the continuation is a node number that collides with
                        # a real element id. Rings here reach 37 nodes, so
                        # every patch over 15 was silently malformed.
                        # textpart in elements.f is dimensioned (16), so a
                        # line must not exceed 16 fields; a trailing comma adds
                        # one more and overflows it, after which ccx reads the
                        # continuation as a fresh element and reports
                        # "element N is already defined". Keep to 15 fields and
                        # no trailing comma -- ccx continues on node count.
                        for c in range(0, len(row), 15):
                            out.append(','.join(row[c:c + 15]))
                for es in elsets:
                    if any(g[0] == es for g in bysize):
                        out.append('*SOLID SECTION,ELSET=BBARP_%s,MATERIAL=%s'
                                   % (es, mat_of[es]))
        if iskw(ln) and ln.strip().upper().startswith('*STATIC'):
            # PARDISO by default. The assembled system is SPD, so incomplete-
            # Cholesky PCG is valid -- but the patch rows are ~5x denser than
            # C3D4's and carry a K/mu = 5000 contrast spread across each node
            # patch, which the preconditioner does not like: on LMESH_m0p0240
            # PARDISO finished while ITERATIVE CHOLESKY was still running after
            # 24 minutes. Override with SPAX_BBAR_SOLVER if memory is tighter
            # than time.
            ln = '*STATIC, SOLVER=' + os.environ.get('SPAX_BBAR_SOLVER',
                                                     'PARDISO')
        out.append(ln)

    open(dst, 'w').write('\n'.join(out))
    rsz = [len(ring[k]) for k in keys]
    print('nodalbbar: %s -> %s' % (src, dst))
    print('  %d elements retyped to U5 in ELSET=%s' % (len(els), elset))
    print('  %d U6 patches (one material each), ring %d..%d (mean %.1f), '
          '%d types' % (len(keys), min(rsz), max(rsz),
                        sum(rsz) / len(rsz), len(bysize)))
    for e in elsets:
        print('  %-14s material %-16s K = %.4e'
              % (e, mat_of[e], Kof[e]))
    print('  volumetric term delivered as elements -- no MPCs')


if __name__ == '__main__':
    main()
