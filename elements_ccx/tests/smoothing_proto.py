"""Numerical comparison of strain-smoothing schemes on linear tets.

Written to settle, before any Fortran, whether a FACE-based smoothing fixes
the instability that node-based smoothing (U5+U6) has on brine.

The cell is the failure in miniature: a brine slab normal to x, embedded in
ice, driven by a prescribed affine displacement u = eps.X on the whole
boundary.  Two numbers come out.

  C1111   2U/(eps^2 V), the confined stiffness.  This is the locking-sensitive
          modulus -- it contains K directly.  Too high = locking.
  fluc    max |u - eps.X| over the interior, divided by the applied
          displacement.  This is the mode detector used on the real cells:
          an affine boundary condition admits NO fluctuation in the exact
          solution of a homogeneous body, and in a two-phase body only a
          small one.  U5+U6 reads 19x on BRKB_b280.

Schemes:
  c3d4     standard displacement tet
  ns_vol   element deviatoric + NODE-smoothed volumetric  (this is U5+U6)
  ns_full  node-smoothed full strain                      (textbook NS-FEM)
  fs_full  face-smoothed full strain                      (textbook FS-FEM)
  fs_ns    FACE-smoothed deviatoric + node-smoothed volumetric
           (the selective FS/NS-FEM recommended for near-incompressible)
  fs_vol   face-smoothed volumetric + element deviatoric  (completeness)

Smoothing never crosses a material interface: a face whose two tets have
different materials is treated as a boundary face, and node patches are keyed
by (node, material), exactly as nodalbbar.py does it.
"""
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

M = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
T = np.diag([1.0, 1.0, 1.0, 0.5, 0.5, 0.5])


def dmat(K, G, part='full'):
    """Voigt [xx,yy,zz,gxy,gyz,gzx]; energy = 1/2 eps^T D eps."""
    vol = K * np.outer(M, M)
    dev = 2.0 * G * (T - np.outer(M, M) / 3.0)
    if part == 'vol':
        return vol
    if part == 'dev':
        return dev
    return vol + dev


def mesh_box(n, slab_lo, slab_hi, jitter=0.0, seed=7, geom=None):
    """n^3 hexes, 6 tets each (Freudenthal).  Returns nodes, tets, mat.

    `jitter` displaces interior nodes by up to jitter*h.  A perfectly
    structured mesh is the best case for a nodal-smoothing scheme; the real
    cells are unstructured Delaunay, and the spurious mode is known to be
    mesh-dependent, so the comparison has to be made on a disturbed mesh.
    Nodes on the periodic boundary are moved in matched pairs so that
    periodicity survives."""
    if geom is None:
        geom = GEOM['slab']
    g = np.linspace(0.0, 1.0, n + 1)
    X, Y, Z = np.meshgrid(g, g, g, indexing='ij')
    nodes = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    if jitter:
        rng = np.random.default_rng(seed)
        h = 1.0 / n
        # draw one displacement per PERIODIC IMAGE, so paired faces move alike
        key = np.round(np.mod(nodes * n, n)).astype(int)      # 0..n-1, wrapped
        uniq, inv = np.unique(key, axis=0, return_inverse=True)
        d = (rng.random((len(uniq), 3)) - 0.5) * 2.0 * jitter * h
        nodes = nodes + d[inv]

    def nid(i, j, k):
        return (i * (n + 1) + j) * (n + 1) + k

    # the 6-tet split of a cube, all with positive volume for this corner order
    split = ((0, 1, 3, 7), (0, 1, 7, 5), (0, 5, 7, 4),
             (0, 3, 2, 7), (0, 6, 4, 7), (0, 2, 6, 7))
    tets, mat = [], []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                c = [nid(i + (v & 1), j + ((v >> 1) & 1), k + ((v >> 2) & 1))
                     for v in range(8)]
                xc, yc, zc = (i + 0.5) / n, (j + 0.5) / n, (k + 0.5) / n
                m = geom(xc, yc, zc, slab_lo, slab_hi)      # 1 = brine
                for s in split:
                    tets.append([c[s[0]], c[s[1]], c[s[2]], c[s[3]]])
                    mat.append(m)
    return nodes, np.array(tets), np.array(mat)


# --------------------------------------------------------------------------
# GEOMETRY.  A brine slab spanning the cell is useless as a test: under a
# periodic eps_xx it is a one-dimensional layered problem, the strain is
# uniform inside each phase, so the nodal average of the divergence is EXACT
# and every scheme returns bit-identical answers (measured: 5.9656e+09 for
# c3d4, ns_vol and fs_ns alike).  The spurious mode needs brine nodes that
# touch no ice element at all -- a three-dimensional pocket.
def _slab(x, y, z, lo, hi):
    return 1 if lo <= x < hi else 0


def _sphere(x, y, z, lo, hi):
    r = hi                                   # hi carries the radius
    return 1 if (x - .5) ** 2 + (y - .5) ** 2 + (z - .5) ** 2 < r * r else 0


def _bridged(x, y, z, lo, hi):
    """Brine slab pierced by a square lattice of ice bridges -- the BRKB
    geometry in miniature, and the family where the mode was measured."""
    if not (lo <= x < hi):
        return 0
    b = 0.18                                 # bridge half-width
    if (abs(((y * 3) % 1.0) - .5) < b) and (abs(((z * 3) % 1.0) - .5) < b):
        return 0                             # ice bridge
    return 1


GEOM = {'slab': _slab, 'sphere': _sphere, 'bridged': _bridged}


def grads(nodes, tets):
    """Shape-function gradients and volumes for every tet."""
    p = nodes[tets]                                  # (ne,4,3)
    e = p[:, 1:, :] - p[:, :1, :]                    # (ne,3,3)
    det = np.linalg.det(e)
    vol = det / 6.0
    if np.any(vol <= 0):
        # flip the last two nodes of the inverted ones
        bad = vol <= 0
        tets[bad] = tets[bad][:, [0, 1, 3, 2]]
        p = nodes[tets]
        e = p[:, 1:, :] - p[:, :1, :]
        vol = np.linalg.det(e) / 6.0
    inv = np.linalg.inv(e)                           # (ne,3,3)
    g = np.zeros((len(tets), 4, 3))
    g[:, 1:, :] = np.transpose(inv, (0, 2, 1))
    g[:, 0, :] = -g[:, 1:, :].sum(axis=1)
    return g, vol


def bmat(g):
    """6x12 strain-displacement for one tet from its 4x3 gradients."""
    B = np.zeros((6, 12))
    for a in range(4):
        gx, gy, gz = g[a]
        c = 3 * a
        B[0, c + 0] = gx
        B[1, c + 1] = gy
        B[2, c + 2] = gz
        B[3, c + 0] = gy; B[3, c + 1] = gx
        B[4, c + 1] = gz; B[4, c + 2] = gy
        B[5, c + 0] = gz; B[5, c + 2] = gx
    return B


# --------------------------------------------------------------------------
# BUBBLE ENRICHMENT (bFS-FEM, Techscience CMES v127n2).
#
# u = u_lin + phi_b * u_b per element, phi_b = 256 L1 L2 L3 L4, 3 internal dof.
#
# Over a WHOLE element int grad(phi_b) dV = 0 (the bubble vanishes on every
# face), so element-averaged smoothing cannot see the bubble at all.  It is
# visible only to SUB-CELL averages, which is exactly why the S-FEM smoothing
# domains are built by splitting each tet at its centroid.  For the
# face-based domains the sub-cell is (3 face nodes + centroid), volume V/4.
#
# The bubble dofs are NOT condensed here.  With face smoothing an element's
# bubble couples to its four face-neighbours' bubbles, so the bubble block is
# not block-diagonal and element-local condensation is invalid.  Carrying the
# 3 dofs explicitly is exact and this is a prototype.
QP = np.array([[0.5854102, 0.1381966, 0.1381966, 0.1381966],
               [0.1381966, 0.5854102, 0.1381966, 0.1381966],
               [0.1381966, 0.1381966, 0.5854102, 0.1381966],
               [0.1381966, 0.1381966, 0.1381966, 0.5854102]])


def bub_pattern(gp):
    """6x3 strain matrix of the field phi*e_c given grad(phi) = gp."""
    B = np.zeros((6, 3))
    gx, gy, gz = gp
    B[0, 0] = gx
    B[1, 1] = gy
    B[2, 2] = gz
    B[3, 0] = gy; B[3, 1] = gx
    B[4, 1] = gz; B[4, 2] = gy
    B[5, 0] = gz; B[5, 2] = gx
    return B


def subcell_bubble_grads(nodes, tets, g):
    """Mean grad(phi_b) over each of a tet's 4 face sub-cells -> (ne,4,3)."""
    out = np.zeros((len(tets), 4, 3))
    for e, t in enumerate(tets):
        p = nodes[t]                       # (4,3) parent nodes
        cen = p.mean(axis=0)
        for drop in range(4):              # sub-cell opposite node `drop`
            verts = np.vstack([p[[i for i in range(4) if i != drop]], cen])
            acc = np.zeros(3)
            for w in QP:                   # 4-point degree-2 rule, weights 1/4
                x = w @ verts
                L = np.array([1.0 + g[e, k] @ (x - p[k]) for k in range(4)])
                gp = np.zeros(3)
                for k in range(4):
                    prod = np.prod([L[j] for j in range(4) if j != k])
                    gp += prod * g[e, k]
                acc += 0.25 * 256.0 * gp
            out[e, drop] = acc
    return out


def topology(tets, mat):
    """faces -> the tets sharing them; (node,material) -> its patch of tets."""
    faces = {}
    for e, t in enumerate(tets):
        for drop in range(4):
            f = tuple(sorted(int(t[i]) for i in range(4) if i != drop))
            faces.setdefault(f, []).append(e)
    patch = {}
    for e, t in enumerate(tets):
        for a in t:
            patch.setdefault((int(a), int(mat[e])), []).append(e)
    return faces, patch


class Assembler:
    def __init__(self, ndof):
        self.ndof = ndof
        self.r, self.c, self.v = [], [], []

    def add(self, dofs, Kloc):
        d = np.asarray(dofs)
        self.r.append(np.repeat(d, len(d)))
        self.c.append(np.tile(d, len(d)))
        self.v.append(Kloc.ravel())

    def tocsr(self):
        return sp.coo_matrix(
            (np.concatenate(self.v),
             (np.concatenate(self.r), np.concatenate(self.c))),
            shape=(self.ndof, self.ndof)).tocsr()


def dofs_of(nodelist):
    return np.array([3 * n + c for n in nodelist for c in range(3)])


def expand(B, tet, nodelist):
    """Re-index a tet's 6x12 B onto a wider node list."""
    idx = {n: i for i, n in enumerate(nodelist)}
    out = np.zeros((6, 3 * len(nodelist)))
    for a in range(4):
        j = 3 * idx[int(tet[a])]
        out[:, j:j + 3] = B[:, 3 * a:3 * a + 3]
    return out


def assemble(scheme, nodes, tets, mat, g, vol, faces, patch, props,
             stab=0.0, bubble=False, sbg=None):
    """props[m] = (K, G) for material m.  Bubble dofs follow the nodal ones."""
    nn = len(nodes)
    ndof = 3 * nn + (3 * len(tets) if bubble else 0)
    A = Assembler(ndof)
    B = [bmat(g[e]) for e in range(len(tets))]

    def bdofs(e):
        return np.array([3 * nn + 3 * e + c for c in range(3)])

    # ---- element-wise parts -------------------------------------------
    if scheme in ('c3d4', 'ns_vol', 'fs_vol'):
        for e, t in enumerate(tets):
            K, G = props[mat[e]]
            if scheme == 'c3d4':
                D = dmat(K, G, 'full')
            elif scheme == 'ns_vol':
                D = dmat(stab * K, G, 'full') if stab else dmat(K, G, 'dev')
            else:
                D = dmat(K, G, 'dev')
            A.add(dofs_of(t), vol[e] * B[e].T @ D @ B[e])

    # ---- face-smoothed parts ------------------------------------------
    if scheme in ('fs_full', 'fs_ns', 'fs_vol'):
        part = {'fs_full': 'full', 'fs_ns': 'dev', 'fs_vol': 'vol'}[scheme]
        for f, els in faces.items():
            # never smooth across a material interface
            groups = ([els] if len(els) == 2 and mat[els[0]] == mat[els[1]]
                      else [[e] for e in els])
            for grp in groups:
                nl = sorted(set(int(x) for e in grp for x in tets[e]))
                Vk = sum(vol[e] for e in grp) / 4.0
                cols = [sum((vol[e] / 4.0) * expand(B[e], tets[e], nl)
                            for e in grp) / Vk]
                dl = [dofs_of(nl)]
                if bubble:
                    for e in grp:
                        drop = [i for i in range(4)
                                if int(tets[e][i]) not in f][0]
                        cols.append((vol[e] / 4.0) * bub_pattern(sbg[e, drop])
                                    / Vk)
                        dl.append(bdofs(e))
                Bb = np.hstack(cols)
                K, G = props[mat[grp[0]]]
                A.add(np.concatenate(dl), Vk * Bb.T @ dmat(K, G, part) @ Bb)

    # ---- F-barES-FEM-T4 -------------------------------------------------
    if scheme.startswith('fbar'):
        cyc = int(scheme.split('_')[-1]) if scheme.split('_')[-1].isdigit() else 1
        for m in sorted(set(int(x) for x in mat)):
            sel, idx, ekeys, edges, Vh, Edev, Svol = fbar_es_operators(
                nodes, tets, mat, g, vol, cyc, m)
            K, G = props[m]
            for h, k in enumerate(ekeys):
                els_h = edges[k]
                nl = sorted(set(int(x) for e in els_h for x in tets[e]))
                # deviatoric: edge-smoothed once -- a narrow stencil, so a
                # dense block per edge is fine
                Bd = sum((vol[e] / 6.0 / Vh[h]) * expand(B[e], tets[e], nl)
                         for e in els_h)
                A.add(dofs_of(nl), Vh[h] * Bd.T @ dmat(K, G, 'dev') @ Bd)

            # VOLUMETRIC, as a sparse triple product.
            #
            # E A^c spans a wide element neighbourhood -- roughly 2c+1 rings --
            # so building one dense outer product per edge over that support
            # costs ~(3*300)^2 per edge and ran the machine out of memory.
            # The operator is  K_vol = K * Gv^T diag(V_h) Gv  with
            # Gv = (E A^c) D_div, and both factors are sparse, so the whole
            # thing is one sparse product instead of 13k dense blocks.
            rd, cd, vd = [], [], []
            for j, e in enumerate(sel):
                row = M @ B[e]                       # 1 x 12, the divergence
                for a in range(4):
                    for cc in range(3):
                        w = row[3 * a + cc]
                        if w:
                            rd.append(j)
                            cd.append(3 * int(tets[e][a]) + cc)
                            vd.append(w)
            Ddiv = sp.coo_matrix((vd, (rd, cd)),
                                 shape=(len(sel), ndof)).tocsr()
            Gv = (Svol @ Ddiv).tocsr()
            W = sp.diags(K * Vh)
            Kv = (Gv.T @ W @ Gv).tocoo()
            A.r.append(Kv.row); A.c.append(Kv.col); A.v.append(Kv.data)

    # ---- node-smoothed parts ------------------------------------------
    if scheme in ('ns_vol', 'ns_full', 'fs_ns'):
        part = 'full' if scheme == 'ns_full' else 'vol'
        for (a, m), els in patch.items():
            nl = sorted(set(int(x) for e in els for x in tets[e]))
            Va = sum(vol[e] for e in els) / 4.0
            K, G = props[m]
            if part == 'full':
                Bb = sum((vol[e] / 4.0) * expand(B[e], tets[e], nl)
                         for e in els) / Va
                A.add(dofs_of(nl), Va * Bb.T @ dmat(K, G, 'full') @ Bb)
                continue
            # volumetric only: one scalar row, theta_a
            row = [sum((vol[e] / 4.0) * (M @ expand(B[e], tets[e], nl))
                       for e in els) / Va]
            dl = [dofs_of(nl)]
            if bubble:
                # int_e N_a div(phi_b u_b) dV = -(32 V_e/105) g_a . u_b
                for e in els:
                    la = [i for i in range(4) if int(tets[e][i]) == a][0]
                    row.append(-(32.0 * vol[e] / 105.0) * g[e, la] / Va)
                    dl.append(bdofs(e))
            r = np.concatenate(row)
            A.add(np.concatenate(dl),
                  K * (1.0 - stab) * Va * np.outer(r, r))

    return A.tocsr()


def run(scheme, n, slab, props, eps=1e-3, stab=0.0, bubble=False,
        jitter=0.0, geomname='slab'):
    nodes, tets, mat = mesh_box(n, slab[0], slab[1], jitter,
                                geom=GEOM[geomname])
    g, vol = grads(nodes, tets)
    faces, patch = topology(tets, mat)
    sbg = subcell_bubble_grads(nodes, tets, g) if bubble else None
    K = assemble(scheme, nodes, tets, mat, g, vol, faces, patch, props,
                 stab, bubble, sbg)

    # PERIODIC boundary condition, u(x + L e_i) = u(x) + eps.(L e_i).
    #
    # This is what the real cells use, and it matters: with an affine
    # Dirichlet condition on the whole boundary the brine slab is clamped on
    # every side and the bellows mode simply cannot form -- every scheme then
    # reads a fluctuation of 0.2-0.4 and the test discriminates nothing.
    nn = len(nodes)
    ndof = K.shape[0]
    tol = 1e-9
    # master of each node: its image in [0,1)^3
    base = np.mod(np.round(nodes * 1e9) / 1e9, 1.0)
    base[np.abs(base - 1.0) < tol] = 0.0
    key = np.round(base * 1e7).astype(np.int64)
    uniq, inv = np.unique(key, axis=0, return_inverse=True)
    first = np.full(len(uniq), -1, dtype=int)
    for i in range(nn):
        if first[inv[i]] < 0:
            first[inv[i]] = i
    master = first[inv]
    shift = nodes - nodes[master]                    # 0 or +/-1 per component

    # u = T q + a ; q are the master nodal dofs plus every bubble dof
    nm = len(uniq)
    nq = 3 * nm + (ndof - 3 * nn)
    rows, cols, vals = [], [], []
    for i in range(nn):
        m = inv[i]
        for c in range(3):
            rows.append(3 * i + c); cols.append(3 * m + c); vals.append(1.0)
    for j in range(3 * nn, ndof):
        rows.append(j); cols.append(3 * nm + (j - 3 * nn)); vals.append(1.0)
    T = sp.coo_matrix((vals, (rows, cols)), shape=(ndof, nq)).tocsr()

    a = np.zeros(ndof)
    a[0:3 * nn:3] = eps * shift[:, 0]                # only eps_xx is applied

    Kr = (T.T @ K @ T).tocsr()
    fr = -(T.T @ (K @ a))
    # one master node pinned to kill the rigid translation
    keep = np.ones(nq, dtype=bool)
    keep[0:3] = False
    q = np.zeros(nq)
    q[keep] = spl.spsolve(Kr[keep][:, keep].tocsc(), fr[keep])
    u = T @ q + a

    energy = 0.5 * u @ (K @ u)
    c1111 = 2.0 * energy / (eps ** 2)                # unit volume
    # u is the TOTAL displacement (the jump `a` makes it so), defined up to a
    # constant.  The fluctuation is what is left after the macroscopic affine
    # field and that constant are removed -- the same quantity the real cells
    # report as max|u| against the applied displacement.
    un = u[:3 * nn]
    aff = np.zeros_like(un)
    aff[0::3] = eps * nodes[:, 0]
    fluct = un - aff
    fluct = fluct - fluct.mean()
    fl = np.abs(fluct).max() / eps

    # PRESSURE OSCILLATION, measured rather than assumed.
    #
    # p_e = K div(u)|_e is the element pressure the material actually feels,
    # whatever the scheme used to build the stiffness.  A stable element gives
    # neighbouring brine elements nearly the same pressure, so the face jump
    # is O(h) against the mean.  A checkerboard gives jumps of order the mean
    # or larger.  Reported as mean and max jump over brine-brine faces,
    # normalised by the mean |p| in the brine.
    soft = np.flatnonzero(mat == 1)
    pe = np.zeros(len(tets))
    for e in soft:
        ue = un[dofs_of(tets[e])]
        pe[e] = props[1][0] * (M @ (bmat(g[e]) @ ue))
    pbar = np.abs(pe[soft]).mean()
    jumps = []
    for f, els in faces.items():
        if len(els) == 2 and mat[els[0]] == 1 and mat[els[1]] == 1:
            jumps.append(abs(pe[els[0]] - pe[els[1]]))
    jumps = np.array(jumps) if jumps else np.zeros(1)
    osc = jumps.mean() / pbar if pbar else float('nan')
    oscmax = jumps.max() / pbar if pbar else float('nan')
    return c1111, fl, osc, oscmax, len(nodes), len(tets)


def runR(scheme, n, slab, Kb, Gb, ice, bubble=False, stab=0.0,
         jitter=0.0):
    """R = C1111(undrained brine) / C1111(drained brine), the campaign's own
    ratio, in miniature.  It needs no external reference: R -> 1 means the
    scheme has lost the brine's bulk stiffness entirely, which is exactly how
    MINI failed on the real cells (R 6.29 -> 1.89)."""
    und = {0: ice, 1: (Kb, Gb)}
    drn = {0: ice, 1: (Kb / 1000.0, Gb)}
    cu, fl, osc, oscm, _, _ = run(scheme, n, slab, und, stab=stab,
                                  bubble=bubble, jitter=jitter)
    # THE DENOMINATOR IS ALWAYS PLAIN C3D4, exactly as the campaign builds R:
    # Abaqus substitutes the hybrid element in the undrained cell only, so the
    # drained twin must stay the element both codes share.
    cd, _, _, _, _, _ = run('c3d4', n, slab, drn, jitter=jitter)
    return cu / cd, fl, osc, oscm


def main():
    Gb = 4.4e5
    Ei, ni = 9.37e9, 0.33
    ice = (Ei / (3 * (1 - 2 * ni)), Ei / (2 * (1 + ni)))

    n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    global JIT
    JIT = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
    lo, hi = 0.375, 0.625
    print('ice K/G = %.2f   slab x in [%.3f,%.3f]   n = %d   jitter = %.2f'
          % (ice[0] / ice[1], lo, hi, n, JIT))
    print('periodic cell, macroscopic eps_xx applied')

    exact = ice[0] + 4.0 * ice[1] / 3.0
    print('\n-- patch test, ice only (exact C1111 = %.6e) --' % exact)
    for sc in ['c3d4', 'ns_vol', 'ns_full', 'fs_full', 'fs_ns', 'fs_vol']:
        for bub in (False, True):
            if bub and sc not in ('fs_full', 'fs_ns'):
                continue
            c, fl, _o, _om, _, _ = run(sc, 4, (2.0, 3.0), {0: ice, 1: ice},
                                       bubble=bub)
            print('   %-9s%-8s C1111=%.6e  rel err=%9.2e'
                  % (sc, ' +bubble' if bub else '', c, abs(c / exact - 1)))

    arms = [('c3d4', False), ('ns_vol', False), ('fs_full', False),
            ('fs_ns', False), ('fs_full', True), ('fs_ns', True)]
    print('\n-- R = C1111(undrained)/C1111(drained) vs bulk-to-shear ratio --')
    print('   the paper validates bFS-FEM at K/G = 25..100; brine is 5000')
    hdr = '   %-16s' % 'K/G' + ''.join('%11s' % ('%.0f' % r)
                                       for r in (50, 100, 500, 1000, 5000))
    print(hdr)
    for sc, bub in arms:
        tag = sc + (' +bub' if bub else '')
        row, extra = '   %-16s' % tag, []
        for ratio in (50, 100, 500, 1000, 5000):
            R, fl, osc, oscm = runR(sc, n, (lo, hi), ratio * Gb, Gb, ice,
                                    bubble=bub, jitter=JIT)
            row += '%11.4f' % R
            extra.append((fl, osc))
        print(row)
        print('   %-16s' % '  fluc / p-osc'
              + ''.join('%11s' % ('%.1f/%.2f' % e) for e in extra))


if __name__ == '__main__':
    main()


# --------------------------------------------------------------------------
# F-barES-FEM-T4  (Onishi, Iida & Amaya, IJCM 15(7) 1845003, 2018)
#
# Small-strain reduction of eqs (1)-(11).  Two smoothing paths, deliberately
# different widths:
#
#   DEVIATORIC   edge-based, ONCE.  The smoothing domain of edge h is the set
#                of elements touching h, weighted V_e/6 (a T4 has 6 edges), so
#                sum_h V_h = sum_e V_e.                          [eqs 1-5]
#
#   VOLUMETRIC   J is smoothed c times through a node<->element cycle
#                  node:    J_n = sum_{e at n} J_e V_e/4 / sum V_e/4     (6)
#                  element: J_e = (1/4) sum_{n in e} J_n                 (7)
#                then edge-smoothed once                                 (8)
#                and recombined as F = J^(1/3) F~^iso                    (11)
#
# In small strain J - 1 = div u, so the cycle is a linear smoother on the
# element divergences and the whole thing is
#
#   eps_h = dev(eps~_h) + (1/3) dbar_h I,   dbar = E A^c D u
#
# with A = P Q the one-cycle element->node->element operator.  Energy is then
# G|dev eps~|^2 + (K/2) dbar^2 per edge domain.
#
# c = 0 recovers plain selective ES-FEM (edge deviatoric, edge volumetric).
# Smoothing never crosses a material interface, as everywhere else here.
def fbar_es_operators(nodes, tets, mat, g, vol, cycles, m):
    """Return (edge list, V_h, Bdev rows, dbar rows) for material m."""
    import scipy.sparse as _sp
    sel = np.flatnonzero(mat == m)
    idx = {e: i for i, e in enumerate(sel)}
    ne = len(sel)

    # element -> node (weights V_e/4) and node -> element (1/4)
    rq, cq, vq = [], [], []
    Vn = {}
    for e in sel:
        for a in tets[e]:
            Vn[int(a)] = Vn.get(int(a), 0.0) + vol[e] / 4.0
    nid = {a: i for i, a in enumerate(sorted(Vn))}
    for e in sel:
        for a in tets[e]:
            rq.append(nid[int(a)]); cq.append(idx[e])
            vq.append(vol[e] / 4.0 / Vn[int(a)])
    Q = _sp.coo_matrix((vq, (rq, cq)), shape=(len(nid), ne)).tocsr()
    rp, cp, vp = [], [], []
    for e in sel:
        for a in tets[e]:
            rp.append(idx[e]); cp.append(nid[int(a)]); vp.append(0.25)
    P = _sp.coo_matrix((vp, (rp, cp)), shape=(ne, len(nid))).tocsr()

    # edges of this material
    edges = {}
    for e in sel:
        t = [int(x) for x in tets[e]]
        for i in range(4):
            for j in range(i + 1, 4):
                edges.setdefault((min(t[i], t[j]), max(t[i], t[j])), []).append(e)
    ekeys = sorted(edges)
    Vh = np.array([sum(vol[e] / 6.0 for e in edges[k]) for k in ekeys])
    re, ce, ve = [], [], []
    for h, k in enumerate(ekeys):
        for e in edges[k]:
            re.append(h); ce.append(idx[e]); ve.append(vol[e] / 6.0 / Vh[h])
    E = _sp.coo_matrix((ve, (re, ce)), shape=(len(ekeys), ne)).tocsr()

    A = (P @ Q).tocsr()
    S = E
    for _ in range(cycles):
        S = (S @ A).tocsr()
    return sel, idx, ekeys, edges, Vh, E, S
